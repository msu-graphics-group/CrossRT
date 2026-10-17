#include "FatBVH.h"
#include <iostream>
#include <cassert>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#define OSE_TESTS // enables internal debugging self-tests (usually as 'asserts()')

static const BVHNodeFat dummynode = FatBVH::BVHNodeFatDummy();


/// Element of array to sort in a DESCENDING order. Keeps old index.
struct FloatAndInt
  {
  float x;
  int i;

  FloatAndInt() { x = 0; i = -1; };

  FloatAndInt(float y, int j)
    {
    x = y;
    i = j;
    };

  /// ... and its comparator for sort()
  static bool greater(const FloatAndInt& a, const FloatAndInt& b)
    {
    return a.x > b.x;
    }
  };

/// The number of nodes in the subtree
/// @param[in] a_nodes       - array of nodes, may merge several BVH tree
/// @param[in] offset        - position of the root of the current whole tree in 'a_nodes'
/// @param[in] root_node_idx - index of the root node, or, to be precise, contents of
///                            offs_left / offs_right of the parent node. Counted from
///                            the root of the whole current tree
/// @return The number of nodes in the specified subtree
static int NumSubTreeNodes(const std::vector<BVHNodeFat>& a_nodes, int offset, int root_node_idx)
  {
  if (root_node_idx & LEAF_BIT)
    return 0; // it points not to a node but to a set of triangles

  // In fact EXTRACT_OFFSET() is not needed since (root_node_idx&LEAF_BIT)=0
  const BVHNodeFat& node = a_nodes[EXTRACT_OFFSET(root_node_idx) + offset];

  int num_nodes = 1; // root
  if (!(node.offs_left & LEAF_BIT)) // points to a node
    num_nodes += NumSubTreeNodes(a_nodes, offset, node.offs_left);
  if (!(node.offs_right & LEAF_BIT)) // points to a node
    num_nodes += NumSubTreeNodes(a_nodes, offset, node.offs_right);
  return num_nodes;
  }

/// The number of nodes in the subtree, bounded by the set of nodes
/// @param[in] a_nodes       - array of nodes, may merge several BVH tree
/// @param[in] offset        - position of the root of the current whole tree in 'a_nodes'
/// @param[in] root_node_idx - index of the root node, or, to be precise, contents of
///                            offs_left / offs_right of the parent node. Counted from
///                            the root of the whole current tree
/// @param[in] end           - the set of the (indices of) nodes that are "past the last 
///                            subtree node". That is we descend the tree until the 
///                            child node (index) is in the 'end'. Then we stop, NOT
///                            counting this node. 
///                            Ignored if empty.
/// @return The number of nodes in the subtree (from the root and to the 'end' nodes
///         or the terminal nodes, whichever is closer)
static int NumSubTreeNodes(const std::vector<BVHNodeFat>& a_nodes, int offset, int root_node_idx, const std::vector<int>& end)
  {
  if (root_node_idx & LEAF_BIT)
    return 0; // it points not to a node but to a set of triangles

  // In fact EXTRACT_OFFSET() is not needed since (root_node_idx&LEAF_BIT)=0
  const BVHNodeFat& node = a_nodes[EXTRACT_OFFSET(root_node_idx) + offset];

  int num_nodes = 1; // root
  if (!(node.offs_left & LEAF_BIT)) // points to a node
    {
    if (end.empty() || !std::binary_search(end.begin(), end.end(), node.offs_left)) // check presence
      num_nodes += NumSubTreeNodes(a_nodes, offset, node.offs_left, end);
    }
  if (!(node.offs_right & LEAF_BIT)) // points to a node
    {
    if (end.empty() || !std::binary_search(end.begin(), end.end(), node.offs_right)) // check presence
      num_nodes += NumSubTreeNodes(a_nodes, offset, node.offs_right, end);
    }
  return num_nodes;
  }

/// Enumerates nodes of a subtree, bounded by the set of nodes
/// @param[in] a_nodes       - array of nodes, may merge several BVH tree
/// @param[in] offset        - position of the root of the current whole tree in 'a_nodes'
/// @param[in] root_node_idx - index of the root node, or, to be precise, contents of
///                            offs_left / offs_right of the parent node. Counted from
///                            the root of the whole current tree
/// @param[in] end           - the set of the (indices of) nodes that are "past the last 
///                            subtree node". That is we descend the tree until the 
///                            child node (index) is in the 'end'. Then we stop, NOT
///                            counting this node
///                            HAS NO EFFECT IF EMPTY
/// @param[in,out] indices      - indices of ALL nodes of the sub-tree in the new order
///                               are APPENDED here.
///                               The indices are counted from the root of the whole tree
///                               i.e. they do not include the 'offset', i.e. the 
///                               node in the global array is a_nodes[indices[i]+offset]
/// @return 0 in case of success, -1 if failed
int CollectSubTreeNodes(const std::vector<BVHNodeFat>& a_nodes, int offset, int root_node_idx, const std::vector<int>& end, std::vector<int>& indices)
  {
  if (root_node_idx & LEAF_BIT)
    return 0; // it points not to a node but to a set of triangles

  // In fact EXTRACT_OFFSET() is not needed since (root_node_idx&LEAF_BIT)=0
  const BVHNodeFat& node = a_nodes[EXTRACT_OFFSET(root_node_idx) + offset];

  indices.emplace_back(root_node_idx); // root

  if (!(node.offs_left & LEAF_BIT)) // points to a node
    {
    if (end.empty() || !std::binary_search(end.begin(), end.end(), (int)node.offs_left)) // check presence
      CollectSubTreeNodes(a_nodes, offset, node.offs_left, end, indices);
    }
  if (!(node.offs_right & LEAF_BIT)) // points to a node
    {
    if (end.empty() || !std::binary_search(end.begin(), end.end(), (int)node.offs_right)) // check presence
      CollectSubTreeNodes(a_nodes, offset, node.offs_right, end, indices);
    }
  return 0;
  }

/// Suggests ordering of the nodes of a subtree "by clusters". 
/// 
/// New order of nodes is that they are stored "cluster by cluster", each cluster,
/// in turn, is usually further subdivided by sub-clusters etc. But if it is too small
/// i.e. contain <= 'cluster_size' nodes, it is not further subdivided. This 
/// deepest level of subdivision though can contain smaller number of nodes.
/// 
/// If a cluster is NOT further subdivided, its nodes are ordered so that their
/// criterion "w0" is descending. 
/// If it is further clusterized, then first comes its root cluster, then child 
/// clusters in such an order that first comes the one whose root has the largest 
/// "w0", then the one whose root has the next "w0" and so on. 
///
/// The order of nodes inside each sub-cluster (root or child) obeys the same
/// obey the same rule, i.e. if this sub-cluster is further subdivided, then... else
/// ordered so that their "w0" is descending.
/// 
/// @param[in] a_nodes       - array of nodes, may merge several BVH tree
/// @param[in] offset        - position of the root of the current whole tree in 'a_nodes'
/// @param[in] root_node_idx - index of the root of this subtree. Counted from
///                            the root of the whole current tree
/// @param[in] cluster_size  - max. allowed size of the "deepest level of subdivision" cluster
///                            If cluster is such or smaller it is NOT further clusterised.
///                            Best values are between 7 and 128, better 32. 
/// @param[in]     a_aligned - whether or not to pad the set of cluster nodes in memory
///                            by "dummy" elements to keep the size of the group
///                            exactly the same (as for the ideal tree) = cluster_size
/// @param[in] end           - (indices of) "end nodes" i.e. those BENEATH this subtree
///                            That is, if we descend from the subtree root and a
///                            child node is present in the 'end', then we do NOT
///                            go to this child, and here descending this branch stops.
///                            The indices are counted from the root of the whole tree
///                            i.e. they do not include the 'offset', i.e. the 
///                            node in the global array is a_nodes[end[i]+offset]
/// @param[in] num_nodes     - the number of nodes in the sub-tree, or -1 if unknown
/// @param[in] w0            - weights ("criterion") of all nodes in a_nodes. Can be empty, then approximate
///                            criterion is calculated approximately. 
/// @param[in,out] indices      - indices of ALL nodes of the sub-tree in the new order
///                               are APPENDED here.
///                            The indices are counted from the root of the whole tree
///                            i.e. they do not include the 'offset', i.e. the 
///                            node in the global array is a_nodes[indices[i]+offset]
/// @param[in,out] cluster_roots - here the (indices of) the roots of the LAST LEVEL (smallest, 
///                             those which are not further clusterized) clusters
///                             are APPENDED.
///                             These are indices in the global array @c m_allNodesFat
///                             (that may merge several BVH trees), i.e. do NOT add 'offset'.
/// 
/// @return SUCCESS/FAILURE
static int Clusterize(const std::vector<BVHNodeFat>& a_nodes,
                      int offset,
                      int root_node_idx, 
                      int cluster_size,
                      bool a_aligned,
                      const std::vector<int> &end,
                      int num_nodes,
                      const std::vector<double> &w0,
                      std::vector<int> &indices,
                      std::vector<int>& cluster_roots)
  {
  if (num_nodes < 0)
    num_nodes = NumSubTreeNodes(a_nodes, offset, root_node_idx, end);
#ifdef OSE_TESTS // tests
  else
    {
    if (num_nodes != NumSubTreeNodes(a_nodes, offset, root_node_idx, end))
      assert(false);
    }
#endif

  int iter;
  if (num_nodes <= cluster_size) // if already small, stop further clusterization and collect nodes to the array
    {
    cluster_roots.emplace_back(root_node_idx + offset);
    const size_t num_old = indices.size();
#ifdef OSE_TESTS // tests
    for (size_t i = 0; i < indices.size(); i++)
      {
      if (indices[i] == root_node_idx)
        assert(false);
      }
#endif

    if (!a_aligned)
      assert(num_old + num_nodes <= a_nodes.size());

    // Clusterization does not go further. Put all nodes within this cluster
    // to the target array in the "tree order".
    const int rc = CollectSubTreeNodes(a_nodes, offset, root_node_idx, end, indices);
    assert(indices.size() == num_old + num_nodes);

    if (a_aligned)
      {
      // Pad with dummy nodes, if necessarily, to have exactly 'cluster_size' nodes written to the array
      if (num_nodes < cluster_size)
        indices.insert(indices.end(), cluster_size - num_nodes, dummynode.offs_left); // 0xFFFFFFFF means dummy node
      }
    return rc;
    }

  // Target number of nodes in the NEXT LEVEL cluster (the real value can differ)
  // The root cluster occupies half of the tree levels, or HIGHER so that
  // the clusters "behind" it be SMALLER (because being deeper, they have
  // less regular structure and it looks better to use smaller clusters here)
  const int num_nodes_cluster = (2 << (((int)ceil(0.5 * std::log2(1 + (double)num_nodes))) - 1)) - 1;

  // Continue clusterization.

  // Find the "waistline" to split this cluster at the middle of tree depth.
  // Put nodes which are outside the cluster themselves while their parents are inside,
  // in "childroots" --- these will be the roots of the other clusters.
  std::vector<int> last, boundary;
  std::vector<FloatAndInt> w0_next, child_root_w0;

  last.emplace_back(root_node_idx);
  int num_nodes_real = (int)last.size(); 
  for (iter = 0; iter < 10000; iter++)
    {
    w0_next.clear();
    for (size_t i = 0; i < last.size(); i++)
      {
      const BVHNodeFat &node = a_nodes[offset + last[i]];
      float w0_;
      if ((node.offs_left & LEAF_BIT) == 0)
        {
        if (end.empty() || !std::binary_search(end.begin(), end.end(), (int)node.offs_left)) // check presence
          {
          // The value of criterion "w0" for this node
          if (!w0.empty())
            w0_ = (float)w0[offset + node.offs_left];
          else
            w0_ = GetChildBoxAreaLeft(node) + 1e-10;
          w0_next.emplace_back(FloatAndInt(w0_, node.offs_left));
          }
        else // Boundary knot, it must be added for "smaller" bottom boundary of THIS cluster
          boundary.emplace_back(node.offs_left);
        }
      if ((node.offs_right & LEAF_BIT) == 0)
        {
        if (end.empty() || !std::binary_search(end.begin(), end.end(), (int)node.offs_right)) // check presence
          {
          // The value of criterion "w0" for this node
          if (!w0.empty())
            w0_ = (float)w0[offset + node.offs_right];
          else
            w0_ = GetChildBoxAreaRight(node) + 1e-10;
          w0_next.emplace_back(FloatAndInt(w0_, node.offs_right));
          }
        else // Boundary knot, it must be added for "smaller" bottom boundary of THIS cluster
          boundary.emplace_back(node.offs_right);
        }
      }
    if (w0_next.empty())
      break; // no new nodes are to be added => end

    // Sort in DESCENDING order in "w0" value
    std::sort(w0_next.begin(), w0_next.end(), FloatAndInt::greater);

    double total = 0;
    for (size_t i = 0; i < w0_next.size(); i++)
      total += w0_next[i].x;

    last.clear();

    if (num_nodes_real >= num_nodes_cluster)
      {
      // None of these nodes will go into the cluster => all contain child roots
      child_root_w0.insert(child_root_w0.end(), w0_next.begin(), w0_next.end());
      break;
      }

    // Discard nodes with inessential contribution to the sum of w0 or those adding
    // which will exceed the target cluster size
    const size_t imax = num_nodes_cluster - num_nodes_real;
    double sum = 0;
    int count = -1;
    for (size_t i = 0; i < w0_next.size(); i++)
      {
      if (sum < 0.95 * total && i < imax)///
        {
        count = i + 1;
        sum += w0_next[i].x;
        }
      else // this node shall not be included into the cluster
        child_root_w0.emplace_back(w0_next[i]);
      }
    assert(count > 0 && count <= w0_next.size());
    w0_next.resize(count);

    num_nodes_real += (int)w0_next.size();
    // next scanline:
    last.clear();
    last.resize(w0_next.size());
    for (size_t i = 0; i < w0_next.size(); i++)
      last[i] = w0_next[i].i;

    if (last.empty())
      break; // no new nodes are to be added => end
    }

  // The line "below" the cluster is nothing but the child roots 'child_root_w0' 
  // PLUS those nodes that were in the subtree boundary (they were not included
  // in the child_root_w0).
  last.resize(child_root_w0.size() + boundary.size());
  for (size_t i = 0; i < boundary.size(); i++)
    last[i] = boundary[i];
  for (size_t i = 0; i < child_root_w0.size(); i++)
    last[boundary.size() + i] = child_root_w0[i].i;

  // Now 'last' contains the line "beneath" the root cluster. 
  // Sort it so that we can check fast if the given node (index)
  // is still inside the cluster or not.
  std::sort(last.begin(), last.end());

  if (Clusterize(a_nodes, offset, root_node_idx, cluster_size, a_aligned, last, num_nodes_real, w0, indices, cluster_roots) != 0)
    {
    assert(false);
    return -1;
    }

  // Sort in DESCENDING order in "w0" value
  std::sort(child_root_w0.begin(), child_root_w0.end(), FloatAndInt::greater);

  for (size_t i = 0; i < child_root_w0.size(); i++)
    {
    // child_root_w0 in fact contains additionally the nodes from the 'end' line
    // kept to use this line as the "border" in the previous Clusterize() call.
    // These may not be nodes of root clusters; skip them.
    if (end.empty() || !std::binary_search(end.begin(), end.end(), child_root_w0[i].i)) // check presence
      {
///      const int num_old = indices.size();
      const int num_nodes_ = NumSubTreeNodes(a_nodes, offset, child_root_w0[i].i, end);
      if (Clusterize(a_nodes, offset, child_root_w0[i].i, cluster_size, a_aligned, end, num_nodes_/*-1*/, w0, indices, cluster_roots) != 0)
        {
        assert(false);
        return -1;
        }
      }
    }
  return 0;
  }


/// Changes the order of the tree nodes in memory WHLE RETAINING THE TREE TOPOLOGY
/// and updates the links to the children kept in the nodes accordingly.
/// 
/// New order of nodes is that they are stored "cluster by cluster", each cluster,
/// in turn, is usually further subdivided by sub-clusters etc. But if it is too small
/// i.e. contain <= 'cluster_size' nodes, it is not further subdivided. This 
/// deepest level of subdivision though can contain smaller number of nodes.
/// 
/// If a cluster is NOT further subdivided, its nodes are ordered so that their
/// criterion "w0" is descending. 
/// If it is further clusterized, then first comes its root cluster, then child 
/// clusters in such an order that first comes the one whose root has the largest 
/// "w0", then the one whose root has the next "w0" and so on. 
///
/// The order of nodes inside each sub-cluster (root or child) obeys the same
/// obey the same rule, i.e. if this sub-cluster is further subdivided, then... else
/// ordered so that their "w0" is descending.  
/// 
/// @param[in,out] a_nodes      - array of nodes, may merge several BVH tree
/// @param[in,out] cluster_roots - here the (indices of) the roots of the LAST LEVEL (smallest, 
///                                those which are not further clusterized) clusters
///                                are APPENDED.
///                                These are indices in the global array @c m_allNodesFat
///                                (that may merge several BVH trees), i.e. do NOT add 'offset'.
/// @param[in] offset        - position in 'a_nodes' of the root of the tree to process 
/// @param[in] cluster_size  - max. allowed size of the "deepest level of subdivision" cluster
/// @param[in] a_aligned     - whether or not to pad the set of cluster nodes in memory
///                            by "dummy" elements to keep the size of the group
///                            exactly the same (as for the ideal tree) = cluster_size
/// @param[in] w0            - weights ("criterion") of all nodes in a_nodes. Can be empty, then approximate
///                            criterion is calculated approximately. 
/// @return SUCCESS/FAILURE
int FatBVH::ReorderByClusters(std::vector<BVHNodeFat>& a_nodes, std::vector<int>& cluster_roots, int offset, int cluster_size, bool a_aligned, const std::vector<double>& w0)
  {
  std::vector<int> old2new_indices, new2old_indices;
  std::vector<BVHNodeFat> new_nodes;

#if 0 ///  For test
  FILE* fd = fopen("new2old.indices", "w");
  FILE* fd1 = fopen("old2new.indices", "w");
#endif

  const int num_nodes = NumSubTreeNodes(a_nodes, offset, 0); 

  // New order of the distribution of tree nodes in memory. 
  // The current tree/nodes data itself is not YET changed
  new2old_indices.clear();
  const std::vector<int> end;
  if (Clusterize(a_nodes, offset, 0, cluster_size, a_aligned, end, num_nodes, w0, new2old_indices, cluster_roots) != 0)
    {
    assert(false);
    return -1;
    }

  const int num_nodes_new = new2old_indices.size();
  if (!a_aligned)
    assert(num_nodes_new == num_nodes);
  old2new_indices.resize(num_nodes);

  new_nodes.resize(num_nodes_new);

  // Notice that in new2old_indices old2new_indices the index of element 
  // and its value both relate to the SINGLE CURRENT TREE.

  // Array of nodes in new order:
  for (int i = 0; i < num_nodes_new; i++)
    {
    const int n2o = new2old_indices[i];
    if (n2o == int(dummynode.offs_left)) // 0xFFFFFFFF signals this is dummy node, added for alignment only
      new_nodes[i] = dummynode;
    else
      {
      if (n2o & LEAF_BIT)
        assert(false);

      old2new_indices[n2o] = i;
      new_nodes[i] = a_nodes[offset + n2o];
      }
    }

  const int offset_new = offset; /// TODO:

  // Update links to children so that they fit new array
  for (int i = 0; i < num_nodes_new; i++)
    {
    BVHNodeFat& node = new_nodes[offset_new + i];
    if (node.offs_left != dummynode.offs_left && !(node.offs_left & LEAF_BIT)) 
      node.offs_left = old2new_indices[node.offs_left]; // keeps index of node in OLD array

    if (node.offs_right != dummynode.offs_right && !(node.offs_right & LEAF_BIT)) 
      node.offs_right = old2new_indices[node.offs_right]; // keeps index of node in OLD array
    }
  a_nodes = new_nodes;

  std::vector<int> new_cluster_roots;
  new_cluster_roots.resize(cluster_roots.size());

//  if (a_aligned) 
    {
    // This method works for both aligned and not aligned clusters
    int curr_root_counter = 0, next_root_idx = 0;
    for (int i = 0; i < num_nodes_new; i++)
      {
      const int idx_old = new2old_indices[i];
      if (idx_old == int(dummynode.offs_left))
        continue;

      if (idx_old + offset == next_root_idx)
        {
        // new cluster began
        new_cluster_roots[curr_root_counter] = offset_new + i;
        next_root_idx = (curr_root_counter + 1) < int(cluster_roots.size()) ? int(cluster_roots[curr_root_counter + 1]) : int(cluster_roots.size() + 1);

        curr_root_counter++;
        }
      }
    }
//  else
//    {
//    for (int i = 0; i < cluster_roots.size(); i++)
//      new_cluster_roots[i] = old2new_indices[cluster_roots[i]];
//    }
  cluster_roots = new_cluster_roots;

#if 0 ///  For test
  if (fd != NULL && !w0.empty())
    {
    for (int i = 0; i < num_nodes; i++)
      fprintf(fd, "%15i %15i %15.5f\n", i, new2old_indices[i], w0[offset + i]);
    }
  if (fd1 != NULL && !w0.empty())
    {
    for (int i = 0; i < num_nodes; i++)
      fprintf(fd1, "%15i %15i %15.5f\n", i, old2new_indices[i], w0[offset + old2new_indices[i]]);
    }
#endif

#if 0 ///  For test
  if (fd != NULL)
    fclose(fd);
  if (fd1 != NULL)
    fclose(fd1);
#endif

  return 0;
  }
