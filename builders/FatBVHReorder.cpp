#include "FatBVH.h"
#include <queue>
#include <stack>
#include <iostream>
#include <unordered_set>
#include <cassert>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline float SurfaceAreaOfChildren(uint32_t a_nodeOffset, const std::vector<BVHNodeFat>& a_nodes)
{
  if(a_nodeOffset & LEAF_BIT)
    return 0.0f;

  const auto currNode = a_nodes[a_nodeOffset];

  float summ = 0.0f;
  if(!(currNode.offs_left & LEAF_BIT))
  {
    const auto leftChild = a_nodes[currNode.offs_left];
    summ += GetChildBoxAreaLeft(leftChild);
    summ += GetChildBoxAreaRight(leftChild);
  }

  if(!(currNode.offs_right & LEAF_BIT))
  {
    const auto rightChild = a_nodes[currNode.offs_right];
    summ += GetChildBoxAreaLeft(rightChild);
    summ += GetChildBoxAreaRight(rightChild);
  }

  return summ;
}


struct TraverseDFS
{
  TraverseDFS(const std::vector<BVHNodeFat>& a_nodes, std::vector<BVHNodeFat>& a_newNodes) : nodes(a_nodes), newNodes(a_newNodes){}
  
  const std::vector<BVHNodeFat>& nodes; 
  std::vector<BVHNodeFat>&       newNodes;

  uint32_t VisitNode(uint32_t a_nodeOffset)
  {
    const auto &currNode      = nodes[a_nodeOffset];
    const uint32_t newOffset = uint32_t(newNodes.size());
    newNodes.push_back(currNode);

    if(!(currNode.offs_left & LEAF_BIT))
      newNodes[newOffset].offs_left = VisitNode(currNode.offs_left);

    if(!(currNode.offs_right & LEAF_BIT))
      newNodes[newOffset].offs_right = VisitNode(currNode.offs_right);

    return newOffset;
  }

  uint32_t VisitNodeOrdered(uint32_t a_nodeOffset)
  {
    const auto currNode      = nodes[a_nodeOffset];
    const uint32_t newOffset = uint32_t(newNodes.size());
    newNodes.push_back(currNode);

    const float leftSA  = SurfaceAreaOfChildren(newNodes[newOffset].offs_left, nodes);
    const float rightSA = SurfaceAreaOfChildren(newNodes[newOffset].offs_right, nodes);
    
    if(leftSA >= rightSA)
    {
      if(!(currNode.offs_left & LEAF_BIT))
        newNodes[newOffset].offs_left = VisitNodeOrdered(currNode.offs_left);
  
      if(!(currNode.offs_right & LEAF_BIT))
        newNodes[newOffset].offs_right = VisitNodeOrdered(currNode.offs_right);
    }
    else
    {
      if(!(currNode.offs_right & LEAF_BIT))
        newNodes[newOffset].offs_right = VisitNodeOrdered(currNode.offs_right);

      if(!(currNode.offs_left & LEAF_BIT))
        newNodes[newOffset].offs_left = VisitNodeOrdered(currNode.offs_left);
    }

    return newOffset;
  }
};

void FatBVH::ReorderDFL(std::vector<BVHNodeFat>& a_nodes)
{
  std::vector<BVHNodeFat> newNodes;
  newNodes.reserve(a_nodes.size());
  newNodes.resize(0);

  TraverseDFS trav(a_nodes, newNodes);
  trav.VisitNode(0);
  a_nodes = newNodes;
}

void FatBVH::ReorderODFL(std::vector<BVHNodeFat>& a_nodes)
{
  std::vector<BVHNodeFat> newNodes;
  newNodes.reserve(a_nodes.size());
  newNodes.resize(0);

  TraverseDFS trav(a_nodes, newNodes);
  trav.VisitNodeOrdered(0);
  a_nodes = newNodes;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FatBVH::ReorderBFL(std::vector<BVHNodeFat>& a_nodes)
{
  std::vector<BVHNodeFat> newNodes;
  std::vector<uint32_t>   newOffsetByOld(a_nodes.size());
  newNodes.reserve(a_nodes.size());
  newNodes.resize(0);

  std::queue<uint32_t> treeQueue;
  treeQueue.push(0);

  while (treeQueue.size() != 0) 
  {
    const uint32_t nodeOffset  = treeQueue.front();
    const auto&    currNode    = a_nodes[nodeOffset];
    newOffsetByOld[nodeOffset] = uint32_t(newNodes.size());
    newNodes.push_back(currNode);
    
    if(!(currNode.offs_left & LEAF_BIT))
      treeQueue.push(currNode.offs_left);

    if(!(currNode.offs_right & LEAF_BIT))
      treeQueue.push(currNode.offs_right);

    treeQueue.pop();
  }

  for(size_t i=0;i<newNodes.size();i++)
  {
    auto& currNode = newNodes[i];
    if(!(currNode.offs_left & LEAF_BIT))
      currNode.offs_left = newOffsetByOld[currNode.offs_left];
    if(!(currNode.offs_right & LEAF_BIT))
      currNode.offs_right = newOffsetByOld[currNode.offs_right];
  }
  a_nodes = newNodes;
}

std::vector<uint2> FatBVH::ComputeDepthRanges(const std::vector<BVHNodeFat>& nodes)
{
  std::vector<uint2>    depthRanges;
  std::queue<uint32_t>  treeQueue;
  std::vector<uint32_t> nodesOfCurrentLevel; 
  nodesOfCurrentLevel.reserve(1000);
  
  treeQueue.push(0);
  
  while (treeQueue.size() != 0)
  {
    while (treeQueue.size() != 0)
    {
      nodesOfCurrentLevel.push_back(treeQueue.front());
      treeQueue.pop();
    } 
    
    uint32_t start = 0xFFFFFFFF;
    uint32_t end   = 0;
    for(size_t j=0;j<nodesOfCurrentLevel.size();j++)
    {
      start = std::min(start, nodesOfCurrentLevel[j]);
      end   = std::max(start, nodesOfCurrentLevel[j]);
    }

    depthRanges.push_back(uint2(start, end-start+1));

    for(auto x : nodesOfCurrentLevel)
    { 
      const BVHNodeFat& node = nodes[x];
      if((node.offs_left & LEAF_BIT) == 0)
        treeQueue.push(node.offs_left);
      if((node.offs_right & LEAF_BIT) == 0)
        treeQueue.push(node.offs_right);
    }
    nodesOfCurrentLevel.clear();
  }
  
  return depthRanges;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BVHNodeFat FatBVH::BVHNodeFatDummy()
{
  BVHNodeFat node;
  node.lmin_xyz_rmax_x = LiteMath::float4(0.0f);
  node.lmax_xyz_rmax_y = LiteMath::float4(0.0f);
  node.rmin_xyz_rmax_z = LiteMath::float4(0.0f);
  // To mark the dummy nodes added merely to keep the treelet size constant
  node.offs_left  = 0xFFFFFFFF;///PackOffsetAndSize(0, 0);
  node.offs_right = 0xFFFFFFFF;///PackOffsetAndSize(0,0);
  return node;
}

/// Selects nodes that form treelet of desired size and "starting" from the desired
/// node (treelet's root)
/// 
/// @param[in]     a_nodes - array of nodes of the BVH tree
/// @param[in]     this_treelet_root - index of the root of the then treelet (in a_nodes)
/// @param[in]     a_leavesNum       - the number of the "bottom" (last along the tree)
///                                    treelet nodes. The number of all nodes is then 2*a_leavesNum-1
///                                    The actual number of nodes can be LESS! (if some inside
///                                    treelet are terminal nodes).
/// @param[out]    trVector - indices (related to 'a_nodes', no "offset" is used!!) of the
///                           nodes of treelet, first coming the root.
/// @param[out]    trChild  - APPEDS TO IT indices of (related to 'a_nodes', no "offset" is used!!) 
///                           those children of treelet's nodes which are OUTSIDE treelet. 
/// @param[in,out] trQueue  - workspace, no operation is needed from the caller
/// @param[in,out] wkIdx    - workspace. Must be allocated for a_nodes.size() elements,
///                           and all be set to @c false before the VERY FIRST call to
///                           this method. After that, no action must be needed from the
///                           caller, because on exit from this function array is restored.
/// @return SUCCESS/FAILURE
static int CollectOneTreeletNodes(const std::vector<BVHNodeFat>& a_nodes,
                                   size_t this_treelet_root,
                                   size_t a_leavesNum, 
                                   std::vector<uint32_t>& trVector,
                                   std::vector<uint32_t> &trChild,
                                   std::queue<uint32_t>  &trQueue,
                                   std::vector<bool>& wkIdx)
  {
  if (wkIdx.size() != a_nodes.size())
    {
    assert(false);
    return -1;
    }
  const size_t treeletSize = 2 * a_leavesNum - 1; // number of all nodes in treelet (both inner nodes and leaves)

  trVector.reserve(treeletSize); // +1 is for root node inside treelet
  trVector.resize(0);
  trChild.reserve(2 * treeletSize); // the factor 2 is for safety

  trQueue.push(this_treelet_root);
  while (!trQueue.empty() && trQueue.size() + trVector.size() < treeletSize)
    {
    const BVHNodeFat& currNode = a_nodes[trQueue.front()];
    if (!(currNode.offs_left & LEAF_BIT) && trQueue.size() + trVector.size() < treeletSize)
      trQueue.push(currNode.offs_left);

    if (!(currNode.offs_right & LEAF_BIT) && trQueue.size() + trVector.size() < treeletSize)
      trQueue.push(currNode.offs_right);

    trVector.push_back(trQueue.front());
    trQueue.pop();
    }

  while (!trQueue.empty())
    {
    trVector.push_back(trQueue.front());
    trQueue.pop();
    }

  assert(trVector.size() <= treeletSize);

  // Set wkIdx[idx]=true if 'idx' is index of a treelet element
  for (size_t i = 0; i < trVector.size(); i++)
    {
    const uint32_t& idx = trVector[i];
    if (wkIdx[idx])
      {
      assert(false); // not proper state
      return -1;
      }
    wkIdx[idx] = true;
    }

  // find those children which are outside treelet
  for (size_t i = 0; i < trVector.size(); i++)
    {
    const BVHNodeFat& currNode = a_nodes[trVector[i]];
    if (!(currNode.offs_left & LEAF_BIT))
      {
      if (!wkIdx[currNode.offs_left]) // currNode.offs_left is outside treelet
        trChild.emplace_back(currNode.offs_left);
      }

    if (!(currNode.offs_right & LEAF_BIT))
      {
      if (!wkIdx[currNode.offs_right]) // currNode.offs_left is outside treelet
        trChild.emplace_back(currNode.offs_right);
      }
    }

  // Return wkIdx to its initial state
  for (size_t i = 0; i < trVector.size(); i++)
    wkIdx[trVector[i]] = false;

  return 0;
  }

/// Collects all nodes of a super-treelet into a memory array
/// 
/// The result is array of arrays, whose element (inner array) contains
/// nodes of one treelet. The first element (of the outer array) is the root treelet.
/// 
/// @param[in,out]  a_nodes - array of nodes of the BVH tree. On input, in the original
///                           order. On EXIT, grouped by [[super]super]treelets.
/// @param[in]     a_leavesNum - the number of the "bottom" (last along the tree)
///                              treelet nodes. The number of all nodes is then 2*a_leavesNum-1
///                              The actual number of nodes can be LESS! (if some inside
///                              treelet are terminal nodes).
/// @param[out]    trVectors   - trVectors[i] conveys the indices (related to 'a_nodes', 
///                              no "offset" is used!!) of the i-th treelet ([0] is the
///                              root treelet, then come the child treelets). Each
///                              treelet begins with it root, then the child nodes in BFS.
/// @param[out]    trGrandChild - here APPENDED are the (indices of) such nodes OUSIDE
///                               the supertreelet but whose parens belong to the 
///                               super-treelet. 
///                               These indices are for the OLD order (i.e. relate to a_nodes). 
/// @param[in,out] trChild  - workspace, no operation is needed from the caller
/// @param[in,out] trQueue  - workspace, no operation is needed from the caller
/// @param[in,out] trVector - workspace, no operation is needed from the caller
/// @param[in,out] wkIdx    - workspace. Must be allocated for a_nodes.size() elements,
///                           and all be set to @c false before the VERY FIRST call to
///                           this method. After that, no action must be needed from the
///                           caller, because on exit from this function array is restored.
static void CollectOneSuperTreeletNodes(const std::vector<BVHNodeFat>& a_nodes,
                                        size_t this_treelet_root,
                                        size_t a_leavesNum, 
                                        std::vector< std::vector<uint32_t> > &trVectors,
                                        std::vector<uint32_t> &trGrandChild,
                                        std::vector<uint32_t> &trChild,
                                        std::queue<uint32_t>  &trQueue,
                                        std::vector<uint32_t> &trVector,
                                        std::vector<bool>& wkIdx)
  {
  //const size_t treeletSize = 2 * a_leavesNum - 1; // number of all nodes in treelet (both inner nodes and leaves)
  
  trVectors.resize(2 * a_leavesNum + 1); // 2*a_leavesNum child treelets + 1 root treelet
  trVectors.clear();


  trChild.clear();
  CollectOneTreeletNodes(a_nodes, this_treelet_root, a_leavesNum, trVector, trChild, trQueue, wkIdx);
  trVectors.emplace_back(trVector);

  // process the child treelets of the current and place them next to it in memory
  for (size_t j = 0; j < trChild.size(); j++)
    {
    CollectOneTreeletNodes(a_nodes, trChild[j], a_leavesNum, trVector, trGrandChild, trQueue, wkIdx);
    trVectors.emplace_back(trVector);
    }
  }


/// Appends a set of nodes to newNodes "by horizontal tree sections".
/// 
/// The set can be a single treelet, several treelets or an arbitrary collection
/// 
/// @param[in]     a_nodes - array of nodes of the BVH tree
/// @param[in]     trVector - array of positions of nodes in 'a_nodes'
/// @param[in,out] newNodes - accumulated array of tree nodes in the NEW ORDER.
///                           The nodes of the treelet created are APPENDED to it,
///                           first comes the treelet's root node.
/// @param[in,out] newOffsetByOld - lookup from the index of a node in old order and
///                                 the same node in the new order, i.e.
///                                 newNodes[newOffsetByOld[i]] = a_nodes[i].
///                                 
///                                 Before the very first call to ProcessOneTreelet() 
///                                 this array must be allocated (by a_nodes.size()) by the caller!!! 
///                                 And filled with 0xFFFFFFFF.
/// 
///                                 The data for nodes of the current set is written
///                                 into the related elements, the rest remain unchanged.
static void AppendNodes(const std::vector<BVHNodeFat>& a_nodes,
                        const std::vector<uint32_t>& trVector,
                        std::vector<BVHNodeFat>& newNodes,
                        std::vector<uint32_t>& newOffsetByOld)
  {
  const BVHNodeFat dummyNode = FatBVH::BVHNodeFatDummy();

  for (size_t j = 0; j < trVector.size();j++)
    {
    const uint32_t nodeOffset = trVector[j];
    if (nodeOffset == uint32_t(-1))
      newNodes.push_back(dummyNode); // dummyNode is leaf with zero primitives, so we don't need to change its offsets
    else
      {
      newOffsetByOld[nodeOffset] = uint32_t(newNodes.size());
      newNodes.push_back(a_nodes[nodeOffset]);
      }
    }
  }


/// Selects nodes that form treelet of desired size and "starting" from the desired
/// node (treelet's root)
/// 
/// That is, writes treelet nodes (starting from the root) to 'newNodes'
/// "by horizontal tree sections".
/// 
/// @param[in]     a_nodes - array of nodes of the BVH tree
/// @param[in]     this_treelet_root - index of the root of the then treelet (in a_nodes)
/// @param[in]     a_leavesNum       - the number of the "bottom" (last along the tree)
///                                    treelet nodes. The number of all nodes is then 2*a_leavesNum-1
///                                    The actual number of nodes can be LESS! (if some inside
///                                    treelet are terminal nodes).
/// @param[in]     a_aligned         - whether or not to pad the set of treelet nodes in memory
///                                    by "dummy" elements to keep the size of the group
///                                    exactly the same (as for the ideal tree)
/// @param[in,out] root_idx - accumulated array of the (indices of) treelet roots.
///                           The indices relate to the NEW ORDER of nodes, i.e.
///                           the root of the i-th treelet is newNodes[root_idx[i]].
///                           The root of the created treelet is APPENDED to this array.
/// @param[in,out] newNodes - accumulated array of tree nodes in the NEW ORDER.
///                           The nodes of the treelet created are APPENDED to it,
///                           first comes the treelet's root node.
/// @param[in,out] newOffsetByOld - lookup from the index of a node in old order and
///                                 the same node in the new order, i.e.
///                                 newNodes[newOffsetByOld[i]] = a_nodes[i].
///                                 
///                                 Before the very first call to ProcessOneTreelet() 
///                                 this array must be allocated (by a_nodes.size()) by the caller!!! 
///                                 And filled with 0xFFFFFFFF.
/// 
///                                 The daya for nodes of the current treelet is written
///                                 into the related elements, the rest remain unchanged.
/// @param[out]    trChild - (indices of) those children of treelet's nodes which are
///                          OUTSIDE treelet are APPENDED to this array. 
///                          These indices are for the OLD order (i.e. relate to a_nodes). 
/// @param[in,out] trQueue  - workspace, no operation is needed from the caller
/// @param[in,out] trVector - workspace, no operation is needed from the caller
static void ProcessOneTreelet(const std::vector<BVHNodeFat>& a_nodes,
                              size_t this_treelet_root,
                              size_t a_leavesNum, 
                              bool a_aligned,
                              std::vector<int>& root_idx, 
                              std::vector<BVHNodeFat> &newNodes,
                              std::vector<uint32_t> &newOffsetByOld,
                              std::vector<uint32_t> &trChild,
                              std::queue<uint32_t>  &trQueue,
                              std::vector<uint32_t> &trVector)
  {
  const size_t treeletSize = 2 * a_leavesNum - 1; // number of all nodes in treelet (both inner nodes and leaves)
  const size_t treeletSizeAligned = treeletSize;///+1;

  const BVHNodeFat dummyNode = FatBVH::BVHNodeFatDummy();

  trVector.reserve(treeletSize); // +1 is for root node inside treelet
  trChild.reserve(2 * treeletSize); // the factor 2 is for safety

  trQueue.push(this_treelet_root);
  trVector.resize(0);
  while (!trQueue.empty() && trQueue.size() + trVector.size() < treeletSize)
    {
    const auto& currNode = a_nodes[trQueue.front()];
    if (!(currNode.offs_left & LEAF_BIT) && trVector.size() + trQueue.size() < treeletSizeAligned)
      trQueue.push(currNode.offs_left);

    if (!(currNode.offs_right & LEAF_BIT) && trVector.size() + trQueue.size() < treeletSizeAligned)
      trQueue.push(currNode.offs_right);

    trVector.push_back(trQueue.front());
    trQueue.pop();
    }

  while (!trQueue.empty())
    {
    trVector.push_back(trQueue.front());
    trQueue.pop();
    }

  if (a_aligned && trVector.size() < treeletSizeAligned)
    {
    while (trVector.size() < treeletSizeAligned)
      trVector.push_back(uint32_t(-1));
    }

  root_idx.emplace_back(newNodes.size()); // position where the root will just be added

  //treelets.push_back(uint2(uint32_t(newNodes.size()), uint32_t(trVector.size()))); // save treelets info if we need
  for (size_t j = 0;j < trVector.size();j++)
    {
    const uint32_t nodeOffset = trVector[j];
    if (nodeOffset == uint32_t(-1))
      newNodes.push_back(dummyNode); // dummyNode is leaf with zero primitives, so we don't need to change its offsets
    else
      {
      newOffsetByOld[nodeOffset] = uint32_t(newNodes.size());
      newNodes.push_back(a_nodes[nodeOffset]);
      }
    }

  for (size_t j = 0;j < trVector.size();j++)
    {
    const uint32_t nodeOffset = trVector[j];
    if (nodeOffset != uint32_t(-1))
      {
      const auto& currNode = a_nodes[nodeOffset];

      if (!(currNode.offs_left & LEAF_BIT))
        {
        if (newOffsetByOld[currNode.offs_left] == 0xFFFFFFFF) // left child was not processed yet
          trChild.emplace_back(currNode.offs_left);
        }

      if (!(currNode.offs_right & LEAF_BIT))
        {
        if (newOffsetByOld[currNode.offs_right] == 0xFFFFFFFF) // right child was not processed yet
          trChild.emplace_back(currNode.offs_right);
        }
      }
    }
  }


/// Writes nodes of "supertreelet" (treelet of treelets) into memory.
/// 
/// That is, writes treelet nodes (starting from the root) to 'newNodes'
/// "by horizontal tree sections". Then writes to it all "child treelets"
/// that is, treelets whose roots are OUTSIDE the "root treelet" while parent
/// of their roots are INSIDE the "root treelet".
/// 
/// @param[in]     a_nodes - array of nodes of the BVH tree
/// @param[in]     this_treelet_root - index of the root of the then treelet (in a_nodes)
/// @param[in]     a_leavesNum       - the number of the "bottom" (last along the tree)
///                                    treelet nodes. The number of all nodes is then 2*a_leavesNum-1
///                                    The actual number of nodes can be LESS! (if some inside
///                                    treelet are terminal nodes).
/// @param[in]     a_aligned         - whether or not to pad the set of treelet nodes in memory
///                                    by "dummy" elements to keep the size of the group
///                                    exactly the same (as for the ideal tree)
/// @param[in,out] root_idx - accumulated array of the (indices of) treelet roots.
///                           The indices relate to the NEW ORDER of nodes, i.e.
///                           the root of the i-th treelet is newNodes[root_idx[i]].
///                           The root of the created treelet is APPENDED to this array.
/// @param[in,out] newNodes - accumulated array of tree nodes in the NEW ORDER.
///                           The nodes of the treelet created are APPENDED to it,
///                           first comes the treelet's root node.
/// @param[in,out] newOffsetByOld - lookup from the index of a node in old order and
///                                 the same node in the new order, i.e.
///                                 newNodes[newOffsetByOld[i]] = a_nodes[i].
///                                 
///                                 Before the very first call to ProcessOneTreelet() 
///                                 this array must be allocated (by a_nodes.size()) by the caller!!! 
///                                 And filled with 0xFFFFFFFF.
/// 
///                                 The data for nodes of the current treelet is written
///                                 into the related elements, the rest remain unchanged.
/// @param[out]    trGrandChild - here APPENDED are the (indices of) such nodes OUSIDE
///                               the supertreelet but whose parens belong to the 
///                               super-treelet. 
///                               These indices are for the OLD order (i.e. relate to a_nodes). 
/// @param[in,out] trChild  - workspace, no operation is needed from the caller
/// @param[in,out] trQueue  - workspace, no operation is needed from the caller
/// @param[in,out] trVector - workspace, no operation is needed from the caller
static void ProcessOneSuperTreelet(const std::vector<BVHNodeFat>& a_nodes,
                                   size_t this_treelet_root,
                                   size_t a_leavesNum,
                                   bool a_aligned,
                                   std::vector<int>& root_idx,
                                   std::vector<BVHNodeFat>& newNodes,
                                   std::vector<uint32_t>& newOffsetByOld,
                                   std::vector<uint32_t>& trGrandChild,
                                   std::vector<uint32_t>& trChild,
                                   std::queue<uint32_t>& trQueue,
                                   std::vector<uint32_t>& trVector)
  {
  trChild.clear();
  ProcessOneTreelet(a_nodes, this_treelet_root, a_leavesNum, a_aligned,
    root_idx, newNodes, newOffsetByOld, trChild, trQueue, trVector);

  // process the child treelets of the current and place them next to it in memory
  for (size_t j = 0; j < trChild.size();j++)
    ProcessOneTreelet(a_nodes, trChild[j], a_leavesNum, a_aligned,
                      root_idx, newNodes, newOffsetByOld, trGrandChild, trQueue, trVector);
  }

/// Reorders the nodes of the BVH tree in memory so that they are grouped.
/// 
/// If the 'level' is TREELETS then the group is all nodes of this treelet. Then
/// next to it (in memory) comes the next treelet (which is to the right of current
/// in the tree graph and so on).
/// 
/// If the 'level' is SUPERTREELETS then the group is all nodes of this supertreelet.
/// Supertreelet = treelet and its child treelets (those whose roots are outside the
/// 1st treelet but their parents are INSIDE it). First come the elements of the
/// "root treelet", then its child treelets (along the horizontal tree slices)
/// The rest is the same as for TREELETS.
/// 
/// If the 'level' is SUPERSUPERTREELETS then the group is all nodes of this supersupertreelet.
/// Supersupertreelet = supertreelet and its child supertreelets (those whose roots are outside the
/// 1st treelet but their parents are INSIDE it). First come the elements of the
/// "root supertreelet", then its child supertreelets (along the horizontal tree slices).
/// Each supertreelet (root or child) is processed as for the SUPETREELETS mode.
/// The rest is the same as for TREELETS.
/// 
/// @param[in,out]  a_nodes - array of nodes of the BVH tree. On input, in the original
///                           order. On EXIT, grouped by [[super]super]treelets.
/// @param[out]    root_idx - indices of treelet roots, related to the re-arranged array (new positions)
/// @param[out]    super_root_idx - indices of supertreelet roots, related to the re-arranged array (new positions)
/// @param[in]     a_leavesNum - the number of the "bottom" (last along the tree)
///                              treelet nodes. The number of all nodes is then 2*a_leavesNum-1
///                              The actual number of nodes can be LESS! (if some inside
///                              treelet are terminal nodes).
/// @param[in]     a_aligned   - whether or not to pad the set of treelet nodes in memory
///                              by "dummy" elements to keep the size of the group
///                              exactly the same (as for the ideal tree)
/// @param[in]  level -          the level of treelet-ization.
void FatBVH::ReorderTRB(std::vector<BVHNodeFat>& a_nodes, std::vector<int>& root_idx, std::vector<int>& super_root_idx, size_t a_leavesNum, bool a_aligned, TreeletizLevel level)
{
  const size_t treeletSize = 2*a_leavesNum-1; // number of all nodes in treelet (both inner nodes and leaves)

  //const BVHNodeFat dummyNode = BVHNodeFatDummy();

  //std::vector<uint2> treelets;
  //treelets.reserve(a_nodes.size()/treeletSize);

  std::vector<BVHNodeFat> newNodes;
  newNodes.reserve(a_nodes.size());
  newNodes.resize(0);

  std::vector<uint32_t> newOffsetByOld(a_nodes.size());
  for(auto& index : newOffsetByOld) // mark nodes are not processed yet
    index = 0xFFFFFFFF;

  std::queue<uint32_t>  trQueue;
  std::vector<uint32_t> trVector;
  trVector.reserve(treeletSize); // +1 is for root node inside treelet

  std::queue<uint32_t> roots;
  roots.push(0);

  std::vector<uint32_t> trChild, trChild2, trGrandChild;


  // nodes of supertreelet, first the root treelet then (maybe regrouped) child treelets
  std::vector< std::vector<uint32_t> > strVectors, strVectorsNew;

  //static int num0 = 0;
  while (roots.size() != 0)
  {
    trChild.clear();
    if (level == TreeletizLevel::TREELETS)
      {
      ProcessOneTreelet(a_nodes, roots.front(), a_leavesNum, a_aligned,
                        root_idx, newNodes, newOffsetByOld, trChild, trQueue, trVector);
      for (size_t j = 0; j < trChild.size();j++)
        roots.push(trChild[j]);
      }
    else
      {
      if (super_root_idx.empty() || int(newNodes.size()) > super_root_idx.back())
        super_root_idx.emplace_back(newNodes.size()); // the first free element => the root of SUPERTREELET goes there
 
      ProcessOneSuperTreelet(a_nodes, roots.front(), a_leavesNum, a_aligned,
                             root_idx, newNodes, newOffsetByOld, trChild, trChild2, trQueue, trVector);

      if (level == TreeletizLevel::SUPERSUPERTREELETS)
        {
        // "three level treeletization" i.e. by super-supertreelets (=the child super-treelets 
        // of the current one are placed next to the root supertreelet in memory; inside
        // each super-treelet first goes the nodes of its root treelet, then its child treelets))
        trGrandChild.clear();
        for (size_t j = 0; j < trChild.size(); j++)
          ProcessOneSuperTreelet(a_nodes, trChild[j], a_leavesNum, a_aligned,
                                 root_idx, newNodes, newOffsetByOld, trGrandChild, trChild2, trQueue, trVector);

        for (size_t k = 0; k < trGrandChild.size(); k++)
          roots.push(trGrandChild[k]);
        }
      else if (level == TreeletizLevel::SUPERTREELETS)
        {
        // "two-level treeletization" (by supertreelets)
        for (size_t j = 0; j < trChild.size();j++)
          roots.push(trChild[j]);
        }
      }
  roots.pop();
  }

  for(size_t i=0;i<newNodes.size();i++)
  {
    auto& currNode = newNodes[i];
    if(!(currNode.offs_left & LEAF_BIT))
      currNode.offs_left = newOffsetByOld[currNode.offs_left];
    if(!(currNode.offs_right & LEAF_BIT))
      currNode.offs_right = newOffsetByOld[currNode.offs_right];
  }

  a_nodes = newNodes;
}


/// Comparator for std::sort() (in descending order)
inline bool greater4uint2(const LiteMath::uint2& a, const LiteMath::uint2& b)
  {
  return a.x > b.x;
  }

#if 0
/// The first node is its root
/// TODO: Replace 'bbox' with ready quantization step?
double ExcessByRounding(const std::vector<BVHNodeFat>& a_nodes, 
                        const std::vector<uint32_t>& set,
                        const AABB &bbox,
                        int max_level)
  {
  int i;
  float s1, s2, ds1, ds2, excess = 0.0;
  // quantization step
  const LiteMath::float3 delta((bbox.mx - bbox.mn) / float(max_level)); /// TODO: replace "/" with multiplication by the inverse

  /// TODO: Likely the greatest excess will be for the "deepest" nodes, so maybe just use them?
  for (i = 0; i < set.size(); i++)
    {
    const BVHNodeFat& node = a_nodes[set[i]];

    /// TODO: better provide and use GetSizeOfChildBoxLeft() etc.!
    const LiteMath::AABB cb1 = GetChildBoxLeft(node);
    const LiteMath::AABB cb2 = GetChildBoxRight(node);
    const LiteMath::float3 size1(cb1.mx - cb1.mn);
    const LiteMath::float3 size2(cb2.mx - cb2.mn);

    // Area of boxes
    s1 = size1.x * size1.y + size1.x * size1.z + size1.y * size1.z;
    s2 = size2.x * size2.y + size2.x * size2.z + size2.y * size2.z;

    // Increase of box area due to encoding this B.B. by integers in [0,max_level]
    // of values such that 0 gives bbox.mn and max_level bbox.mx
    // This is an estimate (assuming that as a result the decoded boundaries move by 
    // half of quantization step = delta/2)
    ds1 = 0.5 * (delta.x * (size1.y + size1.z))
        + 0.5 * (delta.y * (size1.x + size1.z))
        + 0.5 * (delta.z * (size1.x + size1.y));
    ds2 = 0.5 * (delta.x * (size2.y + size2.z))
        + 0.5 * (delta.y * (size2.x + size2.z))
        + 0.5 * (delta.z * (size2.x + size2.y));

    excess = std::max(excess, std::max(ds1 / s1, ds2 / s2));
    }
  return excess;
  }


/// The first node in each group is its root
double PenaltyOfMergeing(const std::vector<BVHNodeFat>& a_nodes, 
                         const std::vector<uint32_t>& set1,
                         const std::vector<uint32_t>& set2)
  {
  const BVHNodeFat& root1 = a_nodes[set1[0]];
  const BVHNodeFat& root2 = a_nodes[set2[0]];
  const AABB cb1_left  = GetChildBoxLeft(root1);
  const AABB cb2_left  = GetChildBoxLeft(root2);
  const AABB cb1_right = GetChildBoxRight(root1);
  const AABB cb2_right = GetChildBoxRight(root2);

  // Common bounding box for both groups
  AABB bbox(cb1_left);
  bbox.Include(cb2_left);
  bbox.Include(cb1_right);
  bbox.Include(cb2_right);

  // max(relative increase of box area) for all nodes due to encoding their B.B. 
  // relative to 'bbox' by integers.
  // This is an estimate (assuming that as a result the decoded boundaries move by 
  // half of quantization step = delta/2)
  const float excess1 = ExcessByRounding(a_nodes, set1, bbox, 63);
  const float excess2 = ExcessByRounding(a_nodes, set2, bbox, 63);
  return std::max(excess1, excess2);
  }
#endif

/// Re-arranges treelets so that they form GROUPS of size maximally close (but <=) 
/// to the target treelet size. Treelets can NOT be split across groups!! I.e. a 
/// group is several adjacent treelets taken as whole
/// 
/// Groups can be later pad with dummy nodes to make group size = target size of treelet.
/// 
/// Re-arranged treelets in trVectorsNew can be grouped starting from 
/// trVectorsNew[0]. The 1st group includes trVectorsNew[0], trVectorsNew[1], ...
/// until the total number of nodes remains <= treeletSize. If adding trVectorsNew[i]
/// makes the total count of nodes to exceed treeletSize, NEW group starts with this
/// trVectorsNew[i] and may include trVectorsNew[i+1], trVectorsNew[i+2], ... again
/// until the total count of nodes is to exceed the treeletSize. Then the 2nd group
/// ends and the next one starts, etc.
/// 
/// The order of nodes in a treelet is not changed (only treelets are permuted).
/// 
/// E.g. if the number of nodes in treelets is {7, 7, 5, 2, 5, 1, 3, 2, 2} then the
/// groups are 7, 7, [5+2], [5+1], [3+2+2]
/// 
/// @param[in]     a_leavesNum - the number of the "bottom" (last along the tree)
///                              treelet nodes. The number of all nodes is then 2*a_leavesNum-1
///                              The actual number of nodes can be LESS! (if some inside
///                              treelet are terminal nodes).
/// @param[out]    trVectors   - source treelets.
///                              trVectors[i] conveys the indices (related to 'a_nodes', 
///                              no "offset" is used!!) of the i-th treelet ([0] is the
///                              root treelet, then come the child treelets). Each
///                              treelet begins with it root, then the child nodes in BFS.
/// @param[out]    trVectorsNew - re-arranged treelets (the root one remaining the first!!)
static void Merge(size_t a_leavesNum,
                  const std::vector< std::vector<uint32_t> > &trVectors,
                  std::vector< std::vector<uint32_t> >       &trVectorsNew)
  {
  const bool root_treelet_first = true; // if the root treelet must remain the first in memory
  const size_t treeletSize = 2 * a_leavesNum - 1; // number of all nodes in treelet (both inner nodes and leaves)

  bool inComplete = false;
  for (size_t j = 0; j < trVectors.size(); j++)
    {
    if (trVectors[j].size() > 0 && trVectors[j].size() < treeletSize)
      {
      inComplete = true;
      break;
      }
    }
  if (!inComplete)
    {
    trVectorsNew = trVectors;
    return;
    }

  std::vector<LiteMath::uint2> child_tr_size2;
  child_tr_size2.reserve(2 * a_leavesNum + 1);
  child_tr_size2.clear();
  for (size_t j = 0; j < trVectors.size(); j++)
    {
    if (trVectors[j].size() > 0)
      child_tr_size2.emplace_back(LiteMath::uint2(trVectors[j].size(), j));
    }

  if (root_treelet_first) // root treelet is present in trVectors
    child_tr_size2[0].x = 1000 + treeletSize; // to "fool" sort()

  // Sort in descending order
  std::sort(child_tr_size2.begin(), child_tr_size2.end(), greater4uint2);

  if (root_treelet_first) // root treelet is present in trVectors
    child_tr_size2[0].x = trVectors[0].size(); // return actual value

  std::vector< std::vector<uint32_t> > group;
  group.resize(child_tr_size2.size());

  for (size_t j = 0; j < child_tr_size2.size(); j++)
    group[j].emplace_back(child_tr_size2[j].y);

  for (size_t j = 0; j < child_tr_size2.size(); j++)
    {
    if (child_tr_size2[j].x >= treeletSize)
      continue;
    for (size_t k = j + 1; k < child_tr_size2.size(); k++)
      {
      if (child_tr_size2[k].x != uint32_t(-1) && child_tr_size2[j].x + child_tr_size2[k].x <= treeletSize)
        {
        const uint32_t from = child_tr_size2[k].y; //to = child_tr_size2[j].y;
        // Check whether the expence of merging is not too high
///        if (PenaltyOfMerging(a_nodes, trVectors[to], trVectors[from]) < 1.25)
        {
        child_tr_size2[j].x += child_tr_size2[k].x;
        child_tr_size2[k].x = -1;

        group[j].emplace_back(from);
        assert(group[k].size() == 1);
        group[k].clear(); // because this group yet had only ONE element
        }
        }
      }
    }

  // Regroup node indices according to the calculated data

  // First count non-empty groups
  int num_groups = 0;
  for (size_t j = 0; j < child_tr_size2.size(); j++)
    num_groups += (child_tr_size2[j].x != uint32_t(-1));

  trVectorsNew.reserve(num_groups);
  trVectorsNew.clear();

  for (size_t j = 0; j < child_tr_size2.size(); j++)
    {
    if (child_tr_size2[j].x == uint32_t(-1))
      continue; // all nodes moved out from it

#ifdef OSE_TESTS
    int nn = 0;
    for (size_t k = 0; k < group[j].size(); k++)
      nn += trVectors[group[j][k]].size();
    assert(nn == child_tr_size2[j].x);
#endif

    for (size_t k = 0; k < group[j].size(); k++)
      {
      const int idx = group[j][k];
      if (idx != -1)
        {
        trVectorsNew.emplace_back(trVectors[idx]);
        assert(trVectors[idx].size() <= treeletSize);
        }
      }
    }
  }


/// Pad array of nodes with 'n' dummy nodes
inline void PadWithDummy(const BVHNodeFat &dummyNode, int n, std::vector<BVHNodeFat>& a_nodes)
  {
  for (int jj = 0; jj < n; jj++)
    a_nodes.emplace_back(dummyNode);
  }

/// Writes a set of treelets in memory. First comes trVectors[0], then trVectors[1], 
/// trVectors[2] etc. 
/// 
/// Determines the boundaries of GROUPS (sets of treelets that have size maximally 
/// close (but <= treeletSize) of treelets, and puts to grStart.
/// 
/// The 1st group includes trVectors[0], trVectors[1], ...
/// until the total number of nodes remains <= treeletSize. If adding trVectorsNew[i]
/// makes the total count of nodes to exceed treeletSize, NEW group starts with this
/// trVectors[i] and may include trVectors[i+1], trVectors[i+2], ... again
/// until the total count of nodes is to exceed the treeletSize. Then the 2nd group
/// ends and the next one starts, etc.
/// 
/// E.g. if the number of nodes in treelets is {7, 7, 5, 2, 5, 1, 3, 2, 2} then the
/// groups are 7, 7, [5+2], [5+1], [3+2+2]
/// 
/// @param[in]     a_nodes - array of nodes of the BVH tree
/// @param[out]    trVectors   - trVectors[i] = the indices (related to 'a_nodes', 
///                              no "offset" is used!!) of the nodes in the i-th treelet.
/// @param[in]     a_leavesNum - the number of the "bottom" (last along the tree)
///                              treelet nodes. The number of all nodes is then 2*a_leavesNum-1
///                              The actual number of nodes can be LESS! (if some inside
///                              treelet are terminal nodes).
/// @param[in]     a_aligned         - whether or not to pad the set of treelet nodes in memory
///                                    by "dummy" elements to keep the size of the group
///                                    exactly the same (as for the ideal tree)
/// @param[in,out] newNodes - accumulated array of tree nodes in the NEW ORDER.
///                           The nodes of the treelet created are APPENDED to it,
///                           first comes the treelet's root node.
/// @param[in,out] newOffsetByOld - lookup from the index of a node in old order and
///                                 the same node in the new order, i.e.
///                                 newNodes[newOffsetByOld[i]] = a_nodes[i].
///                                 
///                                 Before the very first call to ProcessOneTreelet() 
///                                 this array must be allocated (by a_nodes.size()) by the caller!!! 
///                                 And filled with 0xFFFFFFFF.
/// 
///                                 The data for nodes of the current set is written
///                                 into the related elements, the rest remain unchanged.
/// @param[in,out] root_idx - accumulated array of the (indices of) treelet roots.
///                           The indices relate to the NEW ORDER of nodes, i.e.
///                           the root os the i-th treelet is newNodes[root_idx[i]].
///                           The root of the created treelet is APPENDED to this array.
/// @param[out]    grStart  - APPENDED to it are the indices of the starting nodes of GROUPS, related to newNodes(new positions)
///                           (similar to root_idx)
static void AppendNodes(const std::vector<BVHNodeFat>& a_nodes,
                        const std::vector< std::vector<uint32_t> > &trVectors,
                        int treeletSize,
                        bool a_aligned,
                        std::vector<BVHNodeFat> &newNodes,
                        std::vector<uint32_t> &newOffsetByOld,
                        std::vector<int>& root_idx,
                        std::vector<int>& grStart)
  {
  static const BVHNodeFat dummyNode = FatBVH::BVHNodeFatDummy();

  int nn = 0; // accumulated number of nodes in current group
  for (size_t i = 0; i < trVectors.size(); i++)
    {
    // After trVectorsNew[i] the size of the group would become:
    const int groupSizeAfter = nn + trVectors[i].size();
    const int accumulatedGroupSize = newNodes.size() - (grStart.empty() ? 0 : grStart.back());
    assert(accumulatedGroupSize <= treeletSize);
    assert(nn == 0 || accumulatedGroupSize == nn);

    if (nn == 0 || groupSizeAfter > treeletSize)
      {
      // Previous group ends or is absent.
      // The trVectorsNew[i] goes/[will go] into new (yet empty) room
            
      // If nn!=0, trVectorsNew[i] can not fit in the current group, so we shall
      // start a new group for it. This new group naturally begins after 
      // the last current element of newNodes.
      // If nn=0 this is a new group unconditionally.

      if (a_aligned)
        {
        if (!grStart.empty() && accumulatedGroupSize < treeletSize)
          {
          // current group of length accumulatedGroupSize ends. 
          // Pad it with dummy nodes to keep group size = treeletSize
          PadWithDummy(dummyNode, treeletSize - accumulatedGroupSize, newNodes);
          }
        }
      grStart.emplace_back(newNodes.size()); // next group starts at this position
      }

    if (groupSizeAfter < treeletSize)
      nn += trVectors[i].size(); // trVectorsNew[i] goes into CURRENT group

    root_idx.emplace_back(newNodes.size()); // next treelet starts at this position
    AppendNodes(a_nodes, trVectors[i], newNodes, newOffsetByOld);
    assert(newNodes.size() <= grStart.back() + treeletSize);

    if (groupSizeAfter >= treeletSize)
      nn = 0; // start next group
    }
  // If needed, pad the last group with dummy nodes to keep its size = treeletSize
  const int last_gr_size = newNodes.size() - grStart.back();
  if (last_gr_size < treeletSize)
    {
    // Pad with dummy nodes to keep group size = treeletSize
    PadWithDummy(dummyNode, treeletSize - last_gr_size, newNodes);
    }
  }

/// Reorders the nodes of the BVH tree in memory so that they are grouped.
/// 
/// Works in for SUPERTREELETS mode!!
/// 
/// Allows to GROUP treelets (can NOT be split across groups!! I.e. a group may include 
/// some treelets as whole!) in groups of size = target size of treelet. Groups
/// can be pad with dummy nodes to make group size = target size of treelet.
/// 
/// Without the "merge" mode the groups contain one and only one treelet!! I.e.
/// group = treelet then!
/// 
/// Supertreelet = treelet and its child treelets (those whose roots are outside the
/// 1st treelet but their parents are INSIDE it). First come the elements of the
/// "root treelet", then its child treelets (along the horizontal tree slices)
/// These child treelets can be reordered in memory (against their initial order in the tree)
/// so that first come those with > target size of treelet, then smaller ones, arranged
/// so that successive treelets for GROUP of size maximally close (but smaller)
/// to the target size of treelet. Then they can be padded with dummy nodes to make
/// group size = target size of treelet.
/// 
/// Child SUPERtreelets (supertreelets outside this supertreelet but such that 
/// PARENTS of their roots belong to this supertreelet) also can be grouped.
/// Again in memory first come thos child supertreelets whose size >= target size of treelet.
/// Smaller ones are re-arranged so that they can be GROUPED with those next to them
/// in memory. Again each incomplete group can be pad with dummy nodes to keep its size.
/// 
/// NOTICE The root treelet always comes the first!!
/// 
/// @param[in,out]  a_nodes - array of nodes of the BVH tree. On input, in the original
///                           order. On EXIT, grouped by [[super]super]treelets.
/// @param[out]    root_idx - indices of treelet roots, related to the re-arranged array (new positions)
/// @param[out]    super_root_idx - indices of supertreelet roots, related to the re-arranged array (new positions)
/// @param[out]    grStart     - indices of the starting nodes of GROUPS, related to the re-arranged array (new positions)
/// @param[in]     a_leavesNum - the number of the "bottom" (last along the tree)
///                              treelet nodes. The number of all nodes is then 2*a_leavesNum-1
///                              The actual number of nodes can be LESS! (if some inside
///                              treelet are terminal nodes).
/// @param[in]     a_aligned   - whether or not to pad groups by "dummy" elements to keep 
///                              the size of the group = the target (ideal) treelet size
/// @param[in]     a_merge     - if 'false', group = treelet. 
///                              if 'true' several small treelets can be "merged"
///                              into a GROUP. 
void FatBVH::ReorderTRBNew(std::vector<BVHNodeFat>& a_nodes, std::vector<int>& root_idx, std::vector<int>& super_root_idx, std::vector<int>& grStart, size_t a_leavesNum, bool a_aligned, bool a_merge)
{
  const size_t treeletSize        = 2*a_leavesNum-1; // number of all nodes in treelet (both inner nodes and leaves)

  //const BVHNodeFat dummyNode = BVHNodeFatDummy();

  grStart.reserve(2 * a_nodes.size() / treeletSize);
  grStart.clear();

  std::vector<BVHNodeFat> newNodes;
  newNodes.reserve(a_nodes.size());
  newNodes.resize(0);

  std::vector<uint32_t> newOffsetByOld(a_nodes.size());
  for(auto& index : newOffsetByOld) // mark nodes are not processed yet
    index = 0xFFFFFFFF;

  std::vector<bool> wkIdx;
  wkIdx.resize(a_nodes.size());
  for (size_t i = 0; i < wkIdx.size(); i++) // mark nodes are not processed yet
    wkIdx[i] = false;

  std::queue<uint32_t>  trQueue;
  std::vector<uint32_t> trVector;
  trVector.reserve(treeletSize); // +1 is for root node inside treelet

  uint32_t curr_parent = 0xFFFFFFFF;

  // Used for supertreelets, keeps the child root node as .x and 
  // the root node of the parent supertreelet as .y
  std::queue<LiteMath::uint2> roots;
  roots.push(LiteMath::uint2(LiteMath::uint(0), LiteMath::uint(0)));

  std::vector<uint32_t> trChild, trChild2;

  // nodes of supertreelet, first the root treelet then (maybe regrouped) child treelets
  std::vector< std::vector<uint32_t> > strVectors, strVectorsNew;
  std::vector< std::vector<uint32_t> > nodes_of_child_treelets, nodes_of_child_treelets_new;

  while (roots.size() != 0)
    {
    trChild.clear();

    if (super_root_idx.empty() || int(newNodes.size()) > super_root_idx.back())
      super_root_idx.emplace_back(newNodes.size()); // the first free element => the root of SUPERTREELET goes there

    curr_parent = roots.front()[1]; // index of root node of the PARENT supertreelet!
    CollectOneSuperTreeletNodes(a_nodes, roots.front()[0], a_leavesNum,
                                strVectors, trChild, trChild2, trQueue, trVector, wkIdx);

    // Re-arrange treelets in strVectors so that they form GROUPS maximally close
    // is size to treeletSize. The root treelet remains the first!
    if (a_merge)
      Merge(a_leavesNum, strVectors, strVectorsNew);
    else
      strVectorsNew = strVectors;

    // The total number of nodes in the super-treelet
    int num_nodes = 0;
    for (size_t i = 0; i < strVectorsNew.size(); i++)
      num_nodes += strVectorsNew[i].size();

    if (num_nodes >= int(treeletSize))
    {
      // Put new nodes to memory
      AppendNodes(a_nodes, strVectorsNew, treeletSize, a_aligned, newNodes, newOffsetByOld, root_idx, grStart);
    }
    else
    {
      // Too few nodes; remember it for future merging
      for (size_t i = 0; i < strVectorsNew.size(); i++)
        nodes_of_child_treelets.emplace_back(strVectorsNew[i]);
    }

    // Use the index of the root node of THIS supertreelet to identify te PARENT 
    // for its child supertreelets
    for (size_t j = 0; j < trChild.size();j++)
      roots.push(LiteMath::uint2(LiteMath::uint(trChild[j]), LiteMath::uint(roots.front()[0]))); /// remember ID of the parent supertreelet

    roots.pop();

    // Complete the remained SMALL child supertreelets, if some had not yet processed
    if (roots.size() == 0 || (roots.size() > 0 && roots.front()[1] != curr_parent))
      {
      // The next supertreelet will be a child of ANOTHER parent.
      // Finish with the pending SMALL child supertreelets of the same (previous) parent
      if (nodes_of_child_treelets.size() > 0)
        {
        if (a_merge)
          Merge(a_leavesNum, nodes_of_child_treelets, nodes_of_child_treelets_new);
        else
          nodes_of_child_treelets_new = nodes_of_child_treelets;

#ifdef OSE_TESTS
        int num_nodes = 0;
        for (int i = 0; i < nodes_of_child_treelets_new.size(); i++)
          num_nodes += nodes_of_child_treelets_new[i].size();
#endif

        // save them to memory
        if (super_root_idx.empty() || int(newNodes.size()) > super_root_idx.back())
          super_root_idx.emplace_back(int(newNodes.size()));

        AppendNodes(a_nodes, nodes_of_child_treelets_new, treeletSize, a_aligned,
          newNodes, newOffsetByOld, root_idx, grStart);
        }

      nodes_of_child_treelets.clear();

      if (roots.size() > 0)
        curr_parent = roots.front()[1]; // Use index of the root node as ID of the supertreelet
      }
    }

  for (size_t i = 0; i < newNodes.size(); i++)
  {
    auto& currNode = newNodes[i];
    if(!(currNode.offs_left & LEAF_BIT))
      currNode.offs_left = newOffsetByOld[currNode.offs_left];
    if(!(currNode.offs_right & LEAF_BIT))
      currNode.offs_right = newOffsetByOld[currNode.offs_right];
  }

#ifdef OSE_TESTS
  printf("%i groups of size <=1\n", num0);
  for (int i = 0; i < newOffsetByOld.size(); i++)
    {
    if (newOffsetByOld[i] < 0 || newOffsetByOld[i] > newNodes.size() - 1)
      assert(false);
    }


  std::vector<uint32_t> hist;
  hist.resize(treeletSize + 1);
  for (int i = 0; i < grStart.size(); i++)
    {
    const int groupSize = (i + 1 < grStart.size() ? grStart[i + 1] : newNodes.size())
                        - grStart[i];
    assert(groupSize > 0 && groupSize <= treeletSize);
    if (a_aligned && groupSize < treeletSize)
      i = i;
    hist[groupSize]++;
    }
#endif

  if (a_aligned && (newNodes.size() % treeletSize) != 0)
  {
    printf("ERROR %i %i\n", int(newNodes.size()), int(treeletSize));
    assert(false);
  }

  a_nodes = newNodes;
  return;
}


