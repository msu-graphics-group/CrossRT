#include <fstream>
#include <queue>

#include "cbvh.h"
#include "cbvh_core.h"
#include "raytrace_common.h"

using cbvh2::BVHNode;
using LiteMath::float4;
extern double   g_buildTime; // used to get build time from deep underground of code
extern uint64_t g_buildTris; // used to get total tris processed by builder

cbvh2::BuilderPresets cbvh2::BuilderPresetsFromString(const char* a_str)
{
  cbvh2::BuilderPresets presets = cbvh2::BuilderPresets();

  const std::string a_buildName(a_str);
  if(a_buildName.find("cbvh_hq") != std::string::npos)
    presets.quality = cbvh2::BVH_CONSTRUCT_QUALITY;
  else if(a_buildName.find("cbvh_med") != std::string::npos)
    presets.quality = cbvh2::BVH_CONSTRUCT_MEDIUM;
  else if(a_buildName.find("cbvh_embree_fast") != std::string::npos)
    presets.quality = cbvh2::BVH_CONSTRUCT_EMBREE_FAST;
  else if(a_buildName.find("cbvh_embree") != std::string::npos)
    presets.quality = cbvh2::BVH_CONSTRUCT_EMBREE;
  else if(a_buildName.find("nanort") != std::string::npos || a_buildName.find("nano_rt") != std::string::npos || a_buildName.find("NanoRT") != std::string::npos)
    presets.quality = cbvh2::BVH_CONSTRUCT_NANORT;
  
  char symb = a_buildName[a_buildName.size()-1];
  if(std::isdigit(symb))
    presets.primsInLeaf = int(symb) - int('0');
  
  if(a_buildName.find("lbvh_gpu") != std::string::npos || a_buildName.find("cbvh_lbvh_gpu") != std::string::npos)
  {
    presets.quality     = cbvh2::BVH_CONSTRUCT_FAST_GPU;
    presets.primsInLeaf = 1;
  }
  else if(a_buildName.find("lbvh") != std::string::npos || a_buildName.find("cbvh_lbvh") != std::string::npos)
  {
    presets.quality     = cbvh2::BVH_CONSTRUCT_FAST;
    presets.primsInLeaf = 1;
  }

  return presets;
}

static inline void CheckNodesAlign(const std::vector<cbvh::BVHNode>& a_nodes, int a_alignCoeff)
{
  bool isOk = true;
  for(size_t i=0;i<a_nodes.size();i++) 
  {
    const auto& node = a_nodes[i];
    if((node.leftOffset & LEAF_BIT) == 0) 
    {
      if(node.leftOffset % a_alignCoeff != 0)
        isOk = false;
    }
  }

  const char* timeFileName = "z_align_check.csv";
  std::ofstream fout;
  std::ifstream fin(timeFileName);
  if (!fin.is_open())
  {
    fout.open(timeFileName);
    fout << "Align;" << std::endl;
  }
  else
  {
    fin.close();
    fout.open(timeFileName, std::ios::app);
  }

  if(isOk)
    fout << "OK;" << std::endl;
  else
    fout << "FAILED!;" << std::endl;
}

struct TraverseDFS
{
  TraverseDFS(const std::vector<cbvh::BVHNode>& a_nodesOld, std::vector<cbvh::BVHNode>& a_nodesNew) : oldNodes(a_nodesOld), newNodes(a_nodesNew) {}
  
  uint32_t VisitNode(uint32_t a_nodeOffsetLeft, uint32_t a_parentEscapeIndex)
  {
    if((a_nodeOffsetLeft & LEAF_BIT) != 0)
      return a_nodeOffsetLeft;

    const auto& nodeLeft  = oldNodes[a_nodeOffsetLeft+0];
    const auto& nodeRight = oldNodes[a_nodeOffsetLeft+1];

    if(newNodes.size()%2 != 0) {
      newNodes.push_back(cbvh2::DummyNode()); // never happends?
      alignCounter++;
    }
    
    newNodes.push_back(nodeLeft);
    newNodes.push_back(nodeRight);
    
    const uint32_t leftOffsetNew = uint32_t(newNodes.size()-2);

    newNodes[leftOffsetNew+0].escapeIndex = leftOffsetNew+1;
    newNodes[leftOffsetNew+1].escapeIndex = a_parentEscapeIndex;

    newNodes[leftOffsetNew+0].leftOffset  = VisitNode(nodeLeft.leftOffset, leftOffsetNew+1);
    newNodes[leftOffsetNew+1].leftOffset  = VisitNode(nodeRight.leftOffset, a_parentEscapeIndex);
  
    return leftOffsetNew;
  }

  uint32_t alignCounter = 0;

protected:
  const std::vector<cbvh::BVHNode>& oldNodes;
  std::vector<cbvh::BVHNode>&       newNodes;
};

static inline std::vector<cbvh::BVHNode> AlignNodes2(const std::vector<cbvh::BVHNode>& a_nodes)
{
  if(a_nodes.size() == 0)
    return a_nodes;

  std::vector<cbvh::BVHNode> newNodes;
  newNodes.reserve(a_nodes.size()*2);
 
  if(a_nodes.size() == 1)
  {
    newNodes.resize(2);
    newNodes[0] = a_nodes[0];
    newNodes[1] = cbvh2::DummyNode();
    return newNodes;
  }

  TraverseDFS trav(a_nodes, newNodes);
  trav.VisitNode(a_nodes[0].leftOffset, cbvh::ESCAPE_ROOT);
  
  if(trav.alignCounter != 0)
    std::cout << "[AlignNodes2]: " << trav.alignCounter << "nodes was aligned!" << std::endl;
  
  if(newNodes.size()%2 != 0)
    newNodes.push_back(cbvh2::DummyNode());

  return newNodes;
}


struct TraverseDFS4
{
  TraverseDFS4(const std::vector<cbvh::BVHNode>& a_nodesOld, std::vector<cbvh::BVHNode>& a_nodesNew) : oldNodes(a_nodesOld), newNodes(a_nodesNew) {}
  
  uint32_t VisitNode(uint32_t a_nodeOffsetLeft, uint32_t a_parentEscapeIndex)
  {
    if((a_nodeOffsetLeft & LEAF_BIT) != 0)
      return a_nodeOffsetLeft;

    const auto& node0 = oldNodes[a_nodeOffsetLeft+0];
    const auto& node1 = oldNodes[a_nodeOffsetLeft+1];
    const auto& node2 = oldNodes[a_nodeOffsetLeft+2];
    const auto& node3 = oldNodes[a_nodeOffsetLeft+3];

    while(newNodes.size()%4!=0) {
      newNodes.push_back(cbvh2::DummyNode()); // never happends?
      alignCounter++;
    }
    
    newNodes.push_back(node0);
    newNodes.push_back(node1);
    newNodes.push_back(node2);
    newNodes.push_back(node3);
    
    const uint32_t leftOffsetNew = uint32_t(newNodes.size()-4);

    newNodes[leftOffsetNew+0].escapeIndex = leftOffsetNew+1;
    newNodes[leftOffsetNew+1].escapeIndex = leftOffsetNew+2;
    newNodes[leftOffsetNew+2].escapeIndex = leftOffsetNew+3;
    newNodes[leftOffsetNew+3].escapeIndex = a_parentEscapeIndex;

    newNodes[leftOffsetNew+0].leftOffset  = VisitNode(node0.leftOffset, leftOffsetNew+1);
    newNodes[leftOffsetNew+1].leftOffset  = VisitNode(node1.leftOffset, leftOffsetNew+2);
    newNodes[leftOffsetNew+2].leftOffset  = VisitNode(node2.leftOffset, leftOffsetNew+3);
    newNodes[leftOffsetNew+3].leftOffset  = VisitNode(node3.leftOffset, a_parentEscapeIndex);
  
    return leftOffsetNew;
  }

  uint32_t alignCounter = 0;

protected:
  const std::vector<cbvh::BVHNode>& oldNodes;
  std::vector<cbvh::BVHNode>&       newNodes;
};

static inline std::vector<cbvh::BVHNode> AlignNodes4(const std::vector<cbvh::BVHNode>& a_nodes)
{
  if(a_nodes.size() == 0)
    return a_nodes;

  std::vector<cbvh::BVHNode> newNodes;
  newNodes.reserve(a_nodes.size()*2);
 
  TraverseDFS4 trav(a_nodes, newNodes);
  trav.VisitNode(a_nodes[0].leftOffset, cbvh::ESCAPE_ROOT);
  
  if(trav.alignCounter != 0)
    std::cout << "[AlignNodes4]: " << trav.alignCounter << "nodes was aligned!" << std::endl;
  
  while(newNodes.size()%4 != 0)
    newNodes.push_back(cbvh2::DummyNode());

  return newNodes;
}

struct Converter2to4
{
  Converter2to4(const std::vector<cbvh2::BVHNode>& in_nodes, std::vector<cbvh2::BVHNode>& out_nodes) : m_inNodes(in_nodes), m_outNodes(out_nodes) 
  {
    m_outNodes.reserve(m_inNodes.size());
    m_outNodes.resize(0);
  }

  static inline uint32_t EXTRACT_START(uint32_t a_leftOffset)  { return  a_leftOffset & START_MASK; }
  static inline uint32_t EXTRACT_COUNT(uint32_t a_leftOffset)  { return (a_leftOffset & SIZE_MASK) >> 24; }

  uint32_t Visit(uint32_t a_offset) //
  {
    if((a_offset & LEAF_BIT) != 0)
      return a_offset;

    const cbvh2::BVHNode leftNode  = m_inNodes[a_offset + 0];
    const cbvh2::BVHNode rightNode = m_inNodes[a_offset + 1];

    if((leftNode.leftOffset & LEAF_BIT) != 0 && (rightNode.leftOffset & LEAF_BIT) != 0)
    {
      const uint32_t start1 = EXTRACT_START(leftNode.leftOffset);
      const uint32_t count1 = EXTRACT_COUNT(leftNode.leftOffset);

      const uint32_t start2 = EXTRACT_START(rightNode.leftOffset);
      const uint32_t count2 = EXTRACT_COUNT(rightNode.leftOffset);
      
      const uint32_t start = std::min(start1, start2);
      const uint32_t count = count1 + count2;
      
      return PackOffsetAndSize(start, count);
    }

    const uint32_t currSize = m_outNodes.size();
    m_outNodes.resize(currSize + 4);

    if((leftNode.leftOffset & LEAF_BIT) == 0 && (rightNode.leftOffset & LEAF_BIT) == 0) // directly get 4 nodes
    {
      m_outNodes[currSize + 0] = m_inNodes[leftNode.leftOffset + 0];
      m_outNodes[currSize + 1] = m_inNodes[leftNode.leftOffset + 1];
      m_outNodes[currSize + 2] = m_inNodes[rightNode.leftOffset + 0];
      m_outNodes[currSize + 3] = m_inNodes[rightNode.leftOffset + 1];
    }
    else if((leftNode.leftOffset & LEAF_BIT) == 0)
    {
      m_outNodes[currSize + 0] = m_inNodes[leftNode.leftOffset + 0];
      m_outNodes[currSize + 1] = m_inNodes[leftNode.leftOffset + 1];
      m_outNodes[currSize + 2] = rightNode;
      m_outNodes[currSize + 3] = cbvh2::DummyNode();
    }
    else if((rightNode.leftOffset & LEAF_BIT) == 0)
    {
      m_outNodes[currSize + 0] = leftNode;
      m_outNodes[currSize + 1] = m_inNodes[rightNode.leftOffset + 0];
      m_outNodes[currSize + 2] = m_inNodes[rightNode.leftOffset + 1];
      m_outNodes[currSize + 3] = cbvh2::DummyNode();
    }
    
    m_outNodes[currSize + 0].leftOffset = Visit(m_outNodes[currSize + 0].leftOffset);
    m_outNodes[currSize + 1].leftOffset = Visit(m_outNodes[currSize + 1].leftOffset);
    m_outNodes[currSize + 2].leftOffset = Visit(m_outNodes[currSize + 2].leftOffset);
    m_outNodes[currSize + 3].leftOffset = Visit(m_outNodes[currSize + 3].leftOffset);

    return currSize;
  }

  const std::vector<cbvh2::BVHNode>& m_inNodes;
  std::vector<cbvh2::BVHNode>& m_outNodes;
};

static inline std::vector<cbvh2::BVHNode> BVH2ToBVH4(const std::vector<cbvh2::BVHNode>& a_nodes)
{
  std::vector<cbvh2::BVHNode> result;
  if(a_nodes.size() == 2)
  {
    result = a_nodes;
    result.push_back(cbvh2::DummyNode());
    result.push_back(cbvh2::DummyNode());
  }
  else
  {
    Converter2to4 converter(a_nodes, result);
    converter.Visit(0);
  }
  return result;
}


cbvh2::BVHTreeCommon cbvh2::BuildBVH(const float* a_vpos3f,     size_t a_vertNum, size_t a_vByteStride, 
                                     const uint32_t* a_indices, size_t a_indexNum, BuilderPresets a_presets)
{
  const float4* inVertices = (const float4*)a_vpos3f;
  std::vector<float4> vertDataTemp;

  g_buildTris += a_indexNum/3;

  if(a_vByteStride != 16)
  { 
    vertDataTemp.resize(a_vertNum);
    const size_t stride = a_vByteStride/4;
    for(size_t i=0;i<a_vertNum;i++) {
      vertDataTemp[i].x = a_vpos3f[i*stride + 0];
      vertDataTemp[i].y = a_vpos3f[i*stride + 1];
      vertDataTemp[i].z = a_vpos3f[i*stride + 2];
      vertDataTemp[i].w = 1.0f;
    }
    inVertices = vertDataTemp.data();
  }
  
  cbvh::BVHPresets presets1;
  presets1.primsInLeaf = a_presets.primsInLeaf;
  
  bool convert2to4 = false;
  if(a_presets.quality == BVH_CONSTRUCT_EMBREE && a_presets.format == cbvh2::BVH4_LEFT_OFFSET)
  {
    a_presets.format = cbvh2::BVH2_LEFT_OFFSET;
    convert2to4 = true;
  }

  switch (a_presets.format)
  {
  case cbvh2::BVH2_LEFT_OFFSET:
    presets1.childrenNum   = 2;
    presets1.desiredFormat = cbvh::FMT_BVH2Node32_Interval32_Static;
    break;

  case cbvh2::BVH2_LEFT_RIGHT:
    presets1.childrenNum   = 2;
    presets1.desiredFormat = cbvh::FMT_BVH2Node32_Interval32_Dynamic;
    break;

  case cbvh2::BVH4_LEFT_OFFSET:
    presets1.childrenNum   = 4;
    presets1.desiredFormat = cbvh::FMT_BVH4Node32_Interval32_Static;
    break;
  
  default:
    {
      std::cout << "[cbvh2::BuildBVH] ERROR: wrong input format = " << int(a_presets.format) << std::endl;
      return cbvh2::BVHTreeCommon();
    }
    break;
  }

  presets1.btype = cbvh::BVH_BUILDER_TYPE(int(a_presets.quality));  
  auto bvhData   = BuildBVH(inVertices, a_vertNum, a_indices, a_indexNum, presets1);
  for(size_t i=0;i<bvhData.nodes.size();i++) {
    if(int(bvhData.nodes[i].leftOffset) < 0) 
      bvhData.nodes[i].leftOffset = PackOffsetAndSize(bvhData.intervals[i].start, bvhData.intervals[i].count);
  }
  
  if(a_presets.format == cbvh2::BVH4_LEFT_OFFSET) //
  {
    auto nodesAligned = AlignNodes4(bvhData.nodes);
    //CheckNodesAlign(nodesAligned, 4);
    return BVHTreeCommon(nodesAligned, bvhData.indicesReordered);  
  }
  else if(convert2to4)
  {
    auto nodesAligned = BVH2ToBVH4(AlignNodes2(bvhData.nodes));
    return BVHTreeCommon(nodesAligned, bvhData.indicesReordered);
  }
  else
  {
    if(a_presets.format == cbvh2::BVH2_LEFT_RIGHT) // assume we don't have to aligh "BVH2_LEFT_RIGHT" 
    {
      if(a_presets.quality == cbvh2::BVH_CONSTRUCT_FAST || a_presets.quality == cbvh2::BVH_CONSTRUCT_FAST_GPU)  // LBVH builder already store rightOffset in 'escapeIndex'
        return BVHTreeCommon(bvhData.nodes, bvhData.indicesReordered);
      else
      {
        for(auto& node : bvhData.nodes)            // other builders save true escapeIndex in 'escapeIndex' field
          node.escapeIndex = node.leftOffset + 1;
        return BVHTreeCommon(bvhData.nodes, bvhData.indicesReordered);  
      }
    }
    else if(a_presets.format == cbvh2::BVH2_LEFT_OFFSET && (a_presets.quality == cbvh2::BVH_CONSTRUCT_FAST || a_presets.quality == cbvh2::BVH_CONSTRUCT_FAST_GPU)) // already aligned
      return BVHTreeCommon(bvhData.nodes, bvhData.indicesReordered);
    else
    {
      auto nodesAligned = AlignNodes2(bvhData.nodes);
      //CheckNodesAlign(nodesAligned, 2);
      return BVHTreeCommon(nodesAligned, bvhData.indicesReordered);
    }
  }
}

struct TraverseDFS_ESCI2
{
  TraverseDFS_ESCI2(const std::vector<cbvh::BVHNode>& a_nodesOld, 
                          std::vector<uint32_t>&      a_escape) : oldNodes(a_nodesOld), escape(a_escape) {}
  
  void VisitNode(uint32_t a_nodeOffset, uint32_t a_parentESCI)
  {
    if((a_nodeOffset & LEAF_BIT) != 0)
      return;
    escape[a_nodeOffset] = a_parentESCI;
    const auto& currNode = oldNodes[a_nodeOffset];
    VisitNode(currNode.leftOffset, currNode.escapeIndex);
    VisitNode(currNode.escapeIndex, a_parentESCI);
  }

protected:
  const std::vector<cbvh::BVHNode>& oldNodes;
  std::vector<uint32_t>&            escape;
};

static inline std::vector<cbvh::BVHNode> RightOffsToEscapeIndex(const std::vector<cbvh::BVHNode>& a_tree) // todo: make this algorithm non recursive
{
  std::vector<cbvh::BVHNode> res(a_tree.size());
  std::vector<uint32_t>      escape(a_tree.size());

  //for(uint32_t i=0;i<uint32_t(parent.size());i++) 
  //{
  //  parent[i] = uint32_t(-1);
  //  escape[i] = cbvh::ESCAPE_ROOT;
  //}
  //
  //for(uint32_t i=0;i<uint32_t(a_tree.size());i++)
  //{
  //  const auto& currNode  = a_tree[i];
  //  const bool isLeaf = ((currNode.leftOffset & LEAF_BIT) != 0);
  //  if(!isLeaf)
  //  {
  //    parent[currNode.leftOffset]  = i;
  //    parent[currNode.escapeIndex] = i;
  //    escape[currNode.leftOffset]  = currNode.escapeIndex; // escape(left) = right
  //  }
  //}
  
  TraverseDFS_ESCI2 trav(a_tree, escape);
  trav.VisitNode(0, cbvh::ESCAPE_ROOT); // todo: implement tree top-down and bottom-up traversal in a different way in loops

  for(uint32_t i=0;i<uint32_t(a_tree.size());i++)
  {
    auto currNode        = a_tree[i];
    currNode.escapeIndex = escape[i];
    res[i]               = currNode;
  }

  return res;
}



std::vector<cbvh2::BVHNode> cbvh2::BuildBVH(const cbvh2::BVHNode* a_nodes, size_t a_objNum, cbvh2::BuilderPresets a_presets)
{
  cbvh::BVHPresets presets1;
  presets1.primsInLeaf = 1;
  presets1.btype       = cbvh::BVH_BUILDER_TYPE(int(a_presets.quality));

  switch (a_presets.format)
  {
  case cbvh2::BVH2_LEFT_OFFSET:
    {
      presets1.childrenNum   = 2;
      presets1.desiredFormat = cbvh::FMT_BVH2Node32_Interval32_Static;
    }
    break;

  case cbvh2::BVH2_LEFT_RIGHT:
  case cbvh2::BVH2_LEFT_ROPES:
    presets1.childrenNum   = 2;
    presets1.desiredFormat = cbvh::FMT_BVH2Node32_Interval32_Dynamic;
    break;

  case cbvh2::BVH4_LEFT_OFFSET:
    {
      std::cout << "[cbvh2::BuildBVH] ERROR: 'BVH4_LEFT_OFFSET' format for input bboxes in not implemented! " << std::endl;
      presets1.childrenNum   = 4;
      presets1.desiredFormat = cbvh::FMT_BVH4Node32_Interval32_Static;
    }
    break;
  
  default:
    {
      std::cout << "[cbvh2::BuildBVH] ERROR: wrong input format = " << int(a_presets.format) << std::endl;
      return std::vector<cbvh2::BVHNode>();
    }
    break;
  }
  

  const cbvh2::BVHNode* input = a_nodes;
  std::vector<cbvh2::BVHNode> copy; 
  
  if(a_presets.quality != cbvh2::BVH_CONSTRUCT_FAST)  // set indices inside 'box-nodes' for internal builder
  {
    copy = std::vector<cbvh2::BVHNode>(a_nodes, a_nodes+a_objNum);
    for(size_t i=0;i<copy.size();i++)
    {
      copy[i].leftOffset  = uint32_t(i);
      copy[i].escapeIndex = 1;
    }
    input = copy.data();
  }

  auto lbvh  = cbvh::BuildBVH(input, a_objNum, presets1);
  auto nodes = (a_presets.format == cbvh2::BVH2_LEFT_ROPES) ? RightOffsToEscapeIndex(lbvh.nodes) : lbvh.nodes;
  
  if(!lbvh.invalidIntervals) // restricted way for TLAS, single BLAS per leaf
  {
    bool isNotOne = false;
    for(size_t i=0;i<nodes.size();i++) 
    {
      if(int(nodes[i].leftOffset) < 0 && nodes[i].leftOffset != 0x80000000)
      {
        const uint32_t objectId = lbvh.indicesReordered[lbvh.intervals[i].start];
        nodes[i].leftOffset     = PackOffsetAndSize(objectId, lbvh.intervals[i].count);
        if(lbvh.intervals[i].count > 1)
          isNotOne = true;
      }
    }
  
    if(isNotOne)
      std::cout << "[cbvh2::BuildBVH] warning: leaf contains more than one primitive!" << std::endl;
  }

  if(a_presets.format == cbvh2::BVH2_LEFT_OFFSET && a_presets.quality != cbvh2::BVH_CONSTRUCT_FAST) // lbvh (cbvh2::BVH_CONSTRUCT_FAST) already convert nodes inside 
  {
    auto nodesAligned = AlignNodes2(nodes);
    //CheckNodesAlign(nodesAligned, 2);
    return nodesAligned;
  }
  else
    return nodes;
}