#include <iostream>
#include <queue>
#include <stack>
#include <limits>
#include <cassert>
#include <cmath>
#include <cfloat>

#include <algorithm>
#include <unordered_set>
#include <fstream>

#include "cbvh_core.h"
#include "cbvh_internal.h"

using LiteMath::float2;
using LiteMath::float3;
using LiteMath::float4;
using LiteMath::Box4f;
using LiteMath::uint;
using cbvh::Interval;
using cbvh::BVHPresets;
using cbvh2::make_BVHNode;

namespace cbvh_embree
{
  cbvh::BVHTree BuildBVH(const float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, cbvh::BVHPresets a_settings);
};


static cbvh::BVHTree BuildSmallBVH(const float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, BVHPresets a_presets)
{
  Box4f box;
  for(size_t i=0;i<a_vertNum;i++)
    box.include(a_vertices[i]); 
    
  cbvh::BVHTree treeData;
  //treeData.indicesReordered = std::vector<uint32_t>(a_indices, a_indices+a_indexNum);
  treeData.indicesReordered = {0,1};
    
  treeData.nodes.resize(a_presets.childrenNum+1);
  treeData.intervals.resize(a_presets.childrenNum+1);
  treeData.nodes[0].boxMin = to_float3(box.boxMin); // - float3(0.1f, 0.1f, 0.1f);
  treeData.nodes[0].boxMax = to_float3(box.boxMax); // + float3(0.1f, 0.1f, 0.1f);;
  treeData.nodes[0].leftOffset  = 1;
  treeData.nodes[0].escapeIndex = uint32_t(-2);
  treeData.intervals[0] = Interval(0, uint32_t(a_indexNum/3));
    
  treeData.nodes[1].boxMin = to_float3(box.boxMin);
  treeData.nodes[1].boxMax = to_float3(box.boxMax);
  treeData.nodes[1].leftOffset  = uint32_t(-1);
  treeData.nodes[1].escapeIndex = (a_presets.desiredFormat == cbvh::FMT_BVH2Node32_Interval32_Dynamic) ? uint32_t(-1) : 2;
  treeData.intervals[1]         = Interval(0, uint32_t(a_indexNum/3));
  for(int i=2;i<a_presets.childrenNum+1;i++)
  {
    treeData.nodes[i].boxMin      = float3(0.0f);
    treeData.nodes[i].boxMax      = float3(0.0f);
    treeData.nodes[i].leftOffset  = uint32_t(-3); // invalid node
    treeData.nodes[i].escapeIndex = (a_presets.desiredFormat == cbvh::FMT_BVH2Node32_Interval32_Dynamic) ? uint32_t(-1) : i+1;
    treeData.intervals[i]         = Interval(0, 0);
  }
  if(a_presets.desiredFormat == cbvh::FMT_BVH4Node32_Interval32_Static || a_presets.desiredFormat == cbvh::FMT_BVH2Node32_Interval32_Static)
    treeData.nodes[treeData.nodes.size()-1].escapeIndex = uint32_t(-2);
  //treeData.depthRanges
  treeData.format = a_presets.desiredFormat;
  return treeData;
}

cbvh::BVHTree cbvh::BuildBVH(const float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, BVHPresets a_presets)
{
  #ifdef CBVH_USE_EMBREE
  if(a_presets.useEmbreeIfPossiable)
  {
    assert(a_presets.childrenNum == 4);
    cbvh::Init();
    auto bvhTree = cbvh::BuildBVH4(a_vertices, a_vertNum, a_indices, a_indexNum);
    bvhTree.format = FMT_BVH4Node32_Interval32_Static;
    cbvh::Destroy();
    return bvhTree;
  }
  #endif
  if(a_presets.btype == cbvh::BVH_CONSTRUCT_FAST || a_presets.btype == cbvh::BVH_CONSTRUCT_FAST_GPU)
    return cbvh::BuildLBVH(a_vertices, a_vertNum, a_indices, a_indexNum, a_presets.desiredFormat, (a_presets.btype == cbvh::BVH_CONSTRUCT_FAST_GPU)); 
  else if(a_indexNum/3 <= size_t(a_presets.primsInLeaf))
    return BuildSmallBVH(a_vertices, a_vertNum, a_indices, a_indexNum, a_presets);
  else if((a_presets.btype == cbvh::BVH_CONSTRUCT_EMBREE2 || a_presets.btype == cbvh::BVH_CONSTRUCT_EMBREE_FAST) && a_presets.childrenNum == 2)
    return cbvh_embree::BuildBVH(a_vertices, a_vertNum, a_indices, a_indexNum, a_presets);
  else
  {
    cbvh_internal::BVHSettings settings;
    settings.enableESC      = (a_presets.btype == BVH_CONSTRUCT_QUALITY);
    settings.primsInLeaf    = a_presets.primsInLeaf;
    settings.childrenNum    = a_presets.childrenNum;
    settings.desiredFormat  = a_presets.desiredFormat;
    settings.alwaysNotEmpty = a_presets.alwaysNotEmpty;
    settings.enableSS       = a_presets.enableSpatialSplit || (a_presets.btype == cbvh::BVH_CONSTRUCT_QUALITY);
    settings.qualityLevel   = 0;
    if(a_presets.btype == cbvh::BVH_CONSTRUCT_MEDIUM)
      settings.qualityLevel = 1;
    else if (a_presets.btype == cbvh::BVH_CONSTRUCT_QUALITY)
      settings.qualityLevel = 2;
    return cbvh_internal::BuildBVH4(a_vertices, a_vertNum, a_indices, a_indexNum, settings);
  }
}

cbvh::BVHTree cbvh::BuildBVH(const BVHNode* a_nodes, size_t a_objNum, BVHPresets a_presets)
{
  #ifdef CBVH_USE_EMBREE
  if(a_presets.useEmbreeIfPossiable && a_presets.childrenNum == 4)
  {
    cbvh::Init();
    auto bvhTree = cbvh::BuildBVH4(a_nodes, a_objNum);
    bvhTree.format = FMT_BVH4Node32_Interval32_Static;
    cbvh::Destroy();
    return bvhTree;
  }
  #endif
  if(a_presets.btype == cbvh::BVH_CONSTRUCT_FAST) // || a_presets.desiredFormat == cbvh::FMT_BVH2Node32_Interval32_Dynamic
  {
    return cbvh::BuildLBVH((const float4*)a_nodes, a_objNum, nullptr, 0, a_presets.desiredFormat);
  }
  else
  {
    cbvh_internal::BVHSettings settings;
    settings.enableESC     = false;
    settings.primsInLeaf   = a_presets.primsInLeaf;
    settings.childrenNum   = a_presets.childrenNum;
    settings.desiredFormat = a_presets.desiredFormat;
    return cbvh_internal::BuildBVH4((const cbvh_internal::Box4f*)a_nodes, a_objNum, nullptr, 0, settings);
  }
}  

void cbvh::BVHTree::ComputeDepthRanges()
{
  if(format != cbvh::CBVH_FORMATS::FMT_BVH4Node32_Interval32_Static) // other formats are not supported currently
    return;

  depthRanges.clear();
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

    depthRanges.push_back(cbvh::Interval(start, end-start+1));

    for(auto x : nodesOfCurrentLevel)
    {
      const auto& node = nodes[x];
      if(cbvh::IsValid(node))
        for(int i=0;i<4;i++)
          treeQueue.push(node.leftOffset + i);
    }
    nodesOfCurrentLevel.clear();
  }

}


static void PrintNode(uint32_t currNodeOffset, uint32_t currDeep, const cbvh::BVHTree& a_tree, std::ostream& out, std::unordered_set<uint32_t>& a_visited)
{
  auto p = a_visited.find(currNodeOffset);
  if(p == a_visited.end())
    a_visited.insert(currNodeOffset);
  else
  {
    std::cout << "node at " << currNodeOffset << " visited twice!" << std::endl;
    out << "node at " << currNodeOffset << " visited twice!" << std::endl;
    std::cout.flush();
    out.flush();
    return;
  }

  const auto currNode = a_tree.nodes[currNodeOffset];
  const bool isLeaf   = (currNode.leftOffset == cbvh::LEAF_NORMAL || currNode.leftOffset == cbvh::LEAF_EMPTY);
  for(uint32_t i=0;i<currDeep;i++)
    out << " ";
  out << " offset = " << currNodeOffset << "; ";
  if(isLeaf)
    out << " leaf : [" << a_tree.intervals[currNodeOffset].start << ":" << a_tree.intervals[currNodeOffset].count << "]";
  else
    out << " node : [" << a_tree.intervals[currNodeOffset].start << ":" << a_tree.intervals[currNodeOffset].count << "]";
  out << std::endl;

  for(uint32_t i=0;i<currDeep;i++)
    out << " ";
  out << " bbox   = {(" << currNode.boxMin.x << ", " << currNode.boxMin.y << ", " << currNode.boxMin.z << ") --";
  out << " (" << currNode.boxMax.x << ", " << currNode.boxMax.y << ", " << currNode.boxMax.z << ")}" << std::endl;

  if(isLeaf)
    return;

  if(a_tree.format == cbvh::CBVH_FORMATS::FMT_BVH2Node32_Interval32_Static)
  {
    PrintNode(currNode.leftOffset + 0, currDeep + 1, a_tree, out, a_visited);
    PrintNode(currNode.leftOffset + 1, currDeep + 1, a_tree, out, a_visited);
  }
  else if(a_tree.format == cbvh::CBVH_FORMATS::FMT_BVH2Node32_Interval32_Dynamic)
  {
    PrintNode(currNode.leftOffset , currDeep + 1, a_tree, out, a_visited);
    PrintNode(currNode.escapeIndex, currDeep + 1, a_tree, out, a_visited);
  }
  else if(a_tree.format == cbvh::CBVH_FORMATS::FMT_BVH4Node32_Interval32_Static)
  {
    PrintNode(currNode.leftOffset + 0, currDeep + 1, a_tree, out, a_visited);
    PrintNode(currNode.leftOffset + 1, currDeep + 1, a_tree, out, a_visited);
    PrintNode(currNode.leftOffset + 2, currDeep + 1, a_tree, out, a_visited);
    PrintNode(currNode.leftOffset + 3, currDeep + 1, a_tree, out, a_visited);
  }
  else // unsupported format
  {

  }
  
}

void cbvh::BVHTree::Print(std::ostream& out)
{
  std::unordered_set<uint32_t> visitedNodes;
  visitedNodes.reserve(nodes.size());
  visitedNodes.clear();

  auto oldPrec = out.precision(3);
  PrintNode(0,0,*this,out,visitedNodes);
  out.precision(oldPrec);
}


struct TempNode
{
  TempNode() : oldOffset(0), escapeIndex(-2), currentDepth(0), sortIndex(0) {} 
  TempNode(uint32_t in_offset, uint32_t in_escapeIndex, uint32_t in_currentDepth) : oldOffset(in_offset), escapeIndex(in_escapeIndex), currentDepth(in_currentDepth), sortIndex(0){}
    
  uint32_t oldOffset;    //<! old offset of the node in nodes array
  uint32_t escapeIndex;  //<!
  uint32_t currentDepth; //<!
  uint32_t sortIndex;    //<! 
};

struct less_than_key
{
  inline bool operator() (const TempNode& struct1, const TempNode& struct2)
  {
    return (struct1.sortIndex < struct2.sortIndex);
  }
};

void ForceGet4ChildernForDynamicLRFormat(const TempNode& node, const cbvh::BVHTree& a_tree, 
                                         TempNode children[4])
{
  const auto& currNode = a_tree.nodes[node.oldOffset];

  const uint32_t leftOffset  = currNode.leftOffset;
  const uint32_t rightOffset = currNode.escapeIndex;

  if(leftOffset != 0xFFFFFFFF && rightOffset != 0xFFFFFFFF)
  {
    const auto& leftNode  = a_tree.nodes[leftOffset];
    const auto& rightNode = a_tree.nodes[rightOffset];

    if(leftNode.leftOffset == 0xFFFFFFFF)
    {
      assert(leftNode.escapeIndex == 0xFFFFFFFF);      
      if(rightNode.leftOffset == 0xFFFFFFFF)
      {
        assert(rightNode.escapeIndex == 0xFFFFFFFF);
        children[0].oldOffset = leftOffset;            // rightLeftOffset
        children[1].oldOffset = rightOffset;
        children[2].oldOffset = 0xFFFFFFFD;            // empty node
        children[3].oldOffset = 0xFFFFFFFD;            // empty node 
      }
      else
      {
        children[0].oldOffset = leftOffset;            // rightLeftOffset
        children[1].oldOffset = rightNode.leftOffset;  // rightLeftOffset
        children[2].oldOffset = rightNode.escapeIndex; // rightRightOffset   
        children[3].oldOffset = 0xFFFFFFFD;
      }
    }
    else if(rightNode.leftOffset == 0xFFFFFFFF)
    {
      assert(rightNode.escapeIndex == 0xFFFFFFFF);
      if(leftNode.leftOffset == 0xFFFFFFFF)
      {
        assert(leftNode.escapeIndex == 0xFFFFFFFF);
        children[0].oldOffset = leftOffset;
        children[1].oldOffset = rightOffset;          // rightLeftOffset
        children[2].oldOffset = 0xFFFFFFFD;           // empty node
        children[3].oldOffset = 0xFFFFFFFD;           // empty node  
      }
      else
      {
        children[0].oldOffset = rightOffset;          // rightLeftOffset
        children[1].oldOffset = leftNode.leftOffset;  // rightLeftOffset
        children[2].oldOffset = leftNode.escapeIndex; // rightRightOffset 
        children[3].oldOffset = 0xFFFFFFFD; 
      } 
    }
    else
    {
      children[0].oldOffset = leftNode.leftOffset;   // leftLeftOffset
      children[1].oldOffset = leftNode.escapeIndex;  // leftRightOffset
      children[2].oldOffset = rightNode.leftOffset;  // rightLeftOffset
      children[3].oldOffset = rightNode.escapeIndex; // rightRightOffset 
    }  
  }
  else if(leftOffset != 0xFFFFFFFF)
  {
    const auto& leftNode  = a_tree.nodes[leftOffset];
    assert(leftNode.leftOffset != 0xFFFFFFFF && leftNode.escapeIndex != 0xFFFFFFFF);
    children[0].oldOffset = leftNode.leftOffset;   // leftLeftOffset
    children[1].oldOffset = leftNode.escapeIndex;  // leftRightOffset
    children[2].oldOffset = 0xFFFFFFFD;            // empty node
    children[3].oldOffset = 0xFFFFFFFD;            // empty node 
  }
  else if(rightOffset != 0xFFFFFFFF)
  {
    const auto& rightNode = a_tree.nodes[rightOffset];
    assert(rightNode.leftOffset != 0xFFFFFFFF && rightNode.escapeIndex != 0xFFFFFFFF);
    children[0].oldOffset = rightNode.leftOffset;  // rightLeftOffset
    children[1].oldOffset = rightNode.escapeIndex; // rightRightOffset  
    children[2].oldOffset = 0xFFFFFFFD;            // empty node
    children[3].oldOffset = 0xFFFFFFFD;            // empty node  
  }
  else
  {
    assert(false);                      // should never happend
    children[0].oldOffset = 0xFFFFFFFD; // rightLeftOffset
    children[1].oldOffset = 0xFFFFFFFD; // rightRightOffset  
    children[2].oldOffset = 0xFFFFFFFD; // empty node
    children[3].oldOffset = 0xFFFFFFFD; // empty node  
  }

  // sort nodes by their intervals.start be sequential!

  for(int i=0;i<4;i++)
  {
    if(children[i].oldOffset < 0xFFFFFFFD)
      children[i].sortIndex = a_tree.intervals[children[i].oldOffset].start;
    else
      children[i].sortIndex = 0xFFFFFFFF;
  }
  
  std::sort(children, children+4, less_than_key());
}

cbvh::BVHTree cbvh::ConvertBVH2DynamicToBVH4Flat(const BVHTree& a_tree)
{
  assert(a_tree.nodes.size() != 0);
  assert(a_tree.format == cbvh::FMT_BVH2Node32_Interval32_Dynamic);

  cbvh::BVHTree res;
  res.nodes.reserve(a_tree.nodes.size());
  res.intervals.reserve(a_tree.intervals.size());
  res.indicesReordered = a_tree.indicesReordered;
  res.format = FMT_BVH4Node32_Interval32_Static;

  TempNode theRootTmp(0, uint32_t(-2), 0);
  
  std::queue<TempNode> treeQueue;
  treeQueue.push(theRootTmp);

  uint32_t currentLevelDepth = 0;
  uint32_t currentLevelStartIndex = 0;

  while (treeQueue.size() != 0) 
  {
    const auto& node = treeQueue.front();

    // check if the current node has jumped to the next depth level
    if(node.currentDepth > currentLevelDepth)
    {
      Interval currentLevelDepthInterval(currentLevelStartIndex, res.nodes.size() - currentLevelStartIndex);
      // push the interval of the previous depth level
      res.depthRanges.push_back(currentLevelDepthInterval);
      currentLevelDepth += 1;
      currentLevelStartIndex = res.nodes.size();
    }

    const int escapeIndex = node.escapeIndex;
    const int bufferSize  = res.nodes.size() + treeQueue.size();
    
    //if(node.oldOffset == 23199)
    //{
    //  std::cout << "process target node at offset" << node.oldOffset << std::endl;
    //}

    if(node.oldOffset >= 0xFFFFFFFD) // empty node
    {
      float3 boxMin(+1e+9f,+1e+9f,+1e+9f);
      float3 boxMax(-1e+9f,-1e+9f,-1e+9f);

      res.nodes.push_back(make_BVHNode(boxMin, boxMax, 0xFFFFFFFD, escapeIndex));
      res.intervals.push_back(Interval(0,0));
    }
    else
    {
      const auto& oldNode     = a_tree.nodes    [node.oldOffset];
      const auto& oldInterval = a_tree.intervals[node.oldOffset]; 
      
      const bool isLeaf       = (oldNode.leftOffset == 0xFFFFFFFF);
      const uint32_t nextOffs = isLeaf ? 0xFFFFFFFF : bufferSize;

      res.nodes.push_back(make_BVHNode(oldNode.boxMin, oldNode.boxMax, nextOffs, escapeIndex));
      res.intervals.push_back(oldInterval);

      if(!isLeaf) 
      {
        assert(IsValid(oldNode));
        // force get 4 clildren and put them in to queue
        TempNode children[4];
        ForceGet4ChildernForDynamicLRFormat(node, a_tree, children);
         
        for (auto i = 0; i < 3; i++) 
        {
          children[i].escapeIndex  = bufferSize + i + 1;
          children[i].currentDepth = node.currentDepth + 1;
          treeQueue.push(children[i]);
        }

        children[3].escapeIndex  = escapeIndex;
        children[3].currentDepth = node.currentDepth + 1;
        treeQueue.push(children[3]);
      }
      else
      {
        assert(oldNode.escapeIndex == 0xFFFFFFFF); // escapeIndex is rightOffset for oldNode, check if it is -1, for leaf nodes
      }
    }

    treeQueue.pop();
  };

  Interval currentLevelDepthInterval(currentLevelStartIndex, res.nodes.size() - currentLevelStartIndex);  
  res.depthRanges.push_back(currentLevelDepthInterval); 

  //res.ComputeDepthRanges();
  return res;
}

using cbvh::BVHNode;
using cbvh::Interval;
static constexpr uint32_t LEAF_BIT = 0x80000000;

struct Converter2to2
{
  Converter2to2(const std::vector<BVHNode>&  in_nodes,     std::vector<BVHNode>&  out_nodes,
                const std::vector<Interval>& in_intervals, std::vector<Interval>& out_intervals) : m_inNodes(in_nodes), m_outNodes(out_nodes), m_inIntervals(in_intervals), m_outIntervals(out_intervals)  
  {
    m_outNodes.reserve(m_inNodes.size());
    m_outNodes.resize(0);
  }

  uint32_t Visit(uint32_t a_leftOffset, uint32_t a_rightOffset, uint32_t a_parentEscapeIndex) //
  {
    const BVHNode leftNode  = m_inNodes[a_leftOffset];
    const BVHNode rightNode = m_inNodes[a_rightOffset];

    const uint32_t currSize = m_outNodes.size();
    m_outNodes.resize(currSize + 2);
    m_outIntervals.resize(currSize + 2);

    m_outNodes[currSize + 0] = leftNode;
    m_outNodes[currSize + 1] = rightNode;
    
    m_outIntervals[currSize + 0] = m_inIntervals[a_leftOffset];
    m_outIntervals[currSize + 1] = m_inIntervals[a_rightOffset];
    
    m_outNodes[currSize + 0].escapeIndex = currSize + 1;
    m_outNodes[currSize + 1].escapeIndex = a_parentEscapeIndex;

    if((leftNode.leftOffset & LEAF_BIT) == 0)
      m_outNodes[currSize + 0].leftOffset = Visit(leftNode.leftOffset, leftNode.escapeIndex, currSize + 1);
    
    if((rightNode.leftOffset & LEAF_BIT) == 0)
      m_outNodes[currSize + 1].leftOffset = Visit(rightNode.leftOffset, rightNode.escapeIndex, a_parentEscapeIndex);

    return currSize;
  }

  const std::vector<BVHNode>&  m_inNodes;
  std::vector<BVHNode>&        m_outNodes;

  const std::vector<Interval>& m_inIntervals;
  std::vector<Interval>&       m_outIntervals;
};

cbvh::BVHTree cbvh::ConvertBVH2DynamicToBVH2Flat(const BVHTree& a_tree)
{
  assert(a_tree.nodes.size() != 0);
  assert(a_tree.format == cbvh::FMT_BVH2Node32_Interval32_Dynamic);

  //PrintBVH2ForGraphViz_LR(a_tree.nodes, "w_helix.txt");

  cbvh::BVHTree res;
  res.nodes.reserve(a_tree.nodes.size());
  res.intervals.reserve(a_tree.intervals.size());
  res.indicesReordered = a_tree.indicesReordered;
  res.format = FMT_BVH2Node32_Interval32_Static;

  auto root = a_tree.nodes[0];
  if((root.leftOffset & LEAF_BIT) != 0)
  {
    res.nodes.resize(2);
    res.intervals.resize(2);
    
    res.nodes[0] = root;
    res.nodes[0].leftOffset  = 0xFFFFFFFF;
    res.nodes[0].escapeIndex = uint32_t(-2);
    res.nodes[1] = cbvh2::DummyNode();

    res.intervals[0] = a_tree.intervals[0];
    res.intervals[1] = Interval(0,0);
  }
  else
  {
    Converter2to2 converter(a_tree.nodes, res.nodes, a_tree.intervals, res.intervals);
    converter.Visit(root.leftOffset, root.escapeIndex, uint32_t(-2));
  }

  return res;
}
