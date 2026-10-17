#include <iostream>
#include <queue>
#include <stack>
#include <limits>
#include <cassert>
#include <cmath>

#include <algorithm>
#include <unordered_set>
#include <set>
#include <fstream>
#include <chrono>

#include "cbvh.h"
#include "cbvh_core.h"
#include "cbvh_internal.h"
#include "EXT_TriBoxOverlap.h"

using LiteMath::float3;
using LiteMath::float4;
using LiteMath::float4x4;

using LiteMath::int4;
using LiteMath::uint4;

using LiteMath::min;
using LiteMath::max;
using LiteMath::to_uint32;

using LiteMath::length3f;
using LiteMath::uint;

using cbvh_internal::Box4f;
using cbvh_internal::Triangle4f;
using cbvh2::make_BVHNode;

Box4f CalcBBox(const Box4f* a_boxes, size_t a_objNum)
{
  Box4f res;
  for(size_t i=0;i<a_objNum;i++)
    res.include(a_boxes[i]);
  return res;
}

template<int axis>
struct box_sort_axis
{
  inline bool operator()(const Box4f& struct1, const Box4f& struct2)
  {
    const float boxCenter1 = 0.5f*(struct1.boxMin[axis] + struct1.boxMax[axis]);
    const float boxCenter2 = 0.5f*(struct2.boxMin[axis] + struct2.boxMax[axis]);
    return (boxCenter1 < boxCenter2);
    //return (struct1.boxMin[axis] < struct2.boxMin[axis]);
  }
};

using cbvh_internal::SplitInfo;
using cbvh_internal::BVHSettings;

SplitInfo FindObjectSplitSAH(const std::vector<Box4f>& in_boxes) 
{
  Box4f boxLeft, boxRight;
  std::vector<Box4f> cummBoxRight(in_boxes.size());

  for(int i=in_boxes.size()-1;i>=0;i--)
  {
    boxRight.include(in_boxes[i]);
    cummBoxRight[i] = boxRight;
  }
  
  SplitInfo res;
  res.splitPos  = 0;
  res.metricVal = +std::numeric_limits<float>::infinity();
  
  for(size_t i=0;i<in_boxes.size()-1;i++)
  {
    boxLeft.include(in_boxes[i]);
    
    size_t numLeft  = i+1;
    size_t numRight = in_boxes.size() - i - 1;

    //const Box4f overlap = BoxBoxOverlap(boxLeft, cummBoxRight[i]);
    //const float ovlpSAH = overlap.surfaceArea()*1.0f*float(std::min(numLeft, numRight));
    //const float ovlpSAH = overlap.volume()*1.0f*float(std::min(numLeft, numRight));
  
    //numLeft  = std::max(size_t(1),numLeft - numLeft%4);
    //numRight = std::max(size_t(1),numRight - numRight%4);

    const float leftSAH  = boxLeft.surfaceArea() * float(numLeft); 
    const float rightSAH = cummBoxRight[i].surfaceArea() * float(numRight);
    const float currSAH  = leftSAH + rightSAH;

    if(currSAH < res.metricVal)
    {
      res.metricVal = currSAH;
      res.splitPos  = i;
      res.boxL      = boxLeft;
      res.boxR      = cummBoxRight[i];
    }
  }

  return res;
}

cbvh_internal::SplitInfo cbvh_internal::FindObjectSplit(const std::vector<Box4f>& in_boxes, const Box4f& in_bbox, SPLIT_TYPE in_splitType,
                                                        std::vector<Box4f>& out_boxesLeft, std::vector<Box4f>& out_boxesRight)
{ 
  constexpr int NUMSPLITS = 3;
  std::vector<Box4f> boxesSorted[NUMSPLITS] = {in_boxes,in_boxes,in_boxes}; 

  std::sort(boxesSorted[0].begin(), boxesSorted[0].end(), box_sort_axis<0>());
  std::sort(boxesSorted[1].begin(), boxesSorted[1].end(), box_sort_axis<1>());
  std::sort(boxesSorted[2].begin(), boxesSorted[2].end(), box_sort_axis<2>()); 

  SplitInfo splits[NUMSPLITS];
  for(int i=0;i<NUMSPLITS;i++)
    splits[i] = FindObjectSplitSAH(boxesSorted[i]);
  
  SplitInfo minSplit;
  for(int i=0;i<NUMSPLITS;i++)
  {
    if(splits[i].metricVal < minSplit.metricVal)
    {
      minSplit = splits[i];
      minSplit.splitAxis = i;
    }
  }
  
  auto& boxesSelected = boxesSorted[minSplit.splitAxis];
  const auto splitPos = boxesSelected.begin() + minSplit.splitPos + 1;

  out_boxesLeft  = std::vector<Box4f>(boxesSelected.begin(), splitPos);
  out_boxesRight = std::vector<Box4f>(splitPos, boxesSelected.end());

  return minSplit;
}

struct TmpNode
{
  std::vector<Box4f>  plist;              // #TODO: use fast list here
  Box4f               box;
  int                 escapeIndex;
  int                 currentDepth;
  bool                stopNext = false;
  size_t              primsCountFake = 0; // used for spaial split search

  //inline bool AxisAligned(int axis, float split) const { return (box.boxMin[axis] == box.boxMax[axis]) && (box.boxMin[axis]==split); }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using cbvh_internal::OriginalInput;

/**
\brief Perform Spatial Split for triangle
\param  a_boxOriginal - input bounding box of a triangle (may already cover only part of a triangle!)
\param  verts         - input triangle vertices
\param  splitAxis     - input split axis: 0,1,2
\param  splitPos      - input split pos: from a_node.boxMin[splitAxis] to a_node.boxMax[splitAxis]
\param  leftBox       - output left bounding box
\param  rightBox      - output right bounding box  
*/
void SpatialSplitTriangle(const Box4f& a_boxOriginal, const float4 verts[3], int splitAxis, float splitPos,
                          Box4f& leftBox, Box4f& rightBox)
{
  const float4 edges[3][2] = { {verts[0], verts[1]},
                               {verts[2], verts[0]},
                               {verts[1], verts[2]} };    
  for (int i=0;i<3;i++)
  {
    const float v0p = edges[i][0][splitAxis];
    const float v1p = edges[i][1][splitAxis];
    
    // Insert vertex to the boxes it belongs to.
    //
    if(v0p <= splitPos)
      leftBox.include(edges[i][0]);
    
    if(v0p >= splitPos)
      rightBox.include(edges[i][0]);
    
    // Edge intersects the plane => insert intersection to both boxes.
    //
    if ((v0p < splitPos && v1p > splitPos) || (v0p > splitPos && v1p < splitPos))
    {
      const float  t = LiteMath::clamp((splitPos - v0p) / (v1p - v0p), 0.0f, 1.0f);
      const float4 p = lerp(edges[i][0], edges[i][1], t);
    
      leftBox.include(p);
      rightBox.include(p);
    }
  } // end for (int i=0;i<3;i++)
      
  leftBox.boxMax [splitAxis] = splitPos;
  rightBox.boxMin[splitAxis] = splitPos;
  
  // Intersect with original bounds.
  //
  leftBox.intersect(a_boxOriginal);
  rightBox.intersect(a_boxOriginal);

  leftBox.setStart(a_boxOriginal.getStart());
  leftBox.setCount(1);
  rightBox.setStart(a_boxOriginal.getStart());
  rightBox.setCount(1);
}

/**
\brief Perform Spatial Split for the node. 
\param  a_node    - input node
\param  splitAxis - input split axis: 0,1 or 2
\param  splitPos  - input split pos: from a_node.boxMin[splitAxis] to a_node.boxMax[splitAxis]
\param  a_input   - input access to triangle mesh
\param  onlyCount - input flag for split estimation without actual filling output arrays
\param  a_left    - output left  node.  
\param  a_right   - output right node.   
\return Splitted bounding boxes of triangles. Long bad triangles are represented with several boxes.
*/
inline void PerformSpatialSplit(const TmpNode& nodeBox, int splitAxis, float splitPos, const OriginalInput& a_input, bool onlyCount,
                                TmpNode& leftBox, TmpNode& rightBox)
{
  
  const float4* a_vertices  = a_input.vertices;
  const uint32_t* a_indices = a_input.indices;
 
  leftBox.plist.resize(0);
  rightBox.plist.resize(0);

  if(leftBox.plist.capacity() < nodeBox.plist.size() && !onlyCount)
    leftBox.plist.reserve(nodeBox.plist.size());

  if(rightBox.plist.capacity() < nodeBox.plist.size() && !onlyCount)
    rightBox.plist.reserve(nodeBox.plist.size());
  
  leftBox.box  = Box4f();
  rightBox.box = Box4f();
  for(const auto& triBox : nodeBox.plist) 
  {
    if(triBox.boxMax[splitAxis] <= splitPos) 
    {
      if(!onlyCount)
        leftBox.plist.push_back(triBox);
      leftBox.box.include(triBox);
      leftBox.primsCountFake++;
    }
    else if(triBox.boxMin[splitAxis] >= splitPos) 
    {
      if(!onlyCount)
        rightBox.plist.push_back(triBox);
      rightBox.box.include(triBox);
      rightBox.primsCountFake++;
    }
    else
    {
      const uint32_t triId  = triBox.getStart();
      const float4 verts[3] = {a_vertices[a_indices[triId*3+0]], 
                               a_vertices[a_indices[triId*3+1]],
                               a_vertices[a_indices[triId*3+2]]};
      Box4f lBox, rBox;                         
      SpatialSplitTriangle(triBox, verts, splitAxis, splitPos,
                           lBox, rBox);
      
      if(!onlyCount)
      {
        leftBox.plist.push_back(lBox);
        rightBox.plist.push_back(rBox);
      }
      leftBox.box.include(lBox);  
      rightBox.box.include(rBox);        
      leftBox.primsCountFake++;
      rightBox.primsCountFake++;    
    }
  }
}


struct SpatialSplitInfo 
{
  SpatialSplitInfo(){}
  SpatialSplitInfo(int a_axis, float a_split) : axis(a_axis), pos(a_split) {} 
  int axis;
  float pos;
};


SpatialSplitInfo FindBestSpatialSplit(const Box4f& leftSweepBox, const Box4f& rightSweepBox, const TmpNode& nodeBox, const OriginalInput& a_input)
{
  const Box4f ovlp = BoxBoxOverlap(leftSweepBox, rightSweepBox);
  const Box4f nodb = Box4f( min(leftSweepBox.boxMin, rightSweepBox.boxMin), max(leftSweepBox.boxMax, rightSweepBox.boxMax));
  
  const float4 center1 = 0.5f*(ovlp.boxMin + ovlp.boxMax);
  const float4 center2 = 0.5f*(nodb.boxMin + nodb.boxMax);

  // (1) propose candidate splits
  //
  std::vector<SpatialSplitInfo> candidates;
  candidates.reserve(6);
   
  for(int axis = 0; axis < 3; axis++)
  {
    candidates.push_back(SpatialSplitInfo(axis,center1[axis]));
    candidates.push_back(SpatialSplitInfo(axis,center2[axis]));
  }
  
  // (2) find best candidate
  //
  float minSAH = std::numeric_limits<float>::max();
  int minId    = 0;

  for(int i=0;i<int(candidates.size());i++)
  {
    TmpNode leftBox, rightBox;
    PerformSpatialSplit(nodeBox, candidates[i].axis, candidates[i].pos, a_input, true,
                        leftBox, rightBox);

    const float currSAH = leftBox.box.surfaceArea()*float(leftBox.primsCountFake) + rightBox.box.surfaceArea()*float(rightBox.primsCountFake);
    if(currSAH < minSAH)
    {
      minSAH = currSAH;
      minId  = i;
    }
  }

  return candidates[minId];
}


static inline std::set<uint32_t> JoinPrimitives(const std::vector<Box4f>& plist)
{
  std::set<uint32_t> trisIds;
  for(size_t i=0;i<plist.size();i++)
  {
    const uint firstPrim   = plist[i].getStart();
    const uint primsInLeaf = plist[i].getCount();
    for(uint j=firstPrim; j < firstPrim + primsInLeaf;j++)
      trisIds.insert(j);
  }
  return trisIds;
}

cbvh::Interval JoinAllPrimitives(const std::vector<Box4f>& plist, std::vector<uint>& primIdBFS, uint32_t a_indexNum, bool a_dups)
{
  const uint primsStart = uint(primIdBFS.size());
  //if(plist.size() == 1) {
  //  primIdBFS.push_back(j);
  //  return cbvh::Interval(plist[0].getStart(), plist[0].getCount());
  //}
    
  if(a_dups)
  {
    auto trisIds = JoinPrimitives(plist);
    for(auto triIndex : trisIds)
      primIdBFS.push_back(triIndex);
  }
  else
  {
    for(size_t i=0;i<plist.size();i++)
    {
      const uint firstPrim   = plist[i].getStart();
      const uint primsInLeaf = plist[i].getCount();
      for(uint j=firstPrim; j < firstPrim + primsInLeaf;j++)
      {
        assert(j < a_indexNum);
        primIdBFS.push_back(j);
      }
    }
  }

  const uint primsEnd   = uint(primIdBFS.size());
  const uint primsCount = primsEnd - primsStart;

  return cbvh::Interval(primsStart, primsCount);
}

static inline size_t CountPrims(const TmpNode& node, const BVHSettings& a_settings)
{
  if(!a_settings.enableESC || node.plist.size() > cbvh_internal::SPLIT_ANYWAY_THRESHOLD)
    return node.plist.size();
  else
  {
    const auto primIds = JoinPrimitives(node.plist);
    return primIds.size();
  }
}

static inline bool StopByPrimsNum(const TmpNode& node, const BVHSettings& a_settings)
{
  if(!a_settings.enableESC)
    return node.plist.size() < size_t(a_settings.primsInLeaf);
  else if (node.plist.size() > cbvh_internal::SPLIT_ANYWAY_THRESHOLD)
    return false;
  else
  {
    const auto primsNum = CountPrims(node, a_settings);
    return (primsNum <= size_t(a_settings.primsInLeaf));
  }
}

static inline int SelectMaxAxis(float4 a_boxSize)
{
  int axis = 0;
  float maxVal = 0.0f;
  for(int i=0;i<3;i++)
  {
    if(a_boxSize[i] > maxVal)
    {
      axis   = i;
      maxVal = a_boxSize[i];
    }
  }
  return axis;
}

static uint32_t g_spatialSplits=0;
static size_t   g_totalPrims=0;

struct ThresholdValuesSS
{
  uint32_t totalPrims;
  uint32_t totalSplits;
  uint32_t stopPrims;
};

ThresholdValuesSS g_thresholds[4] = { {32,16,2}, {64,32,2}, {256,64,4}, {1024,256,4} };

static inline bool SpatialSplitCriterion(int a_currDepth, float a_overlap, size_t a_primCount, int a_qualityLevel)
{ 
  size_t normalPrimsCount = 256;
  //if(g_totalPrims <= 1024)
  //{
  //  for(int i=0;i<4;i++) {
  //    if(g_totalPrims <= g_thresholds[i].totalPrims) {
  //      normalPrimsCount = size_t(g_thresholds[i].stopPrims);
  //      break;
  //    }
  //  }
  //}

  if(a_qualityLevel >= 2)
  {
    if(a_primCount < normalPrimsCount)
      return false;
    if(a_currDepth <= 2 && a_overlap >= 0.10f)
      return true;   
    if(a_currDepth <= 4 && a_overlap >= 0.15f)
      return true;
    if(a_currDepth <= 6 && a_overlap >= 0.20f)
      return true;
    if(a_overlap >= 0.25f)
      return true;
  }
  else
  {
    if(a_primCount < normalPrimsCount*2)
      return false;
    if(a_currDepth <= 2 && a_overlap >= 0.20f)
      return true;   
    if(a_currDepth <= 4 && a_overlap >= 0.30f)
      return true;
    if(a_currDepth <= 6 && a_overlap >= 0.40f)
      return true;
    if(a_overlap >= 0.50f)
      return true;
  }

  return false;
}

static inline bool EnoughSpatialSplits()
{
  if(g_totalPrims <= 1024)
  {
    uint32_t enoughSplits = 10;
    for(int i=0;i<4;i++) {
      if(g_totalPrims <= g_thresholds[i].totalPrims) {
        enoughSplits = size_t(g_thresholds[i].stopPrims);
        break;
      }
    }
  
    return (g_spatialSplits > enoughSplits);
  }

  return (g_totalPrims > 1000000) && (g_spatialSplits > g_totalPrims/100);
}

bool Subdivide2(const TmpNode& a_node, TmpNode& a_left, TmpNode& a_right, const BVHSettings& a_settings, const OriginalInput& a_input)
{
  if(a_node.plist.size() <= 1)
    return false;
    
  auto split = cbvh_internal::FindObjectSplit(a_node.plist, a_node.box, cbvh_internal::OBJECT_SPLIT_SAH,
                                              a_left.plist, a_right.plist);

  a_left.box           = split.boxL;
  a_right.box          = split.boxR;
  a_left.currentDepth  = a_node.currentDepth+1;
  a_right.currentDepth = a_node.currentDepth+1;

  const float overlap  = cbvh_internal::RelativeBoxOverlap(a_left.box, a_right.box);
  const bool  hqSSplit = (a_settings.qualityLevel >= 2); 
  const bool  enableSS = (a_settings.qualityLevel >= 1 && a_settings.enableSS) || (a_settings.qualityLevel >= 2); 
  const bool  decision = SpatialSplitCriterion(a_node.currentDepth, overlap, a_node.plist.size(), a_settings.qualityLevel);
  const bool  enouhSS  = EnoughSpatialSplits();

  if(enableSS && decision && !enouhSS && a_input.indices != nullptr)
  {
    //std::cout << "before `PerformSpatialSplit`, plist.size() = " << a_node.plist.size() << std::endl;
    g_spatialSplits++;
    if(hqSSplit && a_node.plist.size() < 1000000)
    {
      auto splitInfo = FindBestSpatialSplit(a_left.box, a_right.box, a_node, a_input);
      PerformSpatialSplit(a_node, splitInfo.axis, splitInfo.pos, a_input, false,
                          a_left, a_right);
    }
    else
    {
      const int axis  = SelectMaxAxis(a_node.box.boxMax - a_node.box.boxMin);
      const float pos = 0.5f*(a_node.box.boxMin + a_node.box.boxMax)[axis];
      PerformSpatialSplit(a_node, axis, pos, a_input, false,
                          a_left, a_right);
    }
    //std::cout << "after  `PerformSpatialSplit`, plist.size() = " << a_node.plist.size() << std::endl;
    return true;
  }
  
  if(!a_settings.esc.prebuildForESC)
  {
    if(a_node.plist.size() > cbvh_internal::SPLIT_ANYWAY_THRESHOLD)
      return true;
    
    const size_t primNum  = CountPrims(a_node, a_settings);
    const float overlap2  = cbvh_internal::RelativeBoxOverlap(a_left.box, a_right.box);
    return (overlap2 < 0.75f) || (primNum > size_t(a_settings.primsInLeaf*2));
  }
  else
  {
    if(a_node.plist.size() > cbvh_internal::SPLIT_ANYWAY_THRESHOLD)
      return true;
    
    const float overlap2 = cbvh_internal::RelativeBoxOverlap(a_left.box, a_right.box);
    if(overlap2 >= 0.5f)
      return false;
    else
      return true; //splitSAH < nodeSAH; // always true actually ... 
  }
}

bool Subdivide3(const TmpNode& a_node, TmpNode* a_children, const BVHSettings& a_settings, const OriginalInput& a_input)
{
  TmpNode tmp3, tmp4;
  bool res = Subdivide2(a_node, tmp3, tmp4, a_settings, a_input);

  if(tmp3.plist.size() <= 1)
  {
    a_children[0] = tmp3;
    Subdivide2(tmp4, a_children[1], a_children[2], a_settings, a_input);
  }
  else 
  {
    Subdivide2(tmp3, a_children[0], a_children[1], a_settings, a_input);
    a_children[2] = tmp4;
  }

  return res;
}

bool Subdivide(const TmpNode& a_node, TmpNode* a_children, int a_childrenNum, const BVHSettings& a_settings, const OriginalInput& a_input)
{
  assert(a_childrenNum == 4 || a_childrenNum == 2);  // #TODO: make it working for any number of 'a_childrenNum'

  for(int i=0;i<a_childrenNum;i++)
    a_children[i] = TmpNode();

  if(a_node.plist.size() <= size_t(a_childrenNum)) // too small number of primitives, need special case
  {
    for(size_t i=0;i<a_node.plist.size();i++)
    {
      auto leafBox = a_node.plist[i];
      auto start   = leafBox.getStart();
      auto count   = leafBox.getCount();
      leafBox.intersect(a_node.box);  // due to spatial split triangle boxes could be much larger than a_node.box
      leafBox.setStart(start);
      leafBox.setCount(count);
      a_children[i].box = leafBox;
      a_children[i].plist.push_back(a_node.plist[i]);
    }

    //for(auto i=a_node.plist.size();i<a_childrenNum;i++) //empty nodes
    //  a_children[i].leftOffset = 0xFFFFFFFD;   

    return true;
  }

  if(a_childrenNum == 2)
    return Subdivide2(a_node, a_children[0], a_children[1], a_settings, a_input);

  TmpNode tmp1, tmp2;
  bool splitOk = Subdivide2(a_node, tmp1, tmp2, a_settings, a_input);

  if(tmp1.plist.size() <= 1)
  {
    a_children[0] = tmp1;
    Subdivide3(tmp2, a_children+1, a_settings, a_input);
  }
  else if(tmp2.plist.size() <= 1)
  {
    a_children[3] = tmp2;
    Subdivide3(tmp1, a_children+0, a_settings, a_input);
  }
  else
  {
    Subdivide2(tmp1, a_children[0], a_children[1], a_settings, a_input);
    Subdivide2(tmp2, a_children[2], a_children[3], a_settings, a_input);
  }

  int allowedDeviation = a_settings.primsInLeaf/2;
  bool allChildsAreOk  = !a_settings.esc.prebuildForESC;

  for(auto i=0;i<a_childrenNum;i++)
  {
    a_children[i].box.intersect(a_node.box);  // due to spacial split triangle boxes could be much larger than a_node.box
    allChildsAreOk = allChildsAreOk && (a_children[i].plist.size() <= size_t(a_settings.primsInLeaf + allowedDeviation));
  }

  //if(allChildsAreOk) // apply heuristic stop criteria
  //{
  //  for(auto i=0;i<a_childrenNum;i++)
  //    a_children[i].stopNext = true;
  //}

  return splitOk;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct DFSParamsInternal 
{
  DFSParamsInternal(const std::vector<cbvh::Interval>& a_primsBE, 
                    const std::vector<uint>&           a_primIdBFS,
                    cbvh::BVHTree&                     a_bvhOurs,
                    const size_t                       a_ibSize,
                    const int                          a_childrenNum) : m_primsBegEnd(a_primsBE), 
                                                                        m_primIdBFS(a_primIdBFS),
                                                                        m_bvhOurs(a_bvhOurs),
                                                                        m_indexBufferSize(a_ibSize), m_numChildren(a_childrenNum) {}


  size_t BFSTraversalRec(size_t currOffs, int childId = 0) 
  {
    if(currOffs >= m_bvhOurs.nodes.size()) 
    {
      int a = 2;
    }
    else if(m_bvhOurs.nodes.size() == 871)
      std::cout << "currOffs = " << currOffs << " of child(" << childId << ")" << std::endl;

    const cbvh::BVHNode& node = m_bvhOurs.nodes[currOffs];
    if (!cbvh::IsLeaf(node)) 
    {
      // write down intervals, basically taking the sum from the bottom levels of tree
      //
      for (int i = 0; i < m_numChildren; i++)
      {
        const size_t start = m_bvhOurs.indicesReordered.size();
        const size_t end   = BFSTraversalRec(node.leftOffset + i, i);
        
        m_bvhOurs.intervals[node.leftOffset + i].start = uint(start)      ;
        m_bvhOurs.intervals[node.leftOffset + i].count = uint(end - start);
      }
    } 
    else if (!cbvh::IsEmpty(node)) 
    {
      const cbvh::Interval begEnd = m_primsBegEnd[currOffs];
      for(uint i = begEnd.start; i < begEnd.start + begEnd.count; i++)
        m_bvhOurs.indicesReordered.push_back(m_primIdBFS[i]); 
    } 
    
    return m_bvhOurs.indicesReordered.size();
  }

  const std::vector<cbvh::Interval>& m_primsBegEnd;
  const std::vector<uint>&           m_primIdBFS;
  cbvh::BVHTree&                     m_bvhOurs;
  size_t                             m_indexBufferSize;
  int                                m_numChildren;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern double g_buildTime;

cbvh::BVHTree cbvh_internal::BuildBVH4(const Box4f* a_boxes, size_t a_objNum, const uint* a_indices, size_t a_indexNum, BVHSettings a_settings, OriginalInput a_args)
{
  g_spatialSplits = 0;
  g_totalPrims    = a_objNum;
  auto before     = std::chrono::high_resolution_clock::now();
  
  cbvh::BVHTree res;
  
  if(a_indices == nullptr)
  {
    a_indexNum = a_objNum;      // don't use indices, actual object index is stored inside Box4f
    a_settings.primsInLeaf = 1; // when build bvh for custom objects we force it to single object in leaf

    if(a_objNum == 1)
    {
      res.nodes.resize(1);
      res.intervals.resize(1);
      res.indicesReordered.resize(1);

      res.nodes[0].boxMin      = to_float3(a_boxes[0].boxMin);
      res.nodes[0].boxMax      = to_float3(a_boxes[0].boxMax);
      res.nodes[0].leftOffset  = 0x80000000;
      res.nodes[0].escapeIndex = cbvh::ESCAPE_ROOT;
      res.intervals[0]         = cbvh::Interval(0,1);
      res.indicesReordered[0]  = 0;
      
      return res;
    }
  }

  const size_t nodesApproxNum = a_objNum * 2 + 10;
  res.nodes.reserve(nodesApproxNum);
  res.intervals.reserve(nodesApproxNum);
  res.indicesReordered.reserve(a_objNum);
  res.format = cbvh::FMT_BVH4Node32_Interval32_Static;

  std::vector<uint> primIdBFS;
  std::vector<cbvh::Interval> primsBegEnd;
  primsBegEnd.reserve(nodesApproxNum);
  primIdBFS.reserve(nodesApproxNum);

  //cvex4::set_ftz(); ///
  const  size_t progressGranularity = a_objNum/100; 
  //size_t primsDone = 0;
  //size_t primsDoneOld = 0;
  

  // BFS Build
  //
  int currentLevelDepth      = 0;
  int currentLevelStartIndex = 0;
  
  TmpNode theRootTmp;
  theRootTmp.box          = CalcBBox(a_boxes, a_objNum); // #TODO: pass this via arguments
  theRootTmp.plist        = std::vector<Box4f>(a_boxes, a_boxes + a_objNum);
  theRootTmp.escapeIndex  = uint(-2); 
  theRootTmp.currentDepth = 0;
  
  std::queue<TmpNode> treeQueue;
  treeQueue.push(theRootTmp);
  
  while (treeQueue.size() != 0) 
  {
    // check if the current node has jumped to the next depth level
    if(treeQueue.front().currentDepth > currentLevelDepth)
    {
      Interval currentLevelDepthInterval(currentLevelStartIndex, res.nodes.size() - currentLevelStartIndex);
      // push the interval of the previous depth level
      res.depthRanges.push_back(currentLevelDepthInterval);
      currentLevelDepth += 1;
      currentLevelStartIndex = res.nodes.size();
    }

    const auto& node        = treeQueue.front();
    const int escapeIndex   = node.escapeIndex;
    const size_t bufferSize = res.nodes.size() + treeQueue.size();
    
    // (1) Subdivide primitives array to a target children number (4 by default) 
    //
    std::vector<TmpNode> children(a_settings.childrenNum);
    bool subdivOk = false;
    if(node.plist.size() > size_t(a_settings.primsInLeaf))
      subdivOk = Subdivide(node, children.data(), a_settings.childrenNum, a_settings, a_args);

    const bool stopNow = StopByPrimsNum(node, a_settings);

    if(!stopNow && subdivOk && !node.stopNext)
    {
      const uint32_t leftOffset = node.plist.size() > 0 ? bufferSize : 0xFFFFFFFD;

      res.nodes.push_back(make_BVHNode(to_float3(node.box.boxMin), to_float3(node.box.boxMax), leftOffset, escapeIndex));
      primsBegEnd.push_back(cbvh::Interval());

      // (2) process all children via same queue 
      //
      const auto lastChild = a_settings.childrenNum-1;
      for (auto i = 0; i < lastChild; i++) 
      {
        children[i].escapeIndex  = bufferSize + i + 1;
        children[i].currentDepth = node.currentDepth + 1;
        treeQueue.push(children[i]);
      }
      
      children[lastChild].escapeIndex  = escapeIndex;
      children[lastChild].currentDepth = node.currentDepth + 1;
      treeQueue.push(children[lastChild]);

      treeQueue.pop();
    } 
    else // make leaf
    {
      auto interval = JoinAllPrimitives(node.plist, primIdBFS, a_indexNum, a_settings.enableESC);
     
      //const uint32_t leftOffset = interval.count > 0 ? 0xFFFFFFFF : 0xFFFFFFFD;
      
      //primsDone += interval.count;
      //if(a_objNum > 1000000 && primsDone > primsDoneOld + progressGranularity)
      //{
      //  if(currentLevelDepth == 217 && treeQueue.size() == 1724)
      //  {
      //    int a = 2;
      //  }
      //  primsDoneOld = primsDone;
      //  std::cout << "progress = " << float(primsDone)/float(a_objNum) << "%     \r"; //<< "%; (depth, qsize) = (" << currentLevelDepth << ", " << treeQueue.size() << ")" << std::endl; // << "%     \r";
      //  std::cout.flush();
      //}

      res.nodes.push_back(make_BVHNode(to_float3(node.box.boxMin), to_float3(node.box.boxMax), 0xFFFFFFFF, escapeIndex));

      primsBegEnd.push_back(interval); 
      treeQueue.pop();
    }
  }

  // push the interval of the previous depth level
  //
  Interval currentLevelDepthInterval(currentLevelStartIndex, res.nodes.size() - currentLevelStartIndex);  
  res.depthRanges.push_back(currentLevelDepthInterval); 

  // DFS pass (Second pass with depth-first traversal to form triangle indices)
  //
  DFSParamsInternal dfsParams(primsBegEnd, primIdBFS, res, a_indexNum, a_settings.childrenNum);
  res.intervals.resize(primsBegEnd.size());
  dfsParams.BFSTraversalRec(0);
  res.intervals[0] = cbvh::Interval(0, uint(primIdBFS.size()));
  
  if(a_settings.childrenNum == 2 && a_settings.desiredFormat == cbvh::FMT_BVH2Node32_Interval32_Dynamic) 
  {
    // replace escapeIndex with rightOffset for 'FMT_BVH2Node32_Interval32_Dynamic' format
    //
    for(auto& node : res.nodes)
    {
      if(node.leftOffset == 0xFFFFFFFF)
        node.escapeIndex = 0xFFFFFFFF;
      else
        node.escapeIndex = node.leftOffset+1;
    }

    res.format = cbvh::FMT_BVH2Node32_Interval32_Dynamic;
  }
  else if(a_settings.childrenNum == 2 && a_settings.desiredFormat == cbvh::FMT_BVH2Node32_Interval32_Static)
    res.format = cbvh::FMT_BVH2Node32_Interval32_Static;

  if(res.nodes.size() == 1 && a_settings.alwaysNotEmpty) // this is important for FatBVH-like traversal which directly take children
  {
    BVHNode dummyNode;
    dummyNode.boxMin      = float3(0,0,0);
    dummyNode.boxMax      = float3(0,0,0);
    dummyNode.leftOffset  = 0xFFFFFFFD;
    dummyNode.escapeIndex = 0xFFFFFFFD; 
    
    res.nodes.push_back(res.nodes[0]);
    res.intervals.push_back(res.intervals[0]);
    for(int i=1;i<a_settings.childrenNum;i++) { 
      res.nodes.push_back(dummyNode);
      res.intervals.push_back(Interval(0,0));
    }
  }
  
  if(g_spatialSplits != 0)
    std::cout << "[cbvh]: g_spatialSplits = " << g_spatialSplits << std::endl;

  float time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - before).count()/1000.f;
  g_buildTime += double(time);

  return res;
}

std::vector<Box4f> MakeBoxesOnTriangles(const float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, bool a_discardLines)
{
  const size_t trianglesNum = a_indexNum/3;
  std::vector<Box4f> boxes;
  boxes.reserve(trianglesNum);
  boxes.resize(0);

  for(size_t i=0;i<trianglesNum;i++)
  {
    const float4 A = a_vertices[a_indices[i*3+0]];
    const float4 B = a_vertices[a_indices[i*3+1]];
    const float4 C = a_vertices[a_indices[i*3+2]];

    Box4f box;
    box.include(A);
    box.include(B);
    box.include(C);
    box.setStart(uint(i));
    box.setCount(uint(1));
    
    bool discardNow = false;
    if(a_discardLines)
    {
      const float4 flatNorm = cross(A-B,A-C);
      discardNow = dot3f(flatNorm,flatNorm) == 0.0f;
    }

    const float4 boxSize = box.boxMax - box.boxMin;    // #TODO: skip thin triangles also? make preset for this ??
    if(dot3f(boxSize,boxSize) > 1e-24f && !discardNow) // skip degenerative and may be line triangles
      boxes.push_back(box);                            //
  }                                                 
  
  return boxes;
}

cbvh::BVHTree cbvh_internal::BuildBVH4(const float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, BVHSettings a_settings)
{
  //if(a_settings.enableESC && a_settings.esc.enableClusters)
  //{
  //  auto splitted = cbvh_internal::ESC2(a_vertices, a_vertNum, a_indices, a_indexNum, a_settings.esc);
  //  auto bvhDataB = cbvh_internal::BuildBVH4(splitted.boxes.data(), splitted.boxes.size(), nullptr, 0, a_settings,
  //                                           cbvh_internal::OriginalInput(a_vertices, a_vertNum, a_indices, a_indexNum));
//
  //  std::vector<uint32_t> indicesReordered2(bvhDataB.indicesReordered.size());
  //  for(size_t i=0;i<bvhDataB.indicesReordered.size();i++)
  //    indicesReordered2[i] = splitted.boxes[bvhDataB.indicesReordered[i]].getStart();
  //  bvhDataB.indicesReordered = indicesReordered2;
  //  return bvhDataB;        
  //}
  //else if(a_settings.enableESC)
  //{
  //  auto splitted = cbvh_internal::ESC(a_vertices, a_vertNum, a_indices, a_indexNum, a_settings.esc);
  //  return cbvh_internal::BuildBVH4(splitted.boxes.data(), splitted.boxes.size(), splitted.indexBuffer.data(), splitted.indexBuffer.size(), a_settings,
  //                                  cbvh_internal::OriginalInput(a_vertices, a_vertNum, splitted.indexBuffer.data(), splitted.indexBuffer.size()));
  //}
  //else
  //{
    auto boxes = MakeBoxesOnTriangles(a_vertices, a_vertNum, a_indices, a_indexNum, a_settings.discardLines);
    return cbvh_internal::BuildBVH4(boxes.data(), boxes.size(), a_indices, a_indexNum, a_settings, 
                                    cbvh_internal::OriginalInput(a_vertices, a_vertNum, a_indices, a_indexNum));
  //}
}




static void PrintSpaces(int N, std::ostream& a_out)
{
  for(int i=0;i<N;i++)
    a_out << " ";
}

static void PrintNode(std::ostream& a_outFile, const cbvh::BVHTree a_treeData, uint32_t nodeId, int currDeep)
{
  const auto& currBox = a_treeData.nodes[nodeId];
  const auto& begiEnd = a_treeData.intervals[nodeId];

  PrintSpaces(currDeep*2, a_outFile);
  a_outFile << "node[" << nodeId << "]:" << std::endl;
  PrintSpaces(currDeep*2, a_outFile);
  a_outFile << "{" << std::endl;
  PrintSpaces(currDeep*2, a_outFile);  
  a_outFile << "  bmin = { " << currBox.boxMin[0] << " " << currBox.boxMin[1] << " " << currBox.boxMin[2] << " } | " << currBox.leftOffset  << std::endl;
  PrintSpaces(currDeep*2, a_outFile);  
  a_outFile << "  bmax = { " << currBox.boxMax[0] << " " << currBox.boxMax[1] << " " << currBox.boxMax[2] << " } | " << currBox.escapeIndex << std::endl;
  PrintSpaces(currDeep*2, a_outFile);  
  a_outFile << "  intv = { " << begiEnd.start << " " << begiEnd.count << " }" << std::endl;
  PrintSpaces(currDeep*2, a_outFile);
  a_outFile << "}" << std::endl;
  
  if(currBox.leftOffset < a_treeData.nodes.size())
  {
    for(int i=0;i<4;i++)
      PrintNode(a_outFile, a_treeData, currBox.leftOffset+i, currDeep+1);
  }
}

void cbvh_internal::DebugPrintBVH(const char* a_outFileName, const cbvh::BVHTree a_treeData)
{
  std::ofstream fout(a_outFileName);
  
  PrintNode(fout, a_treeData, 0, 0);

  fout << "bvhnods.size() = " << a_treeData.nodes.size() << std::endl;
  fout << "indices.size() = " << a_treeData.indicesReordered.size() << std::endl;
  fout << "indices = [";
  for(size_t i=0;i< a_treeData.indicesReordered.size();i++)
    fout << a_treeData.indicesReordered[i] << " ";
  fout << "]" << std::endl;
  fout << std::endl;
  fout.close();
}
