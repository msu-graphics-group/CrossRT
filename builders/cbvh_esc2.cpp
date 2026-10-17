#include <iostream>
#include <queue>
#include <stack>
#include <limits>
#include <cassert>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>

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

static inline float BoxMaxSize3f(float4 vmin, float4 vmax) 
{
  const float4 abc = vmax - vmin;
  return std::max(abc[0], std::max(abc[1], abc[2]));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

cbvh_internal::ESC_Clusters cbvh_internal::MakeTriangleClustersForESC(const LiteMath::float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, const ESC_Settings& a_settings)
{
  ESC_Clusters res;
  std::vector<Box4f>& a_clustersForSplit   = res.clustersForSplit;
  std::vector<Box4f>& a_theRestOfTriangles = res.theRestOfTriangles;
  std::vector<uint32_t>& a_indBuffReordered = res.indReord;
  res.sceneBox   = Box4f();
  res.avgTriSize = 0.0f;
  size_t totalTris = 0;

  cbvh_internal::BVHSettings settings;
  settings.primsInLeaf        = 16;
  settings.childrenNum        = 2;
  settings.discardLines       = true;
  settings.enableSS           = false;
  settings.esc.prebuildForESC = true; 
  
  auto preBuildBVH = cbvh_internal::BuildBVH4(a_vertices, a_vertNum, a_indices, a_indexNum, settings);

  std::vector<Box4f> leaves, others; 
  leaves.reserve(preBuildBVH.nodes.size());
  others.reserve(preBuildBVH.nodes.size());

  for(size_t i=0; i<preBuildBVH.nodes.size();i++)
  {
    const cbvh::BVHNode& node = preBuildBVH.nodes[i];
    if(node.leftOffset == uint32_t(-1))
    {
      const auto start = preBuildBVH.intervals[i].start; ///<! actual node id
      const auto count = preBuildBVH.intervals[i].count; ///<! store count here
      
      Box4f box;
      for(int k=0;k<3;k++)
      {
        box.boxMin[k] = node.boxMin[k];
        box.boxMax[k] = node.boxMax[k];
      }

      double metric    = 0.0;
      double metrixMax = 0.0;
      for(uint32_t triId=start; triId < start + count; triId++)
      {
        const uint32_t triIdOld = preBuildBVH.indicesReordered[triId]; 
        const uint32_t iA = a_indices[triIdOld*3+0];
        const uint32_t iB = a_indices[triIdOld*3+1];
        const uint32_t iC = a_indices[triIdOld*3+2];

        const float4 A = a_vertices[iA];
        const float4 B = a_vertices[iB];
        const float4 C = a_vertices[iC];

        Box4f triBox;
        triBox.include(A);
        triBox.include(B);
        triBox.include(C);
        const float triSA = 0.5f*length3f(cross(A-B,A-C));
        const float boxSA = triBox.surfaceArea();
        const double val = double(boxSA/std::max(triSA,1e-10f));
        metric    += val;
        metrixMax = std::max(val,metrixMax);

        res.sceneBox.include(triBox.boxMin);
        res.sceneBox.include(triBox.boxMax);
        res.avgTriSize += length3f(triBox.boxMax - triBox.boxMin);
        totalTris++;
      }
      
      box.boxMin.w = float(metric)*box.surfaceArea(); // sort them in such a way to subdivide large boxes first
      box.boxMax   = packUIntW(box.boxMax, i); 

      if(metrixMax > a_settings.thresholdMax)
        leaves.push_back(box);
      else
        others.push_back(box);
    }
  } 

  res.avgTriSize = res.avgTriSize / float(totalTris);

  std::sort(leaves.begin(), leaves.end(), [](const Box4f& s1, const Box4f& s2){ return s1.boxMin.w > s2.boxMin.w;} );

  a_clustersForSplit.resize(leaves.size());
  for(size_t i=0;i<a_clustersForSplit.size();i++)
  {
    a_clustersForSplit[i] = leaves[i];
    const uint32_t nodeId = extractUIntW(a_clustersForSplit[i].boxMax);
    const auto interval   = preBuildBVH.intervals[nodeId];
    a_clustersForSplit[i].setStart(interval.start);
    a_clustersForSplit[i].setCount(interval.count);
  }

  // put other clusters to 'a_theRestOfTriangles'
  //
  const uint32_t totalTrianglesNum = preBuildBVH.indicesReordered.size()/3;
  a_theRestOfTriangles.resize(0);
  a_theRestOfTriangles.reserve(totalTrianglesNum);

  for(size_t i=0;i<others.size();i++)
  {
    Box4f box             = others[i];
    const uint32_t nodeId = extractUIntW(box.boxMax);
    const auto interval   = preBuildBVH.intervals[nodeId];

    box.setStart(interval.start);
    box.setCount(interval.count);
    a_theRestOfTriangles.push_back(box);
  }

  a_indBuffReordered = preBuildBVH.indicesReordered;
  
  return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cbvh_internal::SplitTriangles(const cbvh_internal::TriCluster& nodeBox, float splitPos, int splitAxis,
                                   cbvh_internal::TriCluster& leftBox, cbvh_internal::TriCluster& rightBox)
{
  if (!nodeBox.AxisAligned(splitAxis, splitPos))
  {
    for(auto& tri : nodeBox.triList)
    {
      const float4 edges[3][2] = { {tri.A, tri.B},
                                   {tri.C, tri.A},
                                   {tri.B, tri.C} };
      
      // grow both boxes with vertices and edge-plane intersections
      //
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

    } // end for(auto& tri : nodeBox.triList)

  }
  //else
  //{
  //  leftBox.boxMin  = nodeBox.boxMin;
  //  leftBox.boxMax  = nodeBox.boxMax;
  //  rightBox.boxMin = nodeBox.boxMin;
  //  rightBox.boxMax = nodeBox.boxMax;
  //}
  
  leftBox.boxMax [splitAxis] = splitPos;
  rightBox.boxMin[splitAxis] = splitPos;

  // Intersect with original bounds.
  //
  leftBox.intersect(nodeBox);
  rightBox.intersect(nodeBox);

  // now process triangles
  //
  leftBox.triList.clear();
  rightBox.triList.clear();

  const float eps = 1e-6f;

  leftBox.escMetric  = 0.0f;
  rightBox.escMetric = 0.0f;

  const float boxSALeft  = Box4f(leftBox.boxMin, leftBox.boxMax).surfaceArea();
  const float boxSARight = Box4f(rightBox.boxMin, rightBox.boxMax).surfaceArea();

  for(auto& tri : nodeBox.triList)
  {
    const float4 leftCenter  = (leftBox.boxMax  + leftBox.boxMin) *0.5f;
    const float4 rightCenter = (rightBox.boxMax + rightBox.boxMin)*0.5f;
    const float4 leftHalfSZ  = (leftBox.boxMax  - leftBox.boxMin) *(0.5f + eps) + float4(eps,eps,eps,eps);
    const float4 rightHalfSZ = (rightBox.boxMax - rightBox.boxMin)*(0.5f + eps) + float4(eps,eps,eps,eps);

	  float triverts[3][3];
	  for(int i=0;i<3;i++)
	  {
	 	  triverts[0][i] = tri.A[i];
	 	  triverts[1][i] = tri.B[i];
	 	  triverts[2][i] = tri.C[i];
	  }
    
    const int overlapLeft  = EXTERNAL_TOOLS::triBoxOverlap( (float*)&leftCenter,  (float*)&leftHalfSZ,  triverts);
    const int overlapRight = EXTERNAL_TOOLS::triBoxOverlap( (float*)&rightCenter, (float*)&rightHalfSZ, triverts);
    
    const float triSA      = std::max(0.5f*length3f(cross(tri.A-tri.B, tri.A-tri.C)), 1e-10f);
 
    if(overlapLeft != 0)
    {
      leftBox.triList.push_back(tri);
      const float metric = boxSALeft/triSA;
      leftBox.escMetric  = std::max(leftBox.escMetric, metric);
    }

    if(overlapRight != 0)
    {
      rightBox.triList.push_back(tri);
      const float metric = boxSARight/triSA;
      rightBox.escMetric = std::max(leftBox.escMetric, metric);
    }
  }
  
  return;
}

using cbvh_internal::TriCluster;

static inline uint32_t EstimateMaxSubdivs(const TriCluster& a_triCluster, const cbvh_internal::ESC_Settings& a_settings, size_t totalTrisNum)
{
  //if(totalTrisNum > 16384) // dynamic estimation don't help for large models unfortunately
  //  return 4;

  const float boxSA = Box4f(a_triCluster.boxMin, a_triCluster.boxMax).surfaceArea();
  float triSAMin    = std::numeric_limits<float>::max();

  for(auto& tri : a_triCluster.triList)
  {
    const float triSA = std::max(0.5f*length3f(cross(tri.A-tri.B, tri.A-tri.C)), 1e-10f);
    triSAMin = std::min(triSAMin, triSA);
  }

  const float relation = boxSA/triSAMin;
  
  if (relation >= 500.0f)
    return 64;
  else if (relation >= 200.0f)
    return 16;
  else 
    return 4;
}

static inline std::vector<TriCluster> SubdivideTriangles(const TriCluster& a_triCluster, const cbvh_internal::ESC_Settings& a_settings, size_t totalTrisNum)
{ 
  //const auto maxSubdivs = a_settings.maxSubdivs;
  const auto maxSubdivs = EstimateMaxSubdivs(a_triCluster, a_settings, totalTrisNum);

  std::vector<TriCluster> subdivs; //
  subdivs.reserve(maxSubdivs*2+10); //
  subdivs.resize(0);                //
  
  std::queue<TriCluster> queue;
  queue.push(a_triCluster);
  
  while(!queue.empty())
  {
    const cbvh_internal::TriCluster curr = queue.front(); queue.pop();
    
    // select split axis
    //
    const float4 boxSize   = curr.boxMax - curr.boxMin;
    const float4 boxCenter = 0.5f*(curr.boxMax + curr.boxMin);

    float splitPos  = 0;
    int   splitAxis = 0;
    float maxSize   = -1;

    for(int i=0;i<3;i++)
    {
      if(maxSize < boxSize[i])
      {
        splitPos  = boxCenter[i];
        splitAxis = i;
        maxSize   = boxSize[i];
      }
    }
    
    TriCluster a,b;
    SplitTriangles(curr, splitPos, splitAxis, 
                   a, b);

    if(a.escMetric > a_settings.thresholdMin && a.triList.size() > 0 && (queue.size() + subdivs.size() < maxSubdivs))
      queue.push(a);
    else if (a.triList.size() > 0)
      subdivs.push_back(a);

    if(b.escMetric > a_settings.thresholdMin && b.triList.size() > 0 && (queue.size() + subdivs.size() < maxSubdivs))
      queue.push(b);
    else if (b.triList.size() > 0)
      subdivs.push_back(b);
  }
  
  return subdivs;
}

TriCluster cbvh_internal::ReadTrianglesFromPrimBox(const Box4f& prim, const LiteMath::float4* a_vertices, size_t a_vertNum, 
                                                   const uint* a_indices, size_t a_indexNum,
                                                   const uint* a_indicesOld, size_t a_indexNumOld)
{
  TriCluster currCluster;
  currCluster.boxMin = prim.boxMin;
  currCluster.boxMax = prim.boxMax;
  currCluster.triList.reserve(10);

  const uint start = extractUIntW(prim.boxMin);
  const uint size  = extractUIntW(prim.boxMax);

  for(uint triIndex=start; triIndex < start + size; triIndex++)
  {
    assert(triIndex < a_indexNumOld);
    const uint triIndexOld = a_indices[triIndex];
    assert(3*triIndexOld+2 < a_indexNumOld);
     
    const uint iA = a_indicesOld[3*triIndexOld+0];
    const uint iB = a_indicesOld[3*triIndexOld+1];
    const uint iC = a_indicesOld[3*triIndexOld+2];
    
    assert(iA < a_vertNum);
    assert(iB < a_vertNum);
    assert(iC < a_vertNum);

    const float4 A = a_vertices[iA];
    const float4 B = a_vertices[iB];
    const float4 C = a_vertices[iC];

    cbvh_internal::Triangle4f tri;
    tri.A = packUIntW(A, iA); 
    tri.B = packUIntW(B, iB); 
    tri.C = packUIntW(C, iC); 
    currCluster.triList.push_back(tri);
  }

  return currCluster;
}

void SplitTriangleClusters(const LiteMath::float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, const uint* a_indicesOld, size_t a_indexNumOld, 
                           const std::vector<Box4f>& a_inputClusters, float4 a_boxMin, float4 a_boxMax, const cbvh_internal::ESC_Settings& a_settings,
                           cbvh_internal::ESC_Result& a_res, bool a_debug)
{
  
  const size_t totalTrisNum   = a_indexNum/3;
  const float  EXPAND_FACTOR  = std::max(a_settings.expandFactor, 1.15f);     // we should not generate more that EXPAND_FACTOR*a_inputClusters.size() primitives in total
  const size_t maxExpandSize  = size_t(std::max(float(totalTrisNum)*(EXPAND_FACTOR - 1.0f), 16384.0f));
  
  a_res.boxes.reserve(maxExpandSize + totalTrisNum + 10);
  a_res.boxes.resize(0);
  a_res.indexBuffer.resize(0);

  bool memoryExausted = false;
 
  for(size_t i = 0; i < a_inputClusters.size(); i++)
  {
    auto currCluster = cbvh_internal::ReadTrianglesFromPrimBox(a_inputClusters[i], a_vertices, a_vertNum, a_indices, a_indexNum, a_indicesOld, a_indexNumOld);

    std::vector<TriCluster> splitted;
    if(memoryExausted)
    {
      if(a_debug)
        break;
      splitted.push_back(currCluster);
    }
    else
      splitted = SubdivideTriangles(currCluster, a_settings, totalTrisNum);

    for(const auto& cluster : splitted)
    {
      for(const auto& tri : cluster.triList)
      {
        const uint32_t iA = extractUIntW(tri.A);
        const uint32_t iB = extractUIntW(tri.B);
        const uint32_t iC = extractUIntW(tri.C);
        
        const uint triBegin = uint(a_res.indexBuffer.size()/3);
        a_res.indexBuffer.push_back(iA);
        a_res.indexBuffer.push_back(iB);
        a_res.indexBuffer.push_back(iC);
  
        Box4f triBox = Box4f();
        triBox.include(a_vertices[iA]);
        triBox.include(a_vertices[iB]);
        triBox.include(a_vertices[iC]);
  
        Box4f pb(cluster.boxMin, cluster.boxMax);
        pb.intersect(triBox);
        pb.setStart(triBegin);
        pb.setCount(1);
        a_res.boxes.push_back(pb);
      }
    }
    
    if(!memoryExausted && a_settings.expandFactor > 0 && a_res.boxes.size() >= maxExpandSize) 
    {
      memoryExausted = true;
      //std::cout << "[esc]: memory exausted at " << int(100.0f*float(i)/float(a_inputClusters.size())) << "%" << std::endl; // #TODO: put std::cout out of library
    } 
  }

}  

void ChangeExpandFactorViaHeuritic(size_t a_indexNum, cbvh_internal::ESC_Settings& a_settings)
{   
  const auto totalTris = a_indexNum/3;

  if(totalTris <= 100000)
    a_settings.expandFactor = 1.5f;  
  else if(totalTris <= 500000)
    a_settings.expandFactor = 1.25f;  
  else
    a_settings.expandFactor = 1.15f;
}

cbvh_internal::ESC_Result cbvh_internal::ESC2(const LiteMath::float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, ESC_Settings a_settings, bool a_debug)
{
  assert(a_indexNum%3 == 0);
  assert(a_indexNum   != 0);
  
  if(!a_debug)
    ChangeExpandFactorViaHeuritic(a_indexNum, a_settings);

  auto clustRes = cbvh_internal::MakeTriangleClustersForESC(a_vertices, a_vertNum, a_indices, a_indexNum, a_settings);
  
  ESC_Result res;
  if(clustRes.clustersForSplit.size() != 0)
  {
    SplitTriangleClusters(a_vertices, a_vertNum, clustRes.indReord.data(), clustRes.indReord.size(), a_indices, a_indexNum, 
                          clustRes.clustersForSplit, clustRes.sceneBox.boxMin, clustRes.sceneBox.boxMax, 
                          a_settings, res, a_debug);
  }
  
  if(a_debug)
    return res;

  // append back unsplitted boxes back
  //
  for(size_t bi = 0; bi < clustRes.theRestOfTriangles.size(); bi++) 
  {
    Box4f box = clustRes.theRestOfTriangles[bi];
    const uint32_t start = box.getStart();
    const uint32_t count = box.getCount();
    
    for(uint32_t triId = start; triId < start + count; triId++)
    { 
      const uint32_t oldTriId = clustRes.indReord[triId];
      const auto iA = a_indices[oldTriId*3 + 0];
      const auto iB = a_indices[oldTriId*3 + 1];
      const auto iC = a_indices[oldTriId*3 + 2];

      //const uint triBegin = uint(res.indexBuffer.size()/3);
      //res.indexBuffer.push_back(iA);
      //res.indexBuffer.push_back(iB);
      //res.indexBuffer.push_back(iC);

      Box4f triBox = Box4f();
      triBox.include(a_vertices[iA]);
      triBox.include(a_vertices[iB]);
      triBox.include(a_vertices[iC]);
      triBox.setStart(triId);
      triBox.setCount(1);
      res.boxes.push_back(triBox);
    }
  }
 
  return res;
} 