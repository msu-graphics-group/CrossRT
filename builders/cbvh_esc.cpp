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

struct TriBox
{
  Box4f      box;
  Triangle4f tri;
};

std::vector<TriBox> SelectBadTrianglesForESC(const LiteMath::float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, const cbvh_internal::ESC_Settings& a_settings, 
                                             std::vector<Box4f>& a_theRestOfTriangles)
{
  const uint32_t a_trisNum = uint32_t(a_indexNum/3);
  
  std::vector<TriBox> selectedTris;
  selectedTris.reserve(a_trisNum);
  a_theRestOfTriangles.reserve(a_trisNum);

  for(uint32_t triId=0; triId < a_trisNum; triId++)
  {
    const uint32_t iA = a_indices[triId*3+0];
    const uint32_t iB = a_indices[triId*3+1];
    const uint32_t iC = a_indices[triId*3+2];
    const float4 A = a_vertices[iA];
    const float4 B = a_vertices[iB];
    const float4 C = a_vertices[iC];

    Box4f triBox;
    triBox.include(A);
    triBox.include(B);
    triBox.include(C);
    triBox.setStart(triId);

    const float triSA = 0.5f*length3f(cross(A-B,A-C));
    const float boxSA = triBox.surfaceArea();
    const float val   = boxSA/std::max(triSA,1e-10f);
    triBox.boxMax.w   = boxSA*val;

    if(val > a_settings.thresholdMax)
    {
      TriBox tb;
      tb.box   = triBox;
      tb.tri.A = packUIntW(A, iA);
      tb.tri.B = packUIntW(B, iB);
      tb.tri.C = packUIntW(C, iC); 
      selectedTris.push_back(tb);
    }
    else
      a_theRestOfTriangles.push_back(triBox);
  }

  return selectedTris;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void SplitTriangle(const Box4f& nodeBox, float splitPos, int splitAxis, const Triangle4f& tri,
                          Box4f& leftBox, Box4f& rightBox)
{
  leftBox  = nodeBox;
  rightBox = nodeBox;
  
  const auto triId = nodeBox.getStart();
  
  if (!nodeBox.isAxisAligned(splitAxis, splitPos))
  {
    const float4 edges[3][2] = {  {tri.A, tri.B},
                                  {tri.C, tri.A},
                                  {tri.B, tri.C} };
     
    // grow both boxes with vertices and edge-plane intersections
    //
    for (int i=0;i<3;i++)
    {
      float v0p = edges[i][0][splitAxis];
      float v1p = edges[i][1][splitAxis];
      
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
  }
  
  leftBox.boxMax [splitAxis] = splitPos;
  rightBox.boxMin[splitAxis] = splitPos;

  // Intersect with original bounds.
  //
  leftBox.intersect(nodeBox);
  rightBox.intersect(nodeBox);
  
  leftBox.setStart(triId);
  leftBox.setCount(1);
  rightBox.setStart(triId);
  rightBox.setCount(1);
}

static inline std::vector<Box4f> SubdivideTriangle(const Box4f& a_triBox, const Triangle4f& a_tri, const int maxSubdivs)
{ 
  std::vector<Box4f> subdivs;      //
  subdivs.reserve(maxSubdivs + 10); //
  subdivs.resize(0);                //
  
  std::queue<Box4f> queue;
  queue.push(a_triBox);
  
  while(!queue.empty())
  {
    const auto curr = queue.front(); queue.pop();
    
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

    Box4f a,b;
    SplitTriangle(curr, splitPos, splitAxis, a_tri,
                  a, b);
    
    const float eps          = 1e-6f;
    const float4 leftCenter  = (a.boxMax  + a.boxMin) *0.5f;
    const float4 rightCenter = (b.boxMax + b.boxMin)*0.5f;
    const float4 leftHalfSZ  = (a.boxMax  - a.boxMin) *(0.5f + eps) + float4(eps,eps,eps,0.0f);
    const float4 rightHalfSZ = (b.boxMax - b.boxMin)*(0.5f + eps) + float4(eps,eps,eps,0.0f);

	  float triverts[3][3];
	  for(int i=0;i<3;i++)
	  {
	 	  triverts[0][i] = a_tri.A[i];
	 	  triverts[1][i] = a_tri.B[i];
	 	  triverts[2][i] = a_tri.C[i];
	  }
    
    const int overlapA = EXTERNAL_TOOLS::triBoxOverlap( (float*)&leftCenter,  (float*)&leftHalfSZ,  triverts);
    const int overlapB = EXTERNAL_TOOLS::triBoxOverlap( (float*)&rightCenter, (float*)&rightHalfSZ, triverts);

    if((queue.size() + subdivs.size() < size_t(maxSubdivs)) && overlapA)
      queue.push(a);
    else if(overlapA)
      subdivs.push_back(a);

    if((queue.size() + subdivs.size() < size_t(maxSubdivs)) && overlapB)
      queue.push(b);
    else if(overlapB)
      subdivs.push_back(b);
  }
  
  return subdivs;
}

static inline int EstimateSubdivs(float relation)
{
  //int maxSubdivs; //= (int)(clamp(1.0f/relation,2.0f,128.0f));
  //
  //if(relation > 512.f)
  //  maxSubdivs = 256;
  //else if (relation > 256.f)
  //  maxSubdivs = 128;
  //else if (relation > 128.f)
  //  maxSubdivs = 64;
  //else if(relation > 64.f)
  //  maxSubdivs = 32;
  //else if (relation > 32.f)
  //  maxSubdivs = 16;
  //else
  //  maxSubdivs = 4;
  //
  //return std::min(std::max(maxSubdivs,4), 16);

  if (relation >= 500.0f)
    return 64;
  else if (relation >= 200.0f)
    return 16;
  else 
    return 4;
}

cbvh_internal::ESC_Result cbvh_internal::ESC(const LiteMath::float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, ESC_Settings a_settings, bool a_debug)
{
  assert(a_indexNum%3 == 0);
  assert(a_indexNum   != 0);

  const size_t numTriangles = a_indexNum/3;
  const size_t maxSplit     = std::max(size_t(float(numTriangles)*a_settings.expandFactor), size_t(16384));

  ESC_Result res;
  res.indexBuffer.reserve(maxSplit*3 + a_indexNum + 3000);
  res.boxes.reserve(maxSplit + numTriangles + 1000);

  std::vector<Box4f> theRestOfTriangles;
  auto selectedTriBoxes = SelectBadTrianglesForESC(a_vertices, a_vertNum, a_indices, a_indexNum, a_settings, 
                                                   theRestOfTriangles);

  std::sort(selectedTriBoxes.begin(), selectedTriBoxes.end(), 
            [](const auto& b1, const auto& b2) { return b1.box.boxMax.w > b2.box.boxMax.w; });
  
  bool memExeeded = false;
  for(size_t triId=0; triId<selectedTriBoxes.size(); triId++)
  {
    const auto& currTriBox = selectedTriBoxes[triId];
    const auto& currBox = currTriBox.box;
    const auto& currTri = currTriBox.tri;
    const float val     = currBox.boxMax.w/currTriBox.box.surfaceArea();

    const uint32_t iA = extractUIntW(currTri.A);
    const uint32_t iB = extractUIntW(currTri.B);
    const uint32_t iC = extractUIntW(currTri.C);

    if(memExeeded)
    {
      const uint32_t triBegin = uint32_t(res.indexBuffer.size()/3);
      res.indexBuffer.push_back(iA);
      res.indexBuffer.push_back(iB);
      res.indexBuffer.push_back(iC);

      auto splitBox = currBox;
      splitBox.setStart(triBegin);
      splitBox.setCount(1);
      res.boxes.push_back(splitBox);
      continue;
    }

    auto splitted = SubdivideTriangle(currBox, currTri, EstimateSubdivs(val));
    for(auto splitBox : splitted)
    {
      const uint32_t triBegin = uint32_t(res.indexBuffer.size()/3);
    
      res.indexBuffer.push_back(iA);
      res.indexBuffer.push_back(iB);
      res.indexBuffer.push_back(iC);
  
      splitBox.setStart(triBegin);
      splitBox.setCount(1);
      res.boxes.push_back(splitBox);
    }

    memExeeded = (res.boxes.size() >= maxSplit);
  }

  if(a_debug)
    return res;

  for(auto otherBox : theRestOfTriangles)
  {
    const uint32_t triIdOld = otherBox.getStart();
    const uint32_t triBegin = uint32_t(res.indexBuffer.size()/3);
    
    const uint32_t iA = a_indices[triIdOld*3+0];
    const uint32_t iB = a_indices[triIdOld*3+1];
    const uint32_t iC = a_indices[triIdOld*3+2];
    res.indexBuffer.push_back(iA);
    res.indexBuffer.push_back(iB);
    res.indexBuffer.push_back(iC);

    otherBox.setStart(triBegin);
    otherBox.setCount(1);
    res.boxes.push_back(otherBox);
  }

  return res;
} 