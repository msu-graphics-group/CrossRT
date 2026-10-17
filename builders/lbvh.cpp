#include "lbvh.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <fstream>
#include <sstream>

using cbvh2::make_BVHNode;
using cbvh2::DummyNode;

using LiteMath::sign;
using LiteMath::dot;
using LiteMath::min;
using LiteMath::max;
using LiteMath::abs;
using LiteMath::uint;
using LiteMath::uint2;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

constexpr static unsigned int ZINDEXSIZE = 1024;

static unsigned int SpreadBits(int x, int offset)
{
  //assert(x >=0      && unsigned(x) < ZINDEXSIZE);
  //assert(offset >=0 && offset < 3);
  x = (x | (x << 10)) & 0x000F801F;
  x = (x | (x <<  4)) & 0x00E181C3;
  x = (x | (x <<  2)) & 0x03248649;
  x = (x | (x <<  2)) & 0x09249249;
  return uint(x) << offset;
}
  
static unsigned int ZIndex3D(int x, int y, int z)
{
  return SpreadBits(x, 0) | SpreadBits(y, 1) | SpreadBits(z, 2);
}

static inline float4 ZIndexMultCoeff(const float4 bMinMesh, const float4 bMaxMesh)
{
  float4 bSizeMeshInv = 1.0f/(bMaxMesh - bMinMesh);
  for(int i=0;i<3;i++)
  {
    if(std::abs(bMaxMesh[i] - bMinMesh[i]) < 1e-20f) // std::numeric_limits<float>::min()
      bSizeMeshInv[i] = 0.0f;
  }
  
  return float4(float(ZINDEXSIZE))*bSizeMeshInv;
}

constexpr float LBVH_INF      = +1e10;
constexpr float LBVH_EPS_BIAS = 1e-4f;

static constexpr uint32_t START_MASK = 0x00FFFFFF;
static constexpr uint32_t END_MASK   = 0xFF000000;
static constexpr uint32_t SIZE_MASK  = 0x7F000000;
static constexpr uint32_t LEAF_BIT   = 0x80000000;

static inline uint32_t PackOffsetAndSize(uint32_t start, uint32_t size)
{
  return LEAF_BIT | ((size << 24) & SIZE_MASK) | (start & START_MASK);
}

static inline uint32_t EXTRACT_START(uint32_t a_leftOffset)  { return  a_leftOffset & START_MASK; }
static inline uint32_t EXTRACT_COUNT(uint32_t a_leftOffset)  { return (a_leftOffset & SIZE_MASK) >> 24; }

size_t LBVHBuilder::roundUpTo(size_t a_size, size_t a_num)
{
  while(a_size % a_num != 0)
    a_size++;
  return a_size;
}

size_t LBVHBuilder::powerOfTwo(size_t a_size)
{
  size_t currSize2 = 2;
  while(currSize2 < a_size) // for reserving runCodesAndIndices always with power of 2 size
    currSize2 *= 2;
  return currSize2;
}

size_t LBVHBuilder::Reserve(size_t a_inPrimsNum, cbvh2::FMT a_treeFormat)
{
  m_treeFormat = uint32_t(a_treeFormat);
  const size_t childs = (a_treeFormat == cbvh2::BVH4_LEFT_OFFSET) ? 4 : 2; 
  size_t currOffs  = 0;
  size_t currSize  = 2;
  size_t peachSize = (a_inPrimsNum <= PEACHES_MAX*4) ? PEACHES_MIN : PEACHES_MAX;
  m_peachIntervals.clear();

  if(a_inPrimsNum == 1)
  {
    m_peachIntervals.push_back(uint2(0,childs));
    currOffs = childs;
  }
  else
  { 
    m_peachIntervals.reserve(32);
    do
    {
      m_peachIntervals.push_back(uint2(currOffs,currSize));
      currOffs += currSize;
      currSize *= childs;
    } while(currSize < peachSize);
    
    m_tempCountPrefix.reserve(currSize + 1024); // for prevent prefix summ (scan) bug on ARM
    m_tempCountPrefix.resize(currSize);
    m_peachIntervals.push_back(uint2(currOffs,currSize));
    currOffs += currSize;
    currSize *= childs;
    currOffs += currSize;                            // for first-level nodes we always have from 'kernel1D_RecursiveConvert'
    currOffs += roundUpTo(2*a_inPrimsNum-1, childs); // for all other nodes that we actually convert recursively
  }

  size_t currSize2 = powerOfTwo(a_inPrimsNum);
  runCodesAndIndices.resize(currSize2);
  
  codesEq.resize(a_inPrimsNum);
  prefixCodeEq.reserve(a_inPrimsNum + 1024); // for prevent prefix summ (scan) bug on ARM
  prefixCodeEq.resize(a_inPrimsNum);
  
  newIndices.reserve(a_inPrimsNum);
  compressedCodes.reserve(a_inPrimsNum);
  m_tempCount.resize(2*a_inPrimsNum-1);
  
  if(a_treeFormat == cbvh2::BVH2_LEFT_RIGHT)
    return 2*a_inPrimsNum-1;
  else
  {
    m_tempNodes.resize(2*a_inPrimsNum-1);
    m_tempAddresses.resize(currOffs);
    m_leftAndEscape.resize(currOffs);
    BuildPeachTree(m_peachIntervals, m_leftAndEscape);
    m_reservedNodes = currOffs;
    return currOffs; // for cbvh2::BVH2_LEFT_OFFSET or cbvh2::BVH4_LEFT_OFFSET (peach trees)
  }
}

static void PushDownEscapeIndex(std::vector<uint2>& a_leftAndEscape, uint32_t a_leftOffset, uint32_t a_parentEscapeIndex)
{
  if(a_leftOffset >= a_leftAndEscape.size())
    return;
  
  const uint2 leftNode  = a_leftAndEscape[a_leftOffset + 0];
  const uint2 rightNode = a_leftAndEscape[a_leftOffset + 1];
  
  a_leftAndEscape[a_leftOffset + 0].y = a_leftOffset + 1;
  a_leftAndEscape[a_leftOffset + 1].y = a_parentEscapeIndex;

  PushDownEscapeIndex(a_leftAndEscape, leftNode.x, a_leftOffset + 1);      //
  PushDownEscapeIndex(a_leftAndEscape, rightNode.x, a_parentEscapeIndex);  //
}

void LBVHBuilder::BuildPeachTree(const std::vector<uint2>& a_peachIntervals, std::vector<uint2>& a_leftAndEscape)
{
  // (1) init all as leaves
  //
  for(size_t i = 0; i < a_leftAndEscape.size();i++)
    a_leftAndEscape[i] = uint2(0xFFFFFFFF,0);

  // (2) init left offsets for non-leaf nodes
  //
  for(size_t level=0; level < a_peachIntervals.size()-1; level++)
  {
    const uint2 curr = a_peachIntervals[level+0];
    const uint2 next = a_peachIntervals[level+1];
    for(uint32_t i = curr.x; i < curr.x + curr.y;i++)
      a_leftAndEscape[i] = uint2(next.x + 2*(i - curr.x), 0);
  }
  
  //static int counter = 0;
  //counter++;
  //std::stringstream strOut;
  //strOut << "z" << counter << "_pass";
  //std::string fileName = strOut.str();
  //std::ofstream fout(fileName.c_str());
  //for(size_t i=0; i<a_leftAndEscape.size(); i++)
  //  fout << a_leftAndEscape[i].x << " " << a_leftAndEscape[i].y << std::endl;
  //
  //for(size_t level=0; level < a_peachIntervals.size(); level++)
  //  std::cout << "(" << a_peachIntervals[level].x << ", " << a_peachIntervals[level].y << ")" << std::endl;

  // (3) make correct escapeIndices
  //
  PushDownEscapeIndex(a_leftAndEscape, 0, uint32_t(-2));
}

void LBVHBuilder::kernel1D_RedunctionFromBoxes(const float4* a_vertices, uint32_t a_size)
{
  bboxMin = float4(+LBVH_INF, +LBVH_INF, +LBVH_INF, +LBVH_INF);
  bboxMax = float4(-LBVH_INF, -LBVH_INF, -LBVH_INF, -LBVH_INF);  

  for (uint32_t boxId = 0; boxId < a_size; boxId++)
  {
    const float4 triBoxMin = a_vertices[2*boxId+0];
    const float4 triBoxMax = a_vertices[2*boxId+1];
    
    bboxMin = min(bboxMin, triBoxMin);
    bboxMax = max(bboxMax, triBoxMax);
  }
}

void LBVHBuilder::kernel1D_RedunctionFromTriangles(const float4* a_vertices, uint32_t a_vertNum, const uint32_t* a_indices, uint32_t a_trisNum)
{
  bboxMin = float4(+LBVH_INF, +LBVH_INF, +LBVH_INF, +LBVH_INF);
  bboxMax = float4(-LBVH_INF, -LBVH_INF, -LBVH_INF, -LBVH_INF);  

  for (uint32_t triId = 0; triId < a_trisNum; triId++)
  {
    const uint32_t A = a_indices[triId*3+0];
    const uint32_t B = a_indices[triId*3+1];
    const uint32_t C = a_indices[triId*3+2];
    
    const float4 v1  = a_vertices[A];
    const float4 v2  = a_vertices[B];
    const float4 v3  = a_vertices[C];
    
    const float4 triBoxMin = min(v1, min(v2, v3));
    const float4 triBoxMax = max(v1, max(v2, v3));
    
    bboxMin = min(bboxMin, triBoxMin);
    bboxMax = max(bboxMax, triBoxMax);
  }
}

void LBVHBuilder::kernel1D_RedunctionFromVertices(const float4* a_vertices, uint32_t a_vertNum)
{
  bboxMin = float4(+LBVH_INF, +LBVH_INF, +LBVH_INF, +LBVH_INF);
  bboxMax = float4(-LBVH_INF, -LBVH_INF, -LBVH_INF, -LBVH_INF);  

  for (uint32_t vertId = 0; vertId < a_vertNum; vertId++)
  {
    const float4 v1 = a_vertices[vertId];  
    bboxMin = min(bboxMin, v1);
    bboxMax = max(bboxMax, v1);
  }
}

void LBVHBuilder::kernel1D_EvalMortonCodesFromBoxes(const float4* a_inBoxes, uint32_t a_boxCount, uint32_t a_powerOf2Size, uint2* outRunCodes)
{
  //zIndexScale = ZIndexMultCoeff(bboxMin - abs(bboxMin)*LBVH_EPS_BIAS, bboxMax + abs(bboxMin)*LBVH_EPS_BIAS);
  zIndexScale = ZIndexMultCoeff(bboxMin - float4(LBVH_EPS_BIAS), bboxMax + float4(LBVH_EPS_BIAS));       
  for (uint32_t i = 0; i < a_powerOf2Size; ++i)
  {
    if(i < a_boxCount)
    {
      const float4 triBoxMin    = a_inBoxes[2*i+0];
      const float4 triBoxMax    = a_inBoxes[2*i+1];
      const float4 scaledVertex = zIndexScale*((triBoxMin + triBoxMax) * 0.5f - bboxMin);
      outRunCodes[i] = uint2(ZIndex3D(scaledVertex.x, scaledVertex.y, scaledVertex.z), i);
    }
    else
      outRunCodes[i] = uint2(0xFFFFFFFF, 0xFFFFFFFF); // due to bitonic sort
  }
}

void LBVHBuilder::kernel1D_EvalMortonCodesFromTriangles(const float4* a_vertices, uint32_t a_vertNum, const uint32_t* a_indices, uint32_t a_trisNum, uint32_t a_powerOf2Size,
                                                        uint2* outRunCodes)
{
  //zIndexScale = ZIndexMultCoeff(bboxMin - abs(bboxMin)*LBVH_EPS_BIAS, bboxMax + abs(bboxMin)*LBVH_EPS_BIAS);
  zIndexScale = ZIndexMultCoeff(bboxMin - float4(LBVH_EPS_BIAS), bboxMax + float4(LBVH_EPS_BIAS));       
  for (uint32_t i = 0; i < a_powerOf2Size; ++i)
  {
    if(i < a_trisNum)
    {
      const uint32_t A = a_indices[i*3+0];
      const uint32_t B = a_indices[i*3+1];
      const uint32_t C = a_indices[i*3+2];
      
      const float4 v1  = a_vertices[A];
      const float4 v2  = a_vertices[B];
      const float4 v3  = a_vertices[C];
      
      const float4 triBoxMin = min(v1, min(v2, v3));
      const float4 triBoxMax = max(v1, max(v2, v3));

      const float4 center       = (a_trisNum < 256) ? 0.33333333f*(v1 + v2 + v3) : (triBoxMin + triBoxMax) * 0.5f; 
      const float4 scaledVertex = zIndexScale * (center - bboxMin);

      outRunCodes[i] = uint2(ZIndex3D(scaledVertex.x, scaledVertex.y, scaledVertex.z), i);
    }
    else
      outRunCodes[i] = uint2(0xFFFFFFFF, 0xFFFFFFFF); // we must fill tailf of the array because of bitonic sort
  }
}

void LBVHBuilder::kernel1D_AppendInit(const uint2* in_runCodes, uint32_t a_boxCount, 
                                      uint32_t* a_codesEq, uint32_t* a_codesPref)
{
  for (uint32_t i = 0; i < a_boxCount; ++i) 
  {
    uint32_t codeEQ   = 1;
    if(i > 0)
      codeEQ = (in_runCodes[i - 1].x != in_runCodes[i].x) ? 1 : 0;
    a_codesEq  [i] = codeEQ;
    a_codesPref[i] = codeEQ;
  }
}

void LBVHBuilder::kernel1D_AppendComplete(uint32_t* a_codesEq, uint32_t* a_codesPref, uint32_t a_boxCount)
{
  const uint32_t leafNumber = a_codesPref[a_boxCount-1] + a_codesEq[a_boxCount - 1];
  newIndices.resize(leafNumber);
  compressedCodes.resize(leafNumber);
  m_leafNumber         = int(leafNumber);
  m_leafNumberMinusOne = int(leafNumber-1);

  for (uint32_t i = 0; i < a_boxCount; i++)
  {
    if (a_codesEq[i] != 0)                                 // Parallel Append, transform it via prefix summ
    {
      const uint32_t writeId   = a_codesPref[i];  
      newIndices     [writeId] = i;                        // newIndices.push_back(i);
      compressedCodes[writeId] = runCodesAndIndices[i].x;  // compressedCodes.push_back(runCodesAndIndices[i].x);
    }
  }
}

void LBVHBuilder::kernel1D_MakeLeavesFromBoxes(const float4* a_inBoxes, uint32_t a_boxCount, 
                                               BVHNode* a_nodes, uint32_t* a_outIndices)
{
  for (uint32_t i = 0; i < newIndices.size(); ++i) //
  {
    const uint32_t start = newIndices[i];
    const uint32_t count = (i < newIndices.size()-1) ? newIndices[i + 1] - start : a_boxCount - newIndices[newIndices.size()-1];
    
    float4 bmin(LBVH_INF, LBVH_INF, LBVH_INF, LBVH_INF);
    float4 bmax(-LBVH_INF, -LBVH_INF, -LBVH_INF, -LBVH_INF);
    
    for (uint32_t j = start; j < (start + count); j++) 
    {
      const uint32_t boxId = runCodesAndIndices[j].y;
      a_outIndices[j] = boxId;
      
      const float4 triBoxMin = a_inBoxes[2*boxId+0];
      const float4 triBoxMax = a_inBoxes[2*boxId+1];

      bmin = min(bmin, triBoxMin);
      bmax = max(bmax, triBoxMax);
    }
    
    uint32_t leafOffset = i + newIndices.size() - 1;
    //a_nodes[leafOffset] = make_BVHNode(to_float3(bmin), to_float3(bmax), PackOffsetAndSize(start, count), 0xFFFFFFFF);                 // correct way for TLAS, multiple BLAS per leaf
    a_nodes[leafOffset] = make_BVHNode(to_float3(bmin), to_float3(bmax), PackOffsetAndSize(runCodesAndIndices[start].y, 1), 0xFFFFFFFF); // restricted way for TLAS, single BLAS per leaf
  }
}                                  

void LBVHBuilder::kernel1D_MakeLeavesFromTriangles(const float4* a_vertices, uint32_t a_vertNum, const uint32_t* a_indices, uint32_t a_indexNum,
                                                   BVHNode* a_nodes, uint32_t* a_outIndices)
{
  for (uint32_t i = 0; i < newIndices.size(); ++i) //
  {
    const uint32_t start = newIndices[i];
    const uint32_t count = (i < newIndices.size()-1) ? newIndices[i + 1] - start : a_indexNum/3 - newIndices[newIndices.size()-1];
    
    float4 bmin(LBVH_INF, LBVH_INF, LBVH_INF, LBVH_INF);
    float4 bmax(-LBVH_INF, -LBVH_INF, -LBVH_INF, -LBVH_INF);
    
    for (uint32_t j = start; j < (start + count); j++) 
    {
      const uint32_t boxId = runCodesAndIndices[j].y;
      a_outIndices[j] = boxId;
      
      const uint32_t A = a_indices[boxId*3+0];
      const uint32_t B = a_indices[boxId*3+1];
      const uint32_t C = a_indices[boxId*3+2];
    
      const float4 v1  = a_vertices[A];
      const float4 v2  = a_vertices[B];
      const float4 v3  = a_vertices[C];
    
      const float4 triBoxMin = min(v1, min(v2, v3));
      const float4 triBoxMax = max(v1, max(v2, v3));

      bmin = min(bmin, triBoxMin);
      bmax = max(bmax, triBoxMax);
    }
    
    uint32_t leafOffset = i + newIndices.size() - 1;
    a_nodes[leafOffset] = make_BVHNode(to_float3(bmin), to_float3(bmax), PackOffsetAndSize(start, count), 0xFFFFFFFF);
  }
}    

int LBVHBuilder::delta(int i, int j)
{
  if (j < 0 || j >= int(compressedCodes.size()))
    return -1;
  uint32_t code1 = compressedCodes[i];
  uint32_t code2 = compressedCodes[j];
  uint32_t x = code1 ^ code2;
  if (x == 0)
    return 32;
  
  uint32_t result2 = 0;
  while ((x & 0x80000000) == 0)
  {
    x <<= 1;
    result2++;
  }
  return int(result2);
}

void LBVHBuilder::kernel1D_Karras12(BVHNode* a_nodes)
{     
  #pragma omp parallel for schedule(static)
  for (int thread_id = 0; thread_id < m_leafNumberMinusOne; ++thread_id) // int(compressedCodes.size() - 1)
  {
    int d = sign(delta(thread_id, thread_id + 1) - delta(thread_id, thread_id - 1));
    int delta_min = delta(thread_id, thread_id - d);
    int lmax = 2;
    while (delta(thread_id, thread_id + lmax * d) > delta_min)
    {
      lmax <<= 1;
    }
    int l = 0;
    for (int t = lmax >> 1; t > 0; t >>= 1)
    {
      if (delta(thread_id, thread_id + (l + t) * d) > delta_min)
      {
        l += t;
      }
    }
    int j = thread_id + l * d;
    int delta_node = delta(thread_id, j);
    int s = 0;
    for (int t = (l + 1) >> 1; t > 0; t = (t == 1 ? 0 : (t + 1) >> 1))
    {
      if (delta(thread_id, thread_id + (s + t) * d) > delta_node)
      {
        s += t;
      }
    }
    int gamma = thread_id + s * d + std::min(d, 0);
    bool leftIsLeaf  = (std::min(thread_id, j) == gamma);
    bool rightIsLeaf = (std::max(thread_id, j) == gamma + 1);
    a_nodes[thread_id].leftOffset  = (leftIsLeaf  ? compressedCodes.size() - 1 + gamma : gamma);
    a_nodes[thread_id].escapeIndex = (rightIsLeaf ? compressedCodes.size() - 1 + gamma + 1 : gamma + 1);
  }
}

void LBVHBuilder::kernel1D_RefitInit(const BVHNode* a_nodes, uint2* a_leftRight)
{
  for(int idx = 0; idx < m_leafNumberMinusOne; idx++) // int(newIndices.size()) - 1
  {
    const BVHNode currNode = a_nodes[idx];
    a_leftRight[idx] = uint2(currNode.leftOffset, currNode.escapeIndex);
    m_tempCount[idx] = 2;
  }
}

void LBVHBuilder::kernel1D_RefitPass(BVHNode* a_nodes, uint2* a_leftRight)
{
  for(int idx = 0; idx < m_leafNumberMinusOne; idx++) // int(newIndices.size()) - 1
  {
    const uint2 lr = a_leftRight[idx];      
    if(lr.x != 0xFFFFFFFF)
    {
      uint2 lr_left  = uint2(0xFFFFFFFF, 0xFFFFFFFF); //
      uint2 lr_right = uint2(0xFFFFFFFF, 0xFFFFFFFF); //
      if(int(lr.x) < int(newIndices.size()) - 1)
        lr_left  =  a_leftRight[lr.x];
      if(int(lr.y) < int(newIndices.size()) - 1) 
        lr_right = a_leftRight[lr.y];
        
      if(lr_left.x == 0xFFFFFFFF && lr_right.x == 0xFFFFFFFF) 
      {
        const BVHNode leftNode  = a_nodes[lr.x];
        const BVHNode rightNode = a_nodes[lr.y];
        const uint32_t leftCount  = m_tempCount[lr.x];
        const uint32_t rightCount = m_tempCount[lr.y];
  
        BVHNode currNode;
        currNode.boxMin = min(leftNode.boxMin, rightNode.boxMin);
        currNode.boxMax = max(leftNode.boxMax, rightNode.boxMax);
        currNode.leftOffset  = lr.x;
        currNode.escapeIndex = lr.y;
      
        a_nodes    [idx] = currNode;
        a_leftRight[idx] = uint2(0xFFFFFFFF, 0xFFFFFFFF);
        m_tempCount[idx] = leftCount + rightCount + 2;
        
        if(idx == 0)  // we have finished the root, stop other passes
          newIndices.resize(0); 
      }
    }
  } 
}


void LBVHBuilder::kernel1D_PeachConvertInit(const BVHNode* a_inNodes, BVHNode* a_nodes, uint32_t a_size)
{
  for(uint32_t i = 0; i < a_size; i++)
  {
    BVHNode root  = a_inNodes[0];
    if(int(root.leftOffset) < 0)
    {
      root.escapeIndex   = uint32_t(-2);
      m_tempAddresses[0] = uint2(0xFFFFFFFF, 0xFFFFFFFF);
      m_tempAddresses[1] = uint2(0xFFFFFFFF, 0xFFFFFFFF);
      a_nodes[0] = root;
      a_nodes[1] = DummyNode();
    }
    else
    {
      const BVHNode left  = a_inNodes[root.leftOffset];
      const BVHNode right = a_inNodes[root.escapeIndex];
      m_tempAddresses[0] = uint2(left.leftOffset,  left.escapeIndex);
      m_tempAddresses[1] = uint2(right.leftOffset, right.escapeIndex);
      
      const uint2 leftLE  = m_leftAndEscape[0];
      const uint2 rightLE = m_leftAndEscape[1];
      BVHNode leftNode  = a_nodes[0];
      BVHNode rightNode = a_nodes[1];

      leftNode.boxMin  = left.boxMin;
      leftNode.boxMax  = left.boxMax;
      rightNode.boxMin = right.boxMin;
      rightNode.boxMax = right.boxMax;

      leftNode.leftOffset   = leftLE.x;
      leftNode.escapeIndex  = leftLE.y;      
      rightNode.leftOffset  = rightLE.x;
      rightNode.escapeIndex = rightLE.y;  

      if(int(left.leftOffset) < 0) // leaf  // 
      {
        leftNode.leftOffset = left.leftOffset;
        m_tempAddresses[0] = uint2(0xFFFFFFFF, 0xFFFFFFFF);
      }    

      if(int(right.leftOffset) < 0) // leaf //
      {
        rightNode.leftOffset = right.leftOffset;
        m_tempAddresses[1]   = uint2(0xFFFFFFFF, 0xFFFFFFFF);
      } 
  
      a_nodes[0] = leftNode;
      a_nodes[1] = rightNode;
    }
  }
}   

void LBVHBuilder::kernel1D_PeachConvertPass(const BVHNode* a_inNodes, uint32_t currLevelStart, uint32_t currLevelSize, uint32_t nextLevelStart, uint32_t nextLevelSize, uint32_t a_lastPass,
                                            BVHNode* a_nodes)
{
  for(uint32_t i = 0; i < currLevelSize; i++)
  {
    const uint2 lr = m_tempAddresses[currLevelStart + i];
    
    bool leftIsLeaf  = true;
    bool rightIsLeaf = true;

    if(lr.x == 0xFFFFFFFF)
    {
      a_nodes        [nextLevelStart + 2*i + 0] = DummyNode();
      m_tempAddresses[nextLevelStart + 2*i + 0] = uint2(0xFFFFFFFF, 0xFFFFFFFF);
    }
    else
    {
      const BVHNode left = a_inNodes[lr.x];
      const uint2 leftLE = m_leftAndEscape[nextLevelStart + 2*i+0];
      leftIsLeaf         = (int(left.leftOffset) < 0);
      
      
      BVHNode leftNode;
      leftNode.boxMin      = left.boxMin;
      leftNode.boxMax      = left.boxMax;
      leftNode.leftOffset  = leftIsLeaf ? left.leftOffset : leftLE.x;
      leftNode.escapeIndex = leftLE.y;    
        
      a_nodes        [nextLevelStart + 2*i + 0] = leftNode;
      m_tempAddresses[nextLevelStart + 2*i + 0] = leftIsLeaf ? uint2(0xFFFFFFFF, 0xFFFFFFFF) : uint2(left.leftOffset, left.escapeIndex);
    }

    if(lr.y == 0xFFFFFFFF)
    {
      a_nodes        [nextLevelStart + 2*i + 1] = DummyNode();
      m_tempAddresses[nextLevelStart + 2*i + 1] = uint2(0xFFFFFFFF, 0xFFFFFFFF);
    }
    else
    {
      const BVHNode right = a_inNodes[lr.y];
      const uint2 rightLE = m_leftAndEscape[nextLevelStart + 2*i+1];
      rightIsLeaf         = (int(right.leftOffset) < 0); 

      BVHNode rightNode;
      rightNode.boxMin      = right.boxMin;
      rightNode.boxMax      = right.boxMax; 
      rightNode.leftOffset  = rightIsLeaf ? right.leftOffset : rightLE.x;
      rightNode.escapeIndex = rightLE.y;   

      a_nodes        [nextLevelStart + 2*i + 1] = rightNode;
      m_tempAddresses[nextLevelStart + 2*i + 1] = rightIsLeaf ? uint2(0xFFFFFFFF, 0xFFFFFFFF) : uint2(right.leftOffset, right.escapeIndex);
    }

    if(a_lastPass == 1)
    {
      uint32_t countLeft = 0;
      uint32_t countRight = 0;

      if(!leftIsLeaf)
        countLeft = m_tempCount[lr.x];
      if(!rightIsLeaf)
        countRight = m_tempCount[lr.y];

      m_tempCountPrefix[2*i+0] = countLeft;
      m_tempCountPrefix[2*i+1] = countRight;
    }
  }
}      

uint32_t LBVHBuilder::Visit(const BVHNode* a_inNodes, BVHNode* a_nodes, uint32_t& a_treeOffset, 
                            uint32_t a_leftOffset, uint32_t a_rightOffset, uint32_t a_parentEscapeIndex) //
{
  if(int(a_leftOffset) < 0)
    return a_leftOffset;

  const BVHNode leftNode  = a_inNodes[a_leftOffset];
  const BVHNode rightNode = a_inNodes[a_rightOffset];

  //const uint32_t totalNodesCount = m_tempCount[a_leftOffset] + m_tempCount[a_rightOffset] + 2;

  const uint32_t currSize = a_treeOffset;
  a_treeOffset += 2;

  a_nodes[currSize + 0] = leftNode;
  a_nodes[currSize + 1] = rightNode;

  a_nodes[currSize + 0].escapeIndex = currSize + 1;
  a_nodes[currSize + 1].escapeIndex = a_parentEscapeIndex;

  if(int(leftNode.leftOffset) > 0)
    a_nodes[currSize + 0].leftOffset = Visit(a_inNodes, a_nodes, a_treeOffset, leftNode.leftOffset, leftNode.escapeIndex, currSize + 1);
  
  if(int(rightNode.leftOffset) > 0)
    a_nodes[currSize + 1].leftOffset = Visit(a_inNodes, a_nodes, a_treeOffset, rightNode.leftOffset, rightNode.escapeIndex, a_parentEscapeIndex);
  
  return currSize;
}

constexpr uint32_t CONVERT_STACK_SIZE = 32;
void LBVHBuilder::kernel1D_RecursiveConvert(const BVHNode* a_inNodes, const uint32_t* a_subtreeOffsets,  uint32_t currLevelStart, uint32_t currLevelSize, 
                                            BVHNode* a_nodes)
{
  for(uint32_t i = 0; i < currLevelSize; i++)
  { 
    const uint2 lr = m_tempAddresses[currLevelStart + i];
    uint32_t subtreeOffset = currLevelStart + currLevelSize;
    if(i > 0)
      subtreeOffset += a_subtreeOffsets[i-1];
    
    if(lr.x != 0xFFFFFFFF) 
    {
      a_nodes[currLevelStart + i].leftOffset = subtreeOffset;
  
      uint4 stack[CONVERT_STACK_SIZE];
      int   top = 0;
      
      uint32_t leftOffset        = lr.x;
      uint32_t rightOffset       = lr.y;
      uint32_t parentEscapeIndex = a_nodes[currLevelStart + i].escapeIndex;
      
      while (top >= 0) 
      { 
        bool needStackPop = (int(leftOffset) < 0);
       
        if(int(leftOffset) > 0)
        {
          const BVHNode leftNode  = a_inNodes[leftOffset];
          const BVHNode rightNode = a_inNodes[rightOffset];  
          const uint32_t currSize = subtreeOffset;
          subtreeOffset += 2;
        
          a_nodes[currSize + 0] = leftNode;
          a_nodes[currSize + 1] = rightNode;
        
          a_nodes[currSize + 0].escapeIndex = currSize + 1;
          a_nodes[currSize + 1].escapeIndex = parentEscapeIndex;
          
          needStackPop = (int(leftNode.leftOffset) < 0) && (int(rightNode.leftOffset) < 0);
  
          if(int(leftNode.leftOffset) > 0 && int(rightNode.leftOffset) > 0)
          {
            stack[top] = uint4(rightNode.leftOffset, rightNode.escapeIndex, parentEscapeIndex, currSize+1);
            top++;
          }

          if(int(leftNode.leftOffset) > 0)
          {
            leftOffset        = leftNode.leftOffset; 
            rightOffset       = leftNode.escapeIndex;
            parentEscapeIndex = currSize + 1;
            a_nodes[currSize + 0].leftOffset = subtreeOffset;
          }
          else if(int(rightNode.leftOffset) > 0)
          {
            leftOffset        = rightNode.leftOffset; 
            rightOffset       = rightNode.escapeIndex;
            // parentEscapeIndex unchanged
            a_nodes[currSize + 1].leftOffset = subtreeOffset;
          }  
        }
  
        if(needStackPop)
        {
          top--;
          uint4 oldData     = stack[std::max(top,0)];

          leftOffset        = oldData.x;
          rightOffset       = oldData.y;
          parentEscapeIndex = oldData.z;
  
          if(int(leftOffset) > 0 && top >= 0)
            a_nodes[oldData.w].leftOffset = subtreeOffset;
        }
      
      }
    }
    
    //if(lr.x != 0xFFFFFFFF)
    //  a_nodes[currLevelStart + i].leftOffset = Visit(a_inNodes, a_nodes, subtreeOffset, lr.x, lr.y, a_nodes[currLevelStart + i].escapeIndex); 
   
    //uint32_t nextSubtreeOffset = (i < currLevelSize-1) ? a_subtreeOffsets[i] + currLevelStart + currLevelSize : m_reservedNodes;
    //assert(subtreeOffset <= nextSubtreeOffset);
  }

}


void LBVHBuilder::BuildFromTriangles(const float4* a_vertOrBoxes, uint32_t a_vertNum, 
                                     const uint32_t* a_indices,   uint32_t a_indexNum,
                                     BVHNode* a_outNodes,         uint32_t a_outMaxNodes, 
                                     uint32_t* a_outIndices)
{
  const uint32_t a_boxCount       = uint32_t(a_indexNum/3);
  const uint32_t boxesNumExtended = uint32_t(powerOfTwo(a_boxCount));

  //Timer timer;  
  //timer.start();

  // (1) input ==> runCodesAndIndices
  //
  //kernel1D_RedunctionFromTriangles(a_vertOrBoxes, a_vertNum, a_indices, a_boxCount);
  kernel1D_RedunctionFromVertices(a_vertOrBoxes, a_vertNum);
  kernel1D_EvalMortonCodesFromTriangles(a_vertOrBoxes, uint32_t(a_vertNum), a_indices, a_boxCount, boxesNumExtended, runCodesAndIndices.data());
  //m_timings[0] = timer.getElaspedMs();
  
  //return; // DEBUG_BUILDER

  // (2) sort runCodesAndIndices by morton code
  //
  //timer.start();
  std::sort(runCodesAndIndices.begin(), runCodesAndIndices.end(), [](uint2 a, uint2 b) { return a.x < b.x; });
  //m_timings[1] = timer.getElaspedMs();

  // (3) make leafes array
  //  
  //timer.start();
  
  kernel1D_AppendInit(runCodesAndIndices.data(), a_boxCount, // runCodesAndIndices ==> (codesEq, prefixCodeEq)
                      codesEq.data(), prefixCodeEq.data()); 
  
  std::exclusive_scan(prefixCodeEq.begin(), prefixCodeEq.end(), prefixCodeEq.begin(), 0);
    
  kernel1D_AppendComplete(codesEq.data(), prefixCodeEq.data(), a_boxCount); // (codesEq, prefixCodeEq) ==> (newIndices, compressedCodes)
  
  if(m_treeFormat == uint32_t(cbvh2::BVH2_LEFT_RIGHT))
  {
    kernel1D_MakeLeavesFromTriangles(a_vertOrBoxes, a_vertNum, a_indices, a_indexNum, a_outNodes, a_outIndices);
  
    //m_timings[2] = timer.getElaspedMs();
    //timer.start();
    kernel1D_Karras12(a_outNodes);                                // compressedCodes ==> a_outNodes
    //m_timings[3] = timer.getElaspedMs();

    //return; // DEBUG_BUILDER

    //timer.start();
    kernel1D_RefitInit(a_outNodes, runCodesAndIndices.data());    // a_outNodes ==> runCodesAndIndices (used as temp buffer)
    for(int pass = 0; pass < 32; pass++)
    {
      kernel1D_RefitPass(a_outNodes, runCodesAndIndices.data());  // (a_outNodes,runCodesAndIndices) ==> (a_outNodes,runCodesAndIndices)
    }
    //m_timings[4] = timer.getElaspedMs();
  }
  else
  {
    kernel1D_MakeLeavesFromTriangles(a_vertOrBoxes, a_vertNum, a_indices, a_indexNum, m_tempNodes.data(), a_outIndices);
  
    //m_timings[2] = timer.getElaspedMs();
    //timer.start();
    kernel1D_Karras12(m_tempNodes.data());                                // compressedCodes ==> a_outNodes
    //m_timings[3] = timer.getElaspedMs();                
  
    //timer.start();
    kernel1D_RefitInit(m_tempNodes.data(), runCodesAndIndices.data());    // m_tempNodes ==> runCodesAndIndices (used as temp buffer)
    for(int pass = 0; pass < 32; pass++)
    {
      kernel1D_RefitPass(m_tempNodes.data(), runCodesAndIndices.data());  // (m_tempNodes,runCodesAndIndices) ==> (m_tempNodes,runCodesAndIndices)
    }
    //m_timings[4] = timer.getElaspedMs();
    
    //timer.start();
    kernel1D_PeachConvertInit(m_tempNodes.data(), a_outNodes, 1);
    for(size_t i = 0; i < m_peachIntervals.size()-1; i++)
    {
      const uint2 currLevel = m_peachIntervals[i+0];
      const uint2 nextLevel = m_peachIntervals[i+1];
      const uint32_t last   = (i == m_peachIntervals.size()-2) ? 1 : 0; 
      kernel1D_PeachConvertPass(m_tempNodes.data(), currLevel.x, currLevel.y, nextLevel.x, nextLevel.y, last,
                                a_outNodes);
    }
    
    std::inclusive_scan(m_tempCountPrefix.begin(), m_tempCountPrefix.end(), m_tempCountPrefix.begin(), std::plus<uint32_t>(), 0);

    if(m_peachIntervals.size() != 1)
    {
      const uint2 lastLevel = m_peachIntervals.back();
      kernel1D_RecursiveConvert(m_tempNodes.data(), m_tempCountPrefix.data(), lastLevel.x, lastLevel.y, 
                                a_outNodes);
    }
    //m_timings[5] = timer.getElaspedMs();
  }

}

void LBVHBuilder::BuildFromBoxes(const float4* a_boxes,    uint32_t a_boxNum,
                                 BVHNode*      a_outNodes, uint32_t a_outMaxNodes, 
                                 uint32_t*     a_outIndices)
{
  const uint32_t a_boxCount = a_boxNum;
  const uint32_t boxesNumExtended = uint32_t(powerOfTwo(a_boxCount));

  //Timer timer;  
  //timer.start();

  // (1) input ==> runCodesAndIndices
  //
  kernel1D_RedunctionFromBoxes(a_boxes, uint32_t(a_boxCount));
  kernel1D_EvalMortonCodesFromBoxes(a_boxes, uint32_t(a_boxCount), boxesNumExtended, runCodesAndIndices.data());
  //m_timings[0] = timer.getElaspedMs();
  
  // (2) sort runCodesAndIndices by morton code
  //
  //timer.start();
  std::sort(runCodesAndIndices.begin(), runCodesAndIndices.end(), [](uint2 a, uint2 b) { return a.x < b.x; });
  //m_timings[1] = timer.getElaspedMs();

  // (3) make leafes array
  //  
  //timer.start();
  
  kernel1D_AppendInit(runCodesAndIndices.data(), a_boxCount, // runCodesAndIndices ==> (codesEq, prefixCodeEq)
                      codesEq.data(), prefixCodeEq.data()); 
  
  std::exclusive_scan(prefixCodeEq.begin(), prefixCodeEq.end(), prefixCodeEq.begin(), 0);
    
  kernel1D_AppendComplete(codesEq.data(), prefixCodeEq.data(), a_boxCount); // (codesEq, prefixCodeEq) ==> (newIndices, compressedCodes)
  
  if(m_treeFormat == uint32_t(cbvh2::BVH2_LEFT_RIGHT))
  {
    kernel1D_MakeLeavesFromBoxes(a_boxes, a_boxCount, a_outNodes, a_outIndices);
  
    //m_timings[2] = timer.getElaspedMs();
    //timer.start();
    kernel1D_Karras12(a_outNodes);                                // compressedCodes ==> a_outNodes
    //m_timings[3] = timer.getElaspedMs();

    //timer.start();
    kernel1D_RefitInit(a_outNodes, runCodesAndIndices.data());    // a_outNodes ==> runCodesAndIndices (used as temp buffer)
    for(int pass = 0; pass < 32; pass++)
    {
      kernel1D_RefitPass(a_outNodes, runCodesAndIndices.data());  // (a_outNodes,runCodesAndIndices) ==> (a_outNodes,runCodesAndIndices)
    }
    //m_timings[4] = timer.getElaspedMs();
  }
  else
  {
    kernel1D_MakeLeavesFromBoxes(a_boxes, a_boxCount, m_tempNodes.data(), a_outIndices);
  
    //m_timings[2] = timer.getElaspedMs();
    //timer.start();
    kernel1D_Karras12(m_tempNodes.data());                                // compressedCodes ==> a_outNodes
    //m_timings[3] = timer.getElaspedMs();                
  
    //timer.start();
    kernel1D_RefitInit(m_tempNodes.data(), runCodesAndIndices.data());    // m_tempNodes ==> runCodesAndIndices (used as temp buffer)
    for(int pass = 0; pass < 32; pass++)
    {
      kernel1D_RefitPass(m_tempNodes.data(), runCodesAndIndices.data());  // (m_tempNodes,runCodesAndIndices) ==> (m_tempNodes,runCodesAndIndices)
    }
    //m_timings[4] = timer.getElaspedMs();
    
    //timer.start();
    kernel1D_PeachConvertInit(m_tempNodes.data(), a_outNodes, 1);
    for(size_t i = 0; i < m_peachIntervals.size()-1; i++)
    {
      const uint2 currLevel = m_peachIntervals[i+0];
      const uint2 nextLevel = m_peachIntervals[i+1];
      const uint32_t last   = (i == m_peachIntervals.size()-2) ? 1 : 0; 
      kernel1D_PeachConvertPass(m_tempNodes.data(), currLevel.x, currLevel.y, nextLevel.x, nextLevel.y, last,
                                a_outNodes);
    }
    
    std::inclusive_scan(m_tempCountPrefix.begin(), m_tempCountPrefix.end(), m_tempCountPrefix.begin(), std::plus<uint32_t>(), 0);

    if(m_peachIntervals.size() != 1)
    {
      const uint2 lastLevel = m_peachIntervals.back();
      kernel1D_RecursiveConvert(m_tempNodes.data(), m_tempCountPrefix.data(), lastLevel.x, lastLevel.y, 
                                a_outNodes);
    }
    //m_timings[5] = timer.getElaspedMs();
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LBVHBuilder::Refit(BVHNode* a_nodes, int idx)
{
  if (int(a_nodes[idx].leftOffset) > 0)
    Refit(a_nodes, a_nodes[idx].leftOffset);

  if (int(a_nodes[idx].escapeIndex) > 0)
    Refit(a_nodes, a_nodes[idx].escapeIndex);
  
  if(int(a_nodes[idx].leftOffset) > 0 && int(a_nodes[idx].escapeIndex) > 0)
  {
    uint32_t leftOffset  = a_nodes[idx].leftOffset;
    uint32_t rightOffset = a_nodes[idx].escapeIndex;
    a_nodes[idx].boxMin  = min(a_nodes[leftOffset].boxMin, a_nodes[rightOffset].boxMin);
    a_nodes[idx].boxMax  = max(a_nodes[leftOffset].boxMax, a_nodes[rightOffset].boxMax);
  }
}

void LBVHBuilder::GetExecutionTime(const char* a_funcName, float a_out[8])
{
  for(int i=0;i<8;i++)
    a_out[i] = m_timings[i];
}