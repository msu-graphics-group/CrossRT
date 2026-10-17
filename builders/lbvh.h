#pragma once
#include "cbvh.h"

#include <cstdint>
#include <chrono>

using LiteMath::float4;
using LiteMath::float3;
using LiteMath::uint4;
using LiteMath::uint3;
using LiteMath::uint2;
using cbvh2::BVHNode;

struct LBVHBuilder
{
  LBVHBuilder() { }
  virtual ~LBVHBuilder(){}

  virtual void CommitDeviceData() {}                                     // will be overriden in generated class
  virtual void GetExecutionTime(const char* a_funcName, float a_out[8]); // will be overriden in generated class

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  constexpr static uint32_t PEACHES_MIN = 256;  // min number of nodes in perfect tree
  constexpr static uint32_t PEACHES_MAX = 2048; // max number of nodes in perfect tree

  /**
   * \brief reserve memory inside builder
   * 
   * \param a_inPrimsNum -- input number of boxes (a_vertNum if a_indexNum == 0) or triangles (a_indexNum/3)
   * \param a_treeFormat -- output tree format, currently supported only 'BVH2_LEFT_RIGHT' 
   */
  virtual size_t Reserve(size_t a_inPrimsNum, cbvh2::FMT a_treeFormat = cbvh2::BVH2_LEFT_RIGHT);
  
  /**
   * \brief Add geometry of type 'Triangles' to 'internal geometry library' of scene object and return geometry id
   * 
   * \param a_vertOrBoxes -- input vertex 
   * \param a_vertNumber  -- input vertex 
   * \param a_triIndices  -- input triangle indices (standart index buffer)
   * \param a_indNumber   -- input number of indices, should be equal to 3*triaglesNum in your mesh; 
   * 
   * \param a_outNodes    -- output bvh nodes; should be at least of size == 2*primsNum-1 // where prims num == a_indexNum/3
   * \param a_outIndices  -- output indices;   should be at least of size == a_indexNum 
   */
  #ifdef KERNEL_SLICER
  virtual void BuildFromTriangles(const float4*   a_vertOrBoxes __attribute__((size("a_vertNum"))),     uint32_t a_vertNum, 
                                  const uint32_t* a_indices     __attribute__((size("a_indexNum"))),    uint32_t a_indexNum,
                                  BVHNode*        a_outNodes    __attribute__((size("a_outMaxNodes"))), uint32_t a_outMaxNodes, 
                                  uint32_t*       a_outIndices  __attribute__((size("a_indexNum/3"))));
  #else
  virtual void BuildFromTriangles(const float4*   a_vertOrBoxes, uint32_t a_vertNum, 
                                  const uint32_t* a_indices    , uint32_t a_indexNum,
                                  BVHNode*        a_outNodes   , uint32_t a_outMaxNodes, 
                                  uint32_t*       a_outIndices);
  #endif

  /**
   * \brief Add geometry of type 'Triangles' to 'internal geometry library' of scene object and return geometry id
   * 
   * \param a_vertOrBoxes -- input boxes data, of total bytesize sizeof(float4)*2*a_vertNumber;
   * \param a_vertNumber  -- input boxes number;
   * 
   * \param a_outNodes    -- output bvh nodes; should be at least of size == 2*primsNum-1 // where prims num == 'a_vertNumber' for boxes or a_indexNum/3 for indices
   * \param a_outIndices  -- output indices;   should be at least of size == a_indexNum 
   */

  #ifdef KERNEL_SLICER
  virtual void BuildFromBoxes(const float4* a_boxes      __attribute__((size("a_boxNum*2"))),    uint32_t a_boxNum, // a_boxNum*2 because one box is two float4
                              BVHNode*      a_outNodes   __attribute__((size("a_outMaxNodes"))), uint32_t a_outMaxNodes, 
                              uint32_t*     a_outIndices __attribute__((size("a_boxNum"))));
  #else
  virtual void BuildFromBoxes(const float4* a_boxes,    uint32_t a_boxNum,
                              BVHNode*      a_outNodes, uint32_t a_outMaxNodes, 
                              uint32_t*     a_outIndices);
  #endif

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  virtual void kernel1D_RedunctionFromBoxes(const float4* a_vertices, uint32_t a_size);
  virtual void kernel1D_RedunctionFromTriangles(const float4* a_vertices, uint32_t a_vertNum, const uint32_t* a_indices, uint32_t a_trisNum);
  virtual void kernel1D_RedunctionFromVertices(const float4* a_vertices, uint32_t a_vertNum);

  virtual void kernel1D_EvalMortonCodesFromBoxes(const float4* a_inBoxes, uint32_t a_boxCount, uint32_t a_powerOf2Size,
                                                 uint2* outRunCodes);

  virtual void kernel1D_EvalMortonCodesFromTriangles(const float4* a_vertices, uint32_t a_vertNum, const uint32_t* a_indices, uint32_t a_trisNum, uint32_t a_powerOf2Size,
                                                     uint2* outRunCodes);

  virtual void kernel1D_AppendInit(const uint2* in_runCodes, uint32_t a_boxCount, 
                                   uint32_t* a_codesEq, uint32_t* a_codesPref);

  virtual void kernel1D_AppendComplete(uint32_t* a_codesEq, uint32_t* a_codesPref, uint32_t a_boxCount);

  virtual void kernel1D_MakeLeavesFromBoxes(const float4* a_inBoxes, uint32_t a_boxCount, 
                                            BVHNode* a_nodes, uint32_t* a_outIndices);

  virtual void kernel1D_MakeLeavesFromTriangles(const float4* a_vertices, uint32_t a_vertNum, const uint32_t* a_indices, uint32_t a_indexNum,
                                                BVHNode* a_nodes, uint32_t* a_outIndices); 

  int delta(int i, int j);
  virtual void kernel1D_Karras12(BVHNode* a_nodes);

  virtual void kernel1D_RefitInit(const BVHNode* a_nodes, uint2* a_leftRight);
  virtual void kernel1D_RefitPass(BVHNode* a_nodes, uint2* a_leftRight);


  virtual void kernel1D_PeachConvertInit(const BVHNode* a_inNodes, BVHNode* a_nodes, uint32_t a_size);
  virtual void kernel1D_PeachConvertPass(const BVHNode* a_inNodes, uint32_t currLevelStart, uint32_t currLevelSize, uint32_t nextLevelStart, uint32_t nextLevelSize, uint32_t a_lastPass,
                                         BVHNode* a_nodes);
  
  virtual void kernel1D_RecursiveConvert(const BVHNode* a_inNodes, const uint32_t* a_subtreeOffsets,  uint32_t currLevelStart, uint32_t currLevelSize, 
                                         BVHNode* a_nodes);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  float4 bboxMin;
  float4 bboxMax;
  float4 zIndexScale;

  uint32_t m_treeFormat;
  uint32_t m_reservedNodes;
  int      m_leafNumber;
  int      m_leafNumberMinusOne;

  inline uint32_t GetLeavesNumber() const { return uint32_t(m_leafNumber); }
  inline uint32_t GetNodesNumber()  const { return uint32_t(std::max(m_leafNumber*2-1,0)); }

  std::vector<uint2>    runCodesAndIndices;
  std::vector<uint32_t> codesEq;
  std::vector<uint32_t> prefixCodeEq;
  std::vector<uint32_t> compressedCodes;
  std::vector<uint32_t> newIndices; 
  
  ///////////////////////////////////////// this is needed only for tree format conversion (a_treeFormat == BVH2_LEFT_OFFSET or BVH4_LEFT_OFFSET)
  std::vector<BVHNode>  m_tempNodes;
  std::vector<uint2>    m_tempAddresses;
  std::vector<uint32_t> m_tempCount;
  std::vector<uint32_t> m_tempCountPrefix;

  std::vector<uint2>    m_peachIntervals;   ///<! start/count for Perfet Implicit Trees ("Peach Trees") 
  std::vector<uint2>    m_leftAndEscape;
  ///////////////////////////////////////// this is needed only for tree format conversion (a_treeFormat == BVH2_LEFT_OFFSET or BVH4_LEFT_OFFSET)
 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void BuildPeachTree(const std::vector<uint2>& a_peachIntervals, std::vector<uint2>& a_leftAndEscape);
  void Refit(BVHNode* nodes, int idx); // recursive version, for CPU only

  uint32_t Visit(const BVHNode* a_inNodes, BVHNode* a_nodes, uint32_t& a_treeOffset, 
                 uint32_t a_leftOffset, uint32_t a_rightOffset, uint32_t a_parentEscapeIndex);

  static size_t roundUpTo(size_t a_size, size_t a_num);
  static size_t powerOfTwo(size_t a_size); 

  float m_timings[8];           
};

//class Timer 
//{
//  std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
//public:
//  Timer() { start(); }
//  void start();
//  float getElaspedMs() const;
//};