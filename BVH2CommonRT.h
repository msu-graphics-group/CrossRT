#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <string>

#include "CrossRT.h"
#include "raytrace_common.h"

using LiteMath::cross;
using LiteMath::dot;
using LiteMath::float2;
using LiteMath::float3;
using LiteMath::float4;
using LiteMath::float4x4;
using LiteMath::int2;
using LiteMath::inverse4x4;
using LiteMath::normalize;
using LiteMath::sign;
using LiteMath::to_float3;
using LiteMath::uint2;
using LiteMath::Box4f;

//#include "aligned_alloc.h"
#include "builders/cbvh.h"
using cbvh2::BVHNode;

struct BVH2CommonRT : public ISceneObject
{
  BVH2CommonRT(const char* a_builderName = "cbvh_embree2") : m_builderName(a_builderName) {}
  ~BVH2CommonRT() override {}

  const char* Name() const override { return "BVH2_LOFT"; }

  void ClearGeom() override;

  uint32_t AddGeom_Triangles3f(const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, uint32_t a_qualityLevel, size_t vByteStride) override;
  void     UpdateGeom_Triangles3f(uint32_t a_geomId, const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, uint32_t a_qualityLevel, size_t vByteStride) override;
  
  uint32_t AddGeom_AABB(uint32_t a_typeId, const CRT_AABB* boxMinMaxF8, size_t a_boxNumber, uint32_t a_buildFlags, void** a_customPrimPtrs, size_t a_customPrimCount) override;
  uint32_t AddCustomGeom_FromFile(const char *geom_type_name, const char *filename, ISceneObject *fake_this) override;


  void ClearScene() override;
  void CommitScene(uint32_t a_qualityLevel) override;

  uint32_t AddInstance(uint32_t a_geomId, const float4x4 &a_matrix) override;
  void     UpdateInstance(uint32_t a_instanceId, const float4x4 &a_matrix) override;

  CRT_Hit RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar) override;
  bool    RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar) override;
  
  // (!) not implemented properly
  //
  uint32_t AddInstanceMotion(uint32_t a_geomId, const LiteMath::float4x4* a_matrices, uint32_t a_matrixNumber) override { return AddInstance(a_geomId, a_matrices[0]); }
  CRT_Hit  RayQuery_NearestHitMotion(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, float time) override { return RayQuery_NearestHit(posAndNear, dirAndFar); }
  bool     RayQuery_AnyHitMotion(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, float time) override { return RayQuery_AnyHit(posAndNear, dirAndFar);  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  const Box4f* GetGeomBoxes() override { return m_geomBoxes.data(); }            /// array  of added geometry objects
  uint32_t     GetGeomNum()   override { return uint32_t(m_geomBoxes.size()); }  /// number of added geometry objects
  
  const Box4f* GetInstBoxes() override { return m_instBoxes.data(); }           /// array  of instances
  uint32_t     GetInstNum()   override { return uint32_t(m_instBoxes.size()); } /// number of instances

  CRT_AABB  GetSceneBoundingBox() const override { return m_sceneBBox; }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  #ifndef KERNEL_SLICER
  using UpdateFuncPtr = void (*)(size_t, size_t);
  UpdateFuncPtr m_pUpdate_m_allNodes = nullptr;
  UpdateFuncPtr m_pUpdate_m_instMatricesInv = nullptr;
  virtual void Update_m_allNodes(size_t a_first, size_t a_offet) {
    if(m_pUpdate_m_allNodes != nullptr)
      m_pUpdate_m_allNodes(a_first, a_offet);
  }
  virtual void Update_m_instMatricesInv(size_t a_first, size_t a_offet) {
    if(m_pUpdate_m_instMatricesInv != nullptr)
      m_pUpdate_m_instMatricesInv(a_first, a_offet);
  }
  #endif
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//protected:
  
  static constexpr unsigned int GEOM_TP_MASK = 0xf0000000; // 16   types of possible geometry
  static constexpr unsigned int GEOM_ID_MASK = 0x0fffffff; // 2^28 max geometric objects is allowed which is high enough
  static constexpr unsigned int GEOM_ID_SHFT = 28;                               

  static constexpr unsigned int GEOM_TYPE_TRIANGLE   = 0; //
  static constexpr unsigned int GEOM_TYPE_SPHERE     = 1;
  static constexpr unsigned int GEOM_TYPE_MANDELBULB = 2; 

  CRT_AABB m_sceneBBox;

  uint32_t IntersectAllPrimitivesInLeaf(float4 rayPosAndNear, float4 rayDirAndFar, CRT_LeafInfo info, CRT_Hit *pHit);  
  void     IntersectionComplete(float4 rayPosAndNear, float4 rayDirAndFar, CRT_Hit *pHit);                             

  virtual size_t AppendTreeData(const std::vector<BVHNode>& a_nodes, const std::vector<uint32_t>& a_indices, 
                                const uint32_t *a_triIndices, size_t a_indNumber);
  
  std::vector<Box4f>    m_geomBoxes; // use boxMin.w as geomSize and boxMin.w as geomTag
  std::vector<Box4f>    m_instBoxes;

  std::vector<float4x4> m_instMatricesInv; ///< inverse instance matrices
  std::vector<float4x4> m_instMatricesFwd; ///< instance matrices

  std::vector<float4>   m_vertPos;
  std::vector<uint32_t> m_indices;
  std::vector<uint32_t> m_primIndices;

  //std::vector<BVHNode, aligned<BVHNode, 64> >  m_allNodes;
  std::vector<BVHNode>                         m_allNodes;
  std::vector<uint32_t>                        m_bvhOffsets;
  uint32_t                                     m_tlasOffset    = 0;
  size_t                                       m_tlasOffsetFirst = 0; // when we call CommitScene() first time, assume BLAS would not be rebuild!

  std::vector<uint2>    m_geomOffsets;
  std::vector<uint32_t> m_geomIdByInstId;

  std::string m_builderName;
  size_t totalTrisMem     = 0;
  size_t totalTrisInst    = 0;
  size_t totalinstNumInst = 0;
};
