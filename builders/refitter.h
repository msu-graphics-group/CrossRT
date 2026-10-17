#pragma once
#include "cbvh.h"
#include "raytrace_common.h"

using cbvh2::BVHNodeFat;
using LiteMath::uint2;

class Refitter
{
public:
  Refitter(){}
  virtual ~Refitter(){}

  virtual void Refit_BVH2Fat32(BVHNodeFat* a_nodes, uint32_t a_nodes_num, 
                               const uint32_t* a_indices, uint32_t a_indexNum, 
                               const uint2* a_depthRanges, uint32_t a_depthRangesNum,
                               const float *a_vpos3f, uint32_t a_vertNumber, const uint32_t *a_triIndices, uint32_t a_indNumber, uint32_t vFloatStride);


protected:

  static constexpr float REFIT_INF = +1e10;

  void getBoundingBoxFromLeaf(uint32_t a_start, uint32_t a_count, const uint32_t* a_indices, uint32_t a_indexNum, 
                              const float *a_vpos3f, uint32_t a_vertNumber, const uint32_t *a_triIndices, uint32_t a_indNumber, uint32_t vFloatStride, 
                              float3* pBoxMin, float3* pBoxMax);

  virtual void kernel1D_UpdateNodesFat32(BVHNodeFat* a_nodes, uint32_t a_start, uint32_t a_size, const uint32_t* a_indices, uint32_t a_indexNum, 
                                         const float *a_vpos3f, uint32_t a_vertNumber, const uint32_t *a_triIndices, uint32_t a_indNumber, uint32_t vFloatStride);

};
