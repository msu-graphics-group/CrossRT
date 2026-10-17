#include "refitter.h"
#include <cassert>

using LiteMath::float3;
using LiteMath::float4;
using LiteMath::to_float3;
using LiteMath::min;
using LiteMath::max;

void Refitter::getBoundingBoxFromLeaf(uint32_t a_start, uint32_t a_count, const uint32_t* a_indices, uint32_t a_indexNum, 
                                      const float *a_vpos3f, uint32_t a_vertNumber, const uint32_t *a_triIndices, uint32_t a_indNumber, uint32_t vFloatStride, 
                                      float3* pBoxMin, float3* pBoxMax)
{
  
  *pBoxMin = float3(+REFIT_INF, +REFIT_INF, +REFIT_INF);
  *pBoxMax = float3(-REFIT_INF, -REFIT_INF, -REFIT_INF);

  for (uint32_t triId = a_start; triId < a_start + a_count; triId++)
  {
    const uint32_t T = triId; //a_indices[triId];
    const uint32_t A = a_triIndices[T*3 + 0];
    const uint32_t B = a_triIndices[T*3 + 1];
    const uint32_t C = a_triIndices[T*3 + 2];

    const float3 A_pos = float3(a_vpos3f[A*vFloatStride + 0], a_vpos3f[A*vFloatStride + 1], a_vpos3f[A*vFloatStride + 2]);
    const float3 B_pos = float3(a_vpos3f[B*vFloatStride + 0], a_vpos3f[B*vFloatStride + 1], a_vpos3f[B*vFloatStride + 2]);
    const float3 C_pos = float3(a_vpos3f[C*vFloatStride + 0], a_vpos3f[C*vFloatStride + 1], a_vpos3f[C*vFloatStride + 2]);

    *pBoxMin = min(min(A_pos, *pBoxMin), min(B_pos, C_pos));
    *pBoxMax = max(max(A_pos, *pBoxMax), max(B_pos, C_pos));
  }
}

void Refitter::kernel1D_UpdateNodesFat32(BVHNodeFat* a_nodes, uint32_t a_start, uint32_t a_size, const uint32_t* a_indices, uint32_t a_indexNum, 
                                         const float *a_vpos3f, uint32_t a_vertNumber, const uint32_t *a_triIndices, uint32_t a_indNumber, uint32_t vFloatStride)
{
  for(uint32_t i = a_start; i < a_start + a_size; i++)
  {
    BVHNodeFat currNode = a_nodes[i];
    
    // recalculate left box
    {
      float3 leftBoxMin = to_float3(currNode.lmin_xyz_rmax_x);
      float3 leftBoxMax = to_float3(currNode.lmax_xyz_rmax_y);

      if((currNode.offs_left & LEAF_BIT) != 0)
      {
        const uint32_t start = EXTRACT_START(currNode.offs_left);
        const uint32_t count = EXTRACT_COUNT(currNode.offs_left);
        getBoundingBoxFromLeaf(start, count, a_indices, a_indexNum,
                               a_vpos3f, a_vertNumber, a_triIndices, a_indNumber, vFloatStride,
                               &leftBoxMin, &leftBoxMax);
      }
      else
      {
        const BVHNodeFat leftNode = a_nodes[currNode.offs_left];
        leftBoxMin = min( to_float3(leftNode.lmin_xyz_rmax_x), to_float3(leftNode.rmin_xyz_rmax_z));
        leftBoxMax = max( to_float3(leftNode.lmax_xyz_rmax_y), float3(leftNode.lmin_xyz_rmax_x.w, leftNode.lmax_xyz_rmax_y.w, leftNode.rmin_xyz_rmax_z.w));
      }

      assert(std::abs(currNode.lmin_xyz_rmax_x.x - leftBoxMin.x) < 1e-6f); // for debugging only when refit is used immediately after build
      assert(std::abs(currNode.lmin_xyz_rmax_x.y - leftBoxMin.y) < 1e-6f);
      assert(std::abs(currNode.lmin_xyz_rmax_x.z - leftBoxMin.z) < 1e-6f);

      assert(std::abs(currNode.lmax_xyz_rmax_y.x - leftBoxMax.x) < 1e-6f);
      assert(std::abs(currNode.lmax_xyz_rmax_y.y - leftBoxMax.y) < 1e-6f);
      assert(std::abs(currNode.lmax_xyz_rmax_y.z - leftBoxMax.z) < 1e-6f);

      currNode.lmin_xyz_rmax_x.x = leftBoxMin.x;
      currNode.lmin_xyz_rmax_x.y = leftBoxMin.y;
      currNode.lmin_xyz_rmax_x.z = leftBoxMin.z;

      currNode.lmax_xyz_rmax_y.x = leftBoxMax.x;
      currNode.lmax_xyz_rmax_y.y = leftBoxMax.y;
      currNode.lmax_xyz_rmax_y.z = leftBoxMax.z;
    }

    // recalculate right box
    {
      float3 rightBoxMin = to_float3(currNode.rmin_xyz_rmax_z);
      float3 rightBoxMax = float3(currNode.lmin_xyz_rmax_x.w, currNode.lmax_xyz_rmax_y.w, currNode.rmin_xyz_rmax_z.w);

      if((currNode.offs_right & LEAF_BIT) != 0)
      {
        const uint32_t start = EXTRACT_START(currNode.offs_right);
        const uint32_t count = EXTRACT_COUNT(currNode.offs_right);
        getBoundingBoxFromLeaf(start, count, a_indices, a_indexNum,
                               a_vpos3f, a_vertNumber, a_triIndices, a_indNumber, vFloatStride,
                               &rightBoxMin, &rightBoxMax);
      }
      else
      {
        const BVHNodeFat rightNode = a_nodes[currNode.offs_right];
        rightBoxMin = min( to_float3(rightNode.lmin_xyz_rmax_x), to_float3(rightNode.rmin_xyz_rmax_z));
        rightBoxMax = max( to_float3(rightNode.lmax_xyz_rmax_y), float3(rightNode.lmin_xyz_rmax_x.w, rightNode.lmax_xyz_rmax_y.w, rightNode.rmin_xyz_rmax_z.w));
      }
      
      if(EXTRACT_COUNT(currNode.offs_right) != 0)
      {
        assert(std::abs(currNode.rmin_xyz_rmax_z.x - rightBoxMin.x) < 1e-6f); // for debugging only when refit is used immediately after build
        assert(std::abs(currNode.rmin_xyz_rmax_z.y - rightBoxMin.y) < 1e-6f); 
        assert(std::abs(currNode.rmin_xyz_rmax_z.z - rightBoxMin.z) < 1e-6f); 
  
        assert(std::abs(currNode.lmin_xyz_rmax_x.w - rightBoxMax.x) < 1e-6f); 
        assert(std::abs(currNode.lmax_xyz_rmax_y.w - rightBoxMax.y) < 1e-6f); 
        assert(std::abs(currNode.rmin_xyz_rmax_z.w - rightBoxMax.z) < 1e-6f); 
      }
      
      currNode.rmin_xyz_rmax_z.x = rightBoxMin.x;
      currNode.rmin_xyz_rmax_z.y = rightBoxMin.y;
      currNode.rmin_xyz_rmax_z.z = rightBoxMin.z;

      currNode.lmin_xyz_rmax_x.w = rightBoxMax.x;
      currNode.lmax_xyz_rmax_y.w = rightBoxMax.y;
      currNode.rmin_xyz_rmax_z.w = rightBoxMax.z;
    }

    a_nodes[i] = currNode;
  }
}

void Refitter::Refit_BVH2Fat32(BVHNodeFat* a_nodes, uint32_t a_nodes_num, 
                               const uint32_t* a_indices, uint32_t a_indexNum, 
                               const uint2* a_depthRanges, uint32_t a_depthRangesNum,
                               const float *a_vpos3f, uint32_t a_vertNumber, const uint32_t *a_triIndices, uint32_t a_indNumber, uint32_t vFloatStride)
{
  // (1) update leaves only; this may not work for BVH2Fat format because one node could be leaf and the other is not ... 
  //
  
  // (2) update down-top and leaves simultaniously
  //
  const uint2 lastInterval = a_depthRanges[a_depthRangesNum-1];
  for(int i = int(a_depthRangesNum-1); i >=0; i--)
  {
    const uint2 currInterval = a_depthRanges[a_depthRangesNum-1];
    kernel1D_UpdateNodesFat32(a_nodes, currInterval.x, currInterval.y, a_indices, a_indexNum, 
                              a_vpos3f, a_vertNumber, a_triIndices, a_indNumber, vFloatStride);
  }
}
