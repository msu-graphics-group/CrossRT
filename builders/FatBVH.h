#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>

#include "LiteMath.h"
//#include "aligned_alloc.h"

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
using LiteMath::uint;
using LiteMath::uint2;
using LiteMath::uint4;
using LiteMath::Box4f;

#include "CrossRT.h"
#include "raytrace_common.h"

#include "builders/cbvh.h"
#include "builders/cbvh_core.h"

using cbvh::Interval;
using cbvh2::BVHNode;
using cbvh2::BVHNodeFat;
using std::isnan;

/// Bounding box of the left child
static inline Box4f GetChildBoxLeft(const BVHNodeFat& node)
{
  Box4f box;
  box.boxMin = float4(node.lmin_xyz_rmax_x.x, node.lmin_xyz_rmax_x.y, node.lmin_xyz_rmax_x.z, 1.0f);
  box.boxMax = float4(node.lmax_xyz_rmax_y.x, node.lmax_xyz_rmax_y.y, node.lmax_xyz_rmax_y.z, 1.0f);
  return box;
}

/// Bounding box of the right child
static inline Box4f GetChildBoxRight(const BVHNodeFat& node)
{
  Box4f box;
  box.boxMin = float4(node.rmin_xyz_rmax_z.x, node.rmin_xyz_rmax_z.y, node.rmin_xyz_rmax_z.z, 1.0f);
  box.boxMax = float4(node.lmin_xyz_rmax_x.w, node.lmax_xyz_rmax_y.w, node.rmin_xyz_rmax_z.w, 1.0f);
  return box;
}

/// Surface area (all 6 sides) of the bounding box of the left child
static inline float GetChildBoxAreaLeft(const BVHNodeFat& node)
{
  // Dimensions of the bounding box
  const LiteMath::float3 size = LiteMath::to_float3(node.lmax_xyz_rmax_y) - LiteMath::to_float3(node.lmin_xyz_rmax_x);
  return 2 * (size.x * size.y + size.x * size.z + size.y * size.z);
}

/// Surface area (all 6 sides) of the bounding box of the right child
static inline float GetChildBoxAreaRight(const BVHNodeFat& node)
{
  // Dimensions of the bounding box
  const LiteMath::float3 size = LiteMath::float3(node.lmin_xyz_rmax_x.w, node.lmax_xyz_rmax_y.w, node.rmin_xyz_rmax_z.w)
                              - LiteMath::to_float3(node.rmin_xyz_rmax_z);

  return 2 * (size.x * size.y + size.x * size.z + size.y * size.z);
}


namespace FatBVH
{
  /// The level (depth) of treelet-ization
  enum TreeletizLevel
    {
    TREELETS           = 1, /// Simple single-level treelets
    SUPERTREELETS      = 2, /// Supertreelts i.e. treelets of treelets
    SUPERSUPERTREELETS = 3, /// Supersupertreelts i.e. supertreelets of supertreelets
    };

  void ReorderDFL(std::vector<BVHNodeFat>& a_nodes);   ///<! Depth First Layout
  void ReorderODFL(std::vector<BVHNodeFat>& a_nodes);  ///<! Ordered Depth First Layout
  void ReorderBFL(std::vector<BVHNodeFat>& a_nodes);   ///<! Breadth First Layout
  std::vector<uint2> ComputeDepthRanges(const std::vector<BVHNodeFat>& nodes); ///<! Compute intervals of data in several tree levels for Breadth First Layout
 
  /// Treelet Based Layout, set 'a_leavesNum' in leaves count -- 2,4,8,16
  void ReorderTRB(std::vector<BVHNodeFat>& a_nodes, std::vector<int>& treelet_root, std::vector<int>& super_root_idx, size_t a_leavesNum, bool a_aligned = false, TreeletizLevel level = TreeletizLevel::TREELETS);
  void ReorderTRBNew(std::vector<BVHNodeFat>& a_nodes, std::vector<int>& treelet_root, std::vector<int>& super_root_idx, std::vector<int>& grStart, size_t a_leavesNum, bool a_aligned = false, bool a_merge = false);

  /// Reorder by recursive clasterisation
  int  ReorderByClusters(std::vector<BVHNodeFat>& a_nodes, std::vector<int>& cluster_roots, int offset, int cluster_size, bool a_aligned, const std::vector<double>& w0 = std::vector<double>());

  void PrintForGraphViz(const std::vector<BVHNodeFat>& a_nodes, const std::vector<int>& a_treeletRoots, const std::vector<int>& a_treeletRootsSuper, const char* a_fileName);

  BVHNodeFat BVHNodeFatDummy();
};
