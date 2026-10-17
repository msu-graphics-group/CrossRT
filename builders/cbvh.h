#pragma once

#include <vector>
#include <cstdint>
#include <cassert>

#include "LiteMath.h"

namespace cbvh2
{
  using LiteMath::float4;
  using LiteMath::float3;
  using LiteMath::uint;
  using LiteMath::uint2;
  using LiteMath::uint4;
  #ifdef HALFFLOAT
  using LiteMath::half4;
  #endif
  
  struct BVHNode 
  {
    float3 boxMin;
    uint   leftOffset; //!< please note that when LEAF_BIT (0x80000000) is set in leftOffset, this node is a leaf
    float3 boxMax;
    uint   escapeIndex;
  };

  static inline BVHNode DummyNode()
  {
    BVHNode dummyNode;
    dummyNode.boxMin      = float3(0,0,0);
    dummyNode.boxMax      = float3(0,0,0);
    dummyNode.leftOffset  = 0x80000000; // or 0xFFFFFFFD according to old cbvh convention
    dummyNode.escapeIndex = 0xFFFFFFFF; // or 0xFFFFFFFD according to old cbvh convention
    return dummyNode;
  }

  static inline BVHNode make_BVHNode(float3 in_boxMin, float3 in_boxMax, uint in_leftOffset, uint in_escapeIndex)
  {
    BVHNode res;
    res.boxMin      = in_boxMin;
    res.leftOffset  = in_leftOffset;
    res.boxMax      = in_boxMax;
    res.escapeIndex = in_escapeIndex;
    return res;
  }

  struct BVHNodeFat
  { 
    float4 lmin_xyz_rmax_x;
    float4 lmax_xyz_rmax_y;
    float4 rmin_xyz_rmax_z; 
  
    uint32_t offs_left;
    uint32_t offs_right;
    uint32_t dummy1;
    uint32_t dummy2;
  };
  
  #ifdef HALFFLOAT
  struct BVHNodeFat16
  { 
    half4 lmin_xyz_rmax_x;
    half4 lmax_xyz_rmax_y;
    half4 rmin_xyz_rmax_z; 
  
    uint32_t offs_left;
    uint32_t offs_right;
  };
  #endif

  enum FMT{ BVH2_LEFT_OFFSET = 1, //!< Children: (leftOffset, leftOffset+1); 'escapeIndex' is a true escapeIndex;
            BVH2_LEFT_RIGHT  = 2, //!< Children: (leftOffset, rightOffset is 'escapeIndex'); actual escapeIndex is not stored;
            BVH4_LEFT_OFFSET = 3, //!< Children: (leftOffset, leftOffset+1,leftOffset+2,leftOffset+3); escapeIndex is a true escapeIndex;
            BVH2_LEFT_ROPES  = 4, //!< Children: (leftOffset, unknown); 'escapeIndex' is a true escapeIndex; This format is ONLY for stackless traversal
            };

  enum BTYPE{ BVH_CONSTRUCT_QUALITY     = 0, ///<! Restricted SBVH + Sweep builder
              BVH_CONSTRUCT_MEDIUM      = 1, ///<! Sweep builder
              BVH_CONSTRUCT_FAST        = 2, ///<! LBVH and e.t.c (CPU, in early versions was GPU by default)
              BVH_CONSTRUCT_FAST_GPU    = 3, ///<! LBVH and e.t.c (GPU)
              BVH_CONSTRUCT_EMBREE      = 4, ///<! implementation via embree, HQ variant only
              BVH_CONSTRUCT_EMBREE_FAST = 5, ///<! implementation via embree, HQ variant only
              BVH_CONSTRUCT_NANORT      = 6, ///<! implementation via nanort, ~ same as sweep builder
              };

  enum TLAYOUT{ LAYOUT_DFS         = 0, ///<! Depth First
                LAYOUT_ODFS        = 1, ///<! Ordered Depth First
                LAYOUT_BFS         = 2, ///<! Breadth First
                LAYOUT_COLBVH_TRB1 = 3, ///<! Treelet: analogue to [AK10] 
                LAYOUT_COLBVH_TRB2 = 4, ///<! Treelets Super: [AK10]+[YM06]
                LAYOUT_COLBVH_TRB3 = 5, ///<! Treelets Super Super: [AK10]+[YM06]
                LAYOUT_CALBVH      = 6, ///<! our CALBVH implementation, please note that 'BuildBVHFat' function don't compress bvh itself!
                LAYOUT_COLBVH_YM06 = 7, ///<! [YM06] clusters approach
              };

  #ifndef KERNEL_SLICER

  struct BuilderPresets
  {
    FMT   format      = BVH2_LEFT_OFFSET;
    BTYPE quality     = BVH_CONSTRUCT_QUALITY;  
    int   primsInLeaf = 2;                      ///<! recomended primitives in leaf
  };

  struct LayoutPresets
  {
    TLAYOUT layout = LAYOUT_COLBVH_TRB3; ///<! 
    int     grSzId = 1;                  ///<! Group Size Id; 0,1,2,3,4; Please note this is just id, rel group size could be evaluated via 'CalcActualGroupSize' function
    int     grSzXX = 0;                  ///<! Group Size Internal; This in internal parameter, don't use it please!
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  struct BVHTreeCommon
  {
    BVHTreeCommon(){}
    BVHTreeCommon(const std::vector<BVHNode>& a_nodes, const std::vector<uint32_t>& a_indices) : nodes(a_nodes), indices(a_indices) { } 
    std::vector<BVHNode>  nodes;
    std::vector<uint32_t> indices;
  };

  struct BVHTreeFat
  {
    BVHTreeFat(){}
    BVHTreeFat(const std::vector<BVHNodeFat>& a_nodes, const std::vector<uint32_t>& a_indices) : nodes(a_nodes), indices(a_indices) { } 
    BVHTreeFat(const std::vector<BVHNodeFat>& a_nodes, const std::vector<uint32_t>& a_indices, const std::vector<uint2>& a_depthRanges) : nodes(a_nodes), indices(a_indices), depthRanges(a_depthRanges) { } 
    
    std::vector<BVHNodeFat>  nodes;
    std::vector<uint32_t>    indices;
    std::vector<uint2>       depthRanges;  ///!< store begin/size pair for nodes offsets per depth level (i.e. depthRanges[1] in the first level and depthRanges[2] its the second); valid only for BreadthFirst layouts
  };
  
  #ifdef HALFFLOAT
  struct BVHTreeFat16
  {
    BVHTreeFat16(){}
    BVHTreeFat16(const std::vector<BVHNodeFat16>& a_nodes, const std::vector<uint32_t>& a_indices) : nodes(a_nodes), indices(a_indices) { } 
    std::vector<BVHNodeFat16>  nodes;
    std::vector<uint32_t>      indices;
  };
  #endif

  struct BVHTreeFatCompressed
  {
    BVHTreeFatCompressed(){}
    BVHTreeFatCompressed(const std::vector<uint4>& a_nodes, const std::vector<uint32_t>& a_indices, int a_log2) : nodes(a_nodes), indices(a_indices), log2_group_size(a_log2) { } 
    std::vector<uint4>    nodes;   // stricly speak, not 'nodes' because some elements are support bounding boxes, not nodes
    std::vector<uint32_t> indices;
    int log2_group_size;           // if log2_group_size == 0, threat single uint4 as single half float node, else => this is CALBVH with treelets
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
  \brief Main Builder C++ interface, build BLAS (Bottom Level Acceleration Structure) for input triangle mesh.

  \param  a_vertices    - input triangle vertices; 
  \param  a_vertNum     - input verices count
  \param  a_vByteStride - input stride in bytes between vertices in 'a_vpos3f' (12 for float3, 16 flor float4 and e.t.c.)

  \param  a_indices     - input index buffer
  \param  a_indexNum    - input indices number
  \param  a_presets     - input builder presets

  \return linearized layout bvh tree in two arrays 
  */
  BVHTreeCommon BuildBVH(const float* a_vpos3f,     size_t a_vertNum, size_t a_vByteStride, 
                         const uint32_t* a_indices, size_t a_indexNum, BuilderPresets a_presets = BuilderPresets());

  /**
  \brief Main Builder C++ interface, Build both TLAS (Top Level Acceleration Structure) and BLAS for custom user input.

  \param  a_nodes   - input bounding boxes with 2 custom fields (leftOffset and escapeIndex) wich are unused in general
  \param  a_objNum  - input bounding boxes number
  \param  a_presets - input builder presets

  \return linearized layout bvh tree

  User may pack his custom data inside two integer fields of BVHNode: (leftOffset and escapeIndex).
  These fields are ignored by the builder. 

  */
  std::vector<BVHNode> BuildBVH(const BVHNode* a_nodes, size_t a_objNum, BuilderPresets a_presets = BuilderPresets());
  
  BuilderPresets BuilderPresetsFromString(const char* a_str);
  LayoutPresets  LayoutPresetsFromString(const char* a_str);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  /**
  \brief Main Builder C++ interface, build BLAS (Bottom Level Acceleration Structure) for input triangle mesh.

  \param  a_vertices    - input triangle vertices; 
  \param  a_vertNum     - input verices count
  \param  a_vByteStride - input stride in bytes between vertices in 'a_vpos3f' (12 for float3, 16 flor float4 and e.t.c.)

  \param  a_indices     - input index buffer
  \param  a_indexNum    - input indices number
  \param  a_presets     - input builder presets
  \param  a_layout      - input tree layout

  \param  pTreelet_root       - output optional internal parameter; don't use it! you may look inside 'FatBVH::PrintForGraphViz' to understrand it
  \param  pSuper_treelet_root - output optional internal parameter; don't use it! you may look inside 'FatBVH::PrintForGraphViz' to understrand it
  \param  pGrStart            - output optional internal parameter; don't use it! you may look inside 'FatBVH::PrintForGraphViz' to understrand it

  \return linearized layout bvh tree in two arrays 
  */
  BVHTreeFat BuildBVHFat(const float* a_vpos3f,     size_t a_vertNum, size_t a_vByteStride, 
                         const uint32_t* a_indices, size_t a_indexNum, BuilderPresets a_presets = BuilderPresets(), LayoutPresets a_layout = LayoutPresets(),
                         std::vector<int>* pTreelet_root = nullptr, std::vector<int>* pSuper_treelet_root = nullptr, std::vector<int>* pGrStart = nullptr);
  
  #ifdef HALFFLOAT
  /**
  \brief Main Builder C++ interface, build BLAS (Bottom Level Acceleration Structure) for input triangle mesh.

  \param  a_vertices    - input triangle vertices; 
  \param  a_vertNum     - input verices count
  \param  a_vByteStride - input stride in bytes between vertices in 'a_vpos3f' (12 for float3, 16 flor float4 and e.t.c.)

  \param  a_indices     - input index buffer
  \param  a_indexNum    - input indices number
  \param  a_presets     - input builder presets
  \param  a_layout      - input tree layout

  \return linearized layout bvh tree in two arrays 
  */

  BVHTreeFat16 BuildBVHFat16(const float* a_vpos3f,     size_t a_vertNum, size_t a_vByteStride, const uint32_t* a_indices, size_t a_indexNum, 
                             BuilderPresets a_presets = BuilderPresets(), LayoutPresets a_layout = LayoutPresets());
  #endif

  /**
  \brief Main Builder C++ interface, build BLAS (Bottom Level Acceleration Structure) for input triangle mesh.

  \param  a_vertices    - input triangle vertices; 
  \param  a_vertNum     - input verices count
  \param  a_vByteStride - input stride in bytes between vertices in 'a_vpos3f' (12 for float3, 16 flor float4 and e.t.c.)

  \param  a_indices     - input index buffer
  \param  a_indexNum    - input indices number
  \param  a_presets     - input builder presets
  \param  a_layout      - input tree layout

  \return linearizedcand compressed layout bvh tree in two arrays 
  */
  BVHTreeFatCompressed BuildBVHFatCompressed(const float* a_vpos3f,     size_t a_vertNum, size_t a_vByteStride, 
                                             const uint32_t* a_indices, size_t a_indexNum, BuilderPresets a_presets = BuilderPresets(), LayoutPresets a_layout = LayoutPresets());
  


  #ifdef HALFFLOAT
  BVHTreeFatCompressed BuildBVHHalf(const float* a_vpos3f,     size_t a_vertNum, size_t a_vByteStride, 
                                    const uint32_t* a_indices, size_t a_indexNum, BuilderPresets a_presets = BuilderPresets());
  #endif

  int CalcInputGroupSize(const LayoutPresets& a_layout);
  int CalcActualGroupSize(const LayoutPresets& a_layout);

  // todo: add conversion from FatBVH to CommonBVH and CommonBVH Builder with fixed layout. I.e. conversion should be plain, with fixed addrs.

  #endif
};

