#pragma once

#include <vector>
#include <limits>
#include <cstdint>
#include <cassert>
#include <memory>
#include <iostream>

#include "cbvh.h"

namespace cbvh
{
  using BVHNode = cbvh2::BVHNode;

  using LiteMath::float3;
  using LiteMath::float4;
  using LiteMath::float4x4;
  using LiteMath::uint4;
  using LiteMath::uint;
  using LiteMath::dot;

  static std::ostream& operator<<(std::ostream& os, const BVHNode& dt)
  {
    os << dt.boxMin.x << "\t" << dt.boxMin.y << "\t" << dt.boxMin.z << "\t|\t" << dt.leftOffset << std::endl;
    os << dt.boxMax.x << "\t" << dt.boxMax.y << "\t" << dt.boxMax.z << "\t|\t" << dt.escapeIndex << std::endl;
    return os;
  }

  constexpr unsigned int LEAF_NORMAL = 0xFFFFFFFF; ///!<
  constexpr unsigned int LEAF_EMPTY  = 0xFFFFFFFD; ///!<
  constexpr unsigned int ESCAPE_ROOT = 0xFFFFFFFE; ///!< when escapeIndex (for stackless traversal) points out of tree

  static inline bool IsLeaf (const BVHNode& a_node) { return (a_node.leftOffset == LEAF_NORMAL); }
  static inline bool IsEmpty(const BVHNode& a_node) { return (a_node.leftOffset == LEAF_EMPTY); }
  static inline bool IsValid(const BVHNode& a_node) { return (a_node.leftOffset <  LEAF_EMPTY); }

  struct Interval
  {
    Interval() : start(0), count(0) {}
    Interval(uint a, uint b) : start(a), count(b) {}

    uint start;
    uint count;
  };

  enum CBVH_FORMATS { FMT_BVH2Node32_Interval32_Static  = 1, //!< Children: (leftOffset, leftOffset+1); escapeIndex is true escapeIndex;
                      FMT_BVH2Node32_Interval32_Dynamic = 2, //!< Children: (leftOffset, rightOffset == escapeIndex); true escapeIndex is not stored;
                      FMT_BVH4Node32_Interval32_Static  = 3, //!< Children: (leftOffset, leftOffset+1,leftOffset+2,leftOffset+3); escapeIndex is true escapeIndex;

                      FMT_BVH_UNKNOWN                   = -1, //!< if passed to builder via 'BVHPresets', builder will select some format according to other settings
                      };

  /**
  \brief Linearized layout bvh tree
  */
  struct BVHTree
  {
    BVHTree(){}
    BVHTree(const std::vector<BVHNode>&  a_nodes, 
            const std::vector<Interval>& a_intervals,
            const std::vector<uint32_t>& a_indices) : nodes(a_nodes), intervals(a_intervals), indicesReordered(a_indices) {}


    std::vector<BVHNode>  nodes;            //!< nodes stored strongly in bread-first order
    std::vector<Interval> intervals;        //!< 
    std::vector<uint32_t> indicesReordered; //!< reordered index buffer. may not used for TLAS and custom user objects, i.e. indicesReordered.size() will be equal to zero
                                            //!< in this case, reordering of objects in memory is not happened
                                            //!< for triangle meshes indicesReordered.size() == numTris*3, for custom objects indicesReordered.size() == a_objNum or zero (if reordering didn't happened).

    std::vector<Interval> depthRanges;      //!< store begin/size pair for nodes offsets per depth level (i.e. depthRanges[1] in the first leve and depthRanges[2] ins the second)
    unsigned int          geomID;           //!< identifier that builder writes, but user could overwrite it if wants

    CBVH_FORMATS          format;           //!< actual format according to specification (you need to treat 'escapeIndex' field differently for 'FMT_BVH4Node32_Interval32_Static' format)
    unsigned int          leavesNumber = 0; //!< filled by the LBVH builder only

    void ComputeDepthRanges();              //!< user don't have to call this function in general 
    void Print(std::ostream& out);          //!< for debug needs

    bool invalidIntervals = false;          //!< for internal usage
  };

  enum BVH_BUILDER_TYPE{ BVH_CONSTRUCT_QUALITY     = 0, 
                         BVH_CONSTRUCT_MEDIUM      = 1, 
                         BVH_CONSTRUCT_FAST        = 2,
                         BVH_CONSTRUCT_FAST_GPU    = 3,
                         BVH_CONSTRUCT_EMBREE2     = 4,
                         BVH_CONSTRUCT_EMBREE_FAST = 5};

  struct BVHPresets
  {
    int  primsInLeaf  = 4;                             ///<! recomended primitives in leaf
    int  childrenNum  = 4;                             ///<! (2,4,8,16,32) 
    int  maxThreads   = 0;                             ///<! 0 means unbounden number of threads; please note that current implementation is single threaded!
    BVH_BUILDER_TYPE btype = BVH_CONSTRUCT_MEDIUM;     ///<! builder type actually;
    CBVH_FORMATS     desiredFormat = FMT_BVH4Node32_Interval32_Static; ///<! desired format of BVH; if required format is not possiable for current presets (not compatible with childrenNum or not implemented for example), different format will be selected.
    bool useEmbreeIfPossiable = true;                  ///<! switch to embree if it was linked via CMake
    bool alwaysNotEmpty       = true;
    bool enableSpatialSplit   = false;
  };

  /**
  \brief Main Builder C++ interface, build BLAS (Bottom Level Acceleration Structure) for input triangle mesh.

  \param  a_vertices - input triangle vertices; MUST BE ALIGNED(16)
  \param  a_vertNum  - input verices count
  \param  a_indices  - input index buffer
  \param  a_indexNum - input indices number
  \param  a_presets  - input builder presets

  \return linearized layout bvh tree
  */
  BVHTree BuildBVH(const float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, BVHPresets a_presets = BVHPresets());

  /**
  \brief Main Builder C++ interface, Build both TLAS (Top Level Acceleration Structure) and BLAS for custom user input.

  \param  a_nodes  - input bounding boxes with 2 custom fields (leftOffset and escapeIndex) wich are unused in general
  \param  a_objNum - input bounding boxes number
  \param  a_presets  - input builder presets

  \return linearized layout bvh tree

  User must pack his custom data inside two integer fields of BVHNode: (leftOffset and escapeIndex).
  These fields are ignored the builder. 

  */
  BVHTree BuildBVH(const BVHNode* a_nodes, size_t a_objNum, BVHPresets a_presets = BVHPresets());  
}

using cbvh::IsEmpty;
using cbvh::IsLeaf;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

