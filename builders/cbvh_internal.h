#ifndef CBVH_CORE_INTERNAL_H
#define CBVH_CORE_INTERNAL_H

#include "cbvh_core.h"

namespace cbvh
{ 
  void Init();
  void Destroy();

  /**
  \brief Build BLAS (Bottom Level Acceleration Structure) for input triangle mesh.

  \param  a_vertices - input triangle vertices
  \param  a_vertNum  - input verices count
  \param  a_indices  - input index buffer
  \param  a_indexNum - input indices number

  \return linearized layout bvh tree
  */
  BVHTree BuildBVH4(const float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum);


  /**
  \brief Build both TLAS (Top Level Acceleration Structure) and BLAS for custom user input.

  \param  a_nodes  - input bounding boxes with 2 custom fields
  \param  a_objNum - input bounding boxes number

  \return linearized layout bvh tree

  User must pack his custom data inside two integer fields of BVHNode: (leftOffset and escapeIndex).
  These fields are ignored the builder. 

  */
  BVHTree BuildBVH4(const BVHNode* a_nodes, size_t a_objNum);

  BVHTree BuildLBVH(const float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, CBVH_FORMATS a_desiredFormat = FMT_BVH4Node32_Interval32_Static, bool onGPU = false);

  /**
  \brief This is slow but reference way to convert 'BVH2Dynamic' format that we get from LBVH builder to 'BVH4Flat' (or so called 'BVH4Static') 
  */
  BVHTree ConvertBVH2DynamicToBVH4Flat(const BVHTree& a_tree);
  BVHTree ConvertBVH2DynamicToBVH2Flat(const BVHTree& a_tree);
};

namespace cbvh_internal
{
  using LiteMath::float3;
  using LiteMath::float4;
  using LiteMath::float4x4;
  using LiteMath::uint4;
  using LiteMath::uint;

  using cbvh::BVHTree;
  using cbvh::BVHNode;
  using cbvh::Interval;

  using LiteMath::Box4f;

  struct Triangle4f
  {
    float4 A; // as_int(A.w) may store index of A also
    float4 B; // as_int(B.w) may store index of B also
    float4 C; // as_int(C.w) may store index of C also
    //uint32_t primId = uint32_t(-1); // used by SplitBVH builder, not used by ESC
  };


  inline float SurfaceArea(float3 in_boxMin, float3 in_boxMax){
    const float3 abc = in_boxMax - in_boxMin;
    return 2.0f*(abc[0]*abc[1] + abc[0]*abc[2] + abc[1]*abc[2]);
  }

  inline float SurfaceArea(float4 in_boxMin, float4 in_boxMax){
    const float4 abc = in_boxMax - in_boxMin;
    return 2.0f*(abc[0]*abc[1] + abc[0]*abc[2] + abc[1]*abc[2]);
  }

  struct ESC_Settings
  {
    float expandFactor   = 1.5f;   ///<! we can not expand total number of primitives more than expandFactor*N
    int   maxSubdivs     = 4;      ///<! each cluster should not be split greater than this value
    bool  enableClusters = true;   ///<! do not change it please, for debug and comparison purposes only
    bool  prebuildForESC = false;  ///<! internal parameter, don't use it please 
    float thresholdMax   = 25.0f;  ///<! internal parameter, don't change it please 
    float thresholdMin   = 10.0f;  ///<! internal parameter, don't change it please
    int   prebuildPIL    = 16;     ///<! internal parameter, don't change it please
  }; 

  struct ESC_Result
  {
    std::vector<Box4f> boxes;       ///<! as_int(boxes[i].boxMin.w)*3 + 0 -- index start; as_int(boxes[i].boxMax.w)*3 -- index count;
    std::vector<uint>  indexBuffer; ///<! reordered index buffer
  };


  /**
  \brief Perform Early Split Clipping (ESC).

  \param  a_vertices - input triangle mesh vertices
  \param  a_vertNum  - input mesh vertex number
  \param  a_indices  - input index buffer 
  \param  a_indexNum - input size of index buffer
  \param  a_settings - input algorithm settings

  \return Splitted bounding boxes of triangles. Long bad triangles are represented with several boxes.

  */
  ESC_Result ESC(const LiteMath::float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, ESC_Settings a_settings = ESC_Settings(), bool a_debug = false);
  ESC_Result ESC2(const LiteMath::float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, ESC_Settings a_settings = ESC_Settings(), bool a_debug = false);

  struct ESC_Clusters
  {
    std::vector<Box4f> clustersForSplit;
    std::vector<Box4f> theRestOfTriangles;
    
    std::vector<uint32_t> indReord;
    Box4f sceneBox;
    float avgTriSize;
  };
  
  ESC_Clusters MakeTriangleClustersForESC(const LiteMath::float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, const ESC_Settings& a_settings);
  
  struct TriCluster
  {
    TriCluster()
    {
      boxMin = float4( +std::numeric_limits<float>::infinity() );
      boxMax = float4( -std::numeric_limits<float>::infinity() ); 
      triList.reserve(16);  
      escMetric = 0.0f;
    }

    inline bool AxisAligned(int axis, float split) const { return (boxMin[axis] == boxMax[axis]) && (boxMin[axis]==split); }

    inline void include(const LiteMath::float4 p) 
    {
      boxMin = LiteMath::min(boxMin, p);
      boxMax = LiteMath::max(boxMax, p);
    } 

    inline void intersect(const TriCluster& a_box) 
    {
      boxMin = LiteMath::max(boxMin, a_box.boxMin);
      boxMax = LiteMath::min(boxMax, a_box.boxMax);
    }

    float4 boxMin;
    float4 boxMax;
    std::vector<Triangle4f> triList;
    float escMetric; 
  };
  TriCluster ReadTrianglesFromPrimBox(const Box4f& prim, const LiteMath::float4* a_vertices, size_t a_vertNum, 
                                      const uint* a_indices, size_t a_indexNum,
                                      const uint* a_indicesOld, size_t a_indexNumOld);

  void SplitTriangles(const cbvh_internal::TriCluster& nodeBox, float splitPos, int splitAxis,
                      cbvh_internal::TriCluster& leftBox, cbvh_internal::TriCluster& rightBox);

  enum SPLIT_TYPE {OBJECT_SPLIT_SAH = 1, OBJECT_SPLIT_EPO = 2};

  struct SplitInfo
  {
    SplitInfo() : splitPos(0), splitAxis(0), metricVal(+std::numeric_limits<float>::infinity()) {}
  
    Box4f  boxL;
    Box4f  boxR;
    size_t splitPos;
    int    splitAxis;
    float  metricVal;
  };

  /**
  \brief Find best object split of input bounding boxes using SAH or EOP. 

  \param  in_boxes       - input primitive boxes
  \param  in_bbox        - input bounding box of all primitive boxes
  \param  in_splitType   - input algorithm type (SAH or EPO)
  \param  out_boxesLeft  - output left part of subdivided boxes array
  \param  out_boxesRight - output right part of subdividex boxes array
  \return                - output other subdivide info
  */
  SplitInfo FindObjectSplit(const std::vector<Box4f>& in_boxes, const Box4f& in_bbox, SPLIT_TYPE in_splitType,
                            std::vector<Box4f>& out_boxesLeft, std::vector<Box4f>& out_boxesRight);

  struct BVHSettings
  {
    int  primsInLeaf  = 4;     ///<! recomended primitives in leaf
    int  childrenNum  = 4;     ///<! (2,4,8,16,32) 
    cbvh::CBVH_FORMATS desiredFormat = cbvh::FMT_BVH2Node32_Interval32_Static;
    bool enableESC      = false; ///<! enable Early Split Clipping
    bool enableSS       = false; ///<! enable spatial split
    bool discardLines   = true;  ///<! discard triangles that are lines actually
    bool perfMeasure    = false; ///<! internal parameter, don't use it!
    bool alwaysNotEmpty = true;  ///<! 
    int  qualityLevel   = 1;     ///<!
    ESC_Settings esc;            ///<! Early Split Clipping settings
  };

  struct OriginalInput
  {
    OriginalInput() : vertices(nullptr), indices(nullptr), vertNum(0), indexNum(0) {}
    OriginalInput(const float4* a_vertices, size_t a_vertNum, const uint32_t* a_indices, size_t a_indexNum) : vertices(a_vertices), indices(a_indices), vertNum(a_vertNum), indexNum(a_indexNum) {}
    const float4*   vertices = nullptr;
    const uint32_t* indices  = nullptr;
    size_t vertNum;
    size_t indexNum;
  };

  /**
  \brief Build BVH

  \param  a_boxes    - input boxes (constructed over user objects); #MANDATORY_INPUT: box[i].boxMin.w = primIndex; box[i].boxMax.w = primsCount;
  \param  a_objNum   - input boxes number
  \param  a_indices  - input optional index buffer (for triangles meshes)
  \param  a_indexNum - input optional index buffer size
  \param  a_settings - input builder settings
  \param  a_args     - original user input for triangle data if there were actually triangles
  \return            - constructed bvh tree
  */
  BVHTree BuildBVH4(const Box4f* a_boxes, size_t a_objNum, const uint32_t* a_indices = nullptr, size_t a_indexNum = 0, BVHSettings a_settings = BVHSettings(), OriginalInput a_args = OriginalInput());

  /**
  \brief Build BLAS (Bottom Level Acceleration Structure) for input triangle mesh.

  \param  a_vertices - input triangle vertices
  \param  a_vertNum  - input verices count
  \param  a_indices  - input index buffer
  \param  a_indexNum - input indices number

  \return linearized layout bvh tree
  */
  BVHTree BuildBVH4(const float4* a_vertices, size_t a_vertNum, const uint32_t* indices, size_t a_indexNum, BVHSettings a_settings = BVHSettings());


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  struct TraversalRes
  {
    TraversalRes() : numBoxTests(0), numPrimTests(0), numLeafesTests(0), primId(0), t(1e38f) {} 
    int   numBoxTests;
    int   numPrimTests;
    int   numLeafesTests;
    int   primId;
    float t;
  };

  TraversalRes BVH4RayTraversalExample(float4 rayPosAndNear, float4 rayDirAndFar,
                                       const cbvh::BVHTree& a_tree, const float4* a_vertices);

  TraversalRes BVH4RayTraversalExample(float4 rayPosAndNear, float4 rayDirAndFar,
                                       const cbvh::BVHTree& a_tree, const float4* a_vertices, std::vector<uint32_t>& a_leaves);

  struct RayTriangleHit
  {
    RayTriangleHit() : t(std::numeric_limits<float>::infinity()), primId(uint(-1)) {}
  
    float t;
    uint  primId;
  };

  RayTriangleHit IntersectAllPrimitivesInLeaf(const float4 rayPosAndNear, const float4 rayDirAndFar,
                                              const uint* a_indices, uint a_start, uint a_count, const float4* a_vert);

  using Ray4f = LiteMath::Ray4f;
  typedef std::vector<Ray4f> RaysArray;

  static inline float RelativeBoxOverlap(const Box4f& box1, const Box4f& box2) //#TODO: opt this with transpose4 funct
  {
    const Box4f ovlp = BoxBoxOverlap(box1, box2);
  
    const float sa1 = box1.surfaceArea();
    const float sa2 = box2.surfaceArea();
    const float sa3 = ovlp.surfaceArea();
  
    return std::min(sa3/std::min(sa1,sa2), 1.0f);
  }

  static constexpr uint32_t SPLIT_ANYWAY_THRESHOLD = 64;
  void DebugPrintBVH(const char* a_outFileName, const cbvh::BVHTree a_treeData);
};

#endif
