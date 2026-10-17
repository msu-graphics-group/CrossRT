#include <iostream>
#include <queue>
#include <stack>
#include <limits>
#include <tuple>
#include <cmath>
#include <chrono>

#include "cbvh_core.h"
#include "cbvh_internal.h"

using cbvh_internal::Box4f;
using cbvh_internal::Triangle4f;

using LiteMath::float3;
using LiteMath::float4;
using LiteMath::float4x4;
using LiteMath::uint;
using LiteMath::min;
using LiteMath::max;

#ifdef USE_EMBREE4
#include "embree4/rtcore.h"
#else
#include "embree3/rtcore.h"
#endif

extern double g_buildTime;

namespace embree
{
  static uint32_t g_recommendedPrimsInLeaf = 4;
  static float    g_overSplitThreshold     = 0.25f;
  
  struct UserData
  {
    const Triangle4f* tris = nullptr;
  }; 

  /* This function is called by the builder to signal progress and to
   * report memory consumption. */
  static  bool memoryMonitor(void* userPtr, ssize_t bytes, bool post) {
    return true;
  }

  static  bool buildProgress (void* userPtr, double f) {
    return true;
  }

  static void splitPrimitive (const RTCBuildPrimitive* prim, unsigned int dim, float pos, RTCBounds* lprim, RTCBounds* rprim, void* userPtr)
  {
    assert(dim < 3);
    assert(prim->geomID == 0);
  
    *(Box4f*)lprim = *(Box4f*) prim;
    *(Box4f*)rprim = *(Box4f*) prim;

    if(false) // stupid implementation of spatial split
    {
      (&lprim->upper_x)[dim] = pos;
      (&rprim->lower_x)[dim] = pos;
    }
    else     // smart implementation of spatial split
    {
      // grow both boxes with vertices and edge-plane intersections
      //
      const auto splitAxis = dim;
      const auto splitPos  = pos;

      Box4f leftBox  = Box4f();
      Box4f rightBox = Box4f();
      
      const UserData* data = (const UserData*)userPtr;
      const auto& tri = data->tris[prim->primID];

      const float4 edges[3][2] = { {tri.A, tri.B},
                                   {tri.C, tri.A},
                                   {tri.B, tri.C} };
      
      for (int i=0;i<3;i++)
      {
        const float v0p = edges[i][0][splitAxis];
        const float v1p = edges[i][1][splitAxis];
        
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
      
      Box4f nodeBox = Box4f(float4(prim->lower_x, prim->lower_y, prim->lower_z, 0.0f), 
                            float4(prim->upper_x, prim->upper_y, prim->upper_z, 0.0f));

      leftBox.boxMax [splitAxis] = splitPos;
      rightBox.boxMin[splitAxis] = splitPos;

      // Intersect with original bounds.
      //
      leftBox.intersect(nodeBox);
      rightBox.intersect(nodeBox);

      lprim->lower_x = leftBox.boxMin.x;
      lprim->lower_y = leftBox.boxMin.y;
      lprim->lower_z = leftBox.boxMin.z;
      lprim->upper_x = leftBox.boxMax.x;
      lprim->upper_y = leftBox.boxMax.y;
      lprim->upper_z = leftBox.boxMax.z;

      rprim->lower_x = rightBox.boxMin.x;
      rprim->lower_y = rightBox.boxMin.y;
      rprim->lower_z = rightBox.boxMin.z;
      rprim->upper_x = rightBox.boxMax.x;
      rprim->upper_y = rightBox.boxMax.y;
      rprim->upper_z = rightBox.boxMax.z;
    }
  }

  struct Node
  {
    Node(){}
    virtual ~Node(){}
    virtual float sah() const = 0;
    virtual bool  isLeaf() const = 0;
    uint32_t primsInNode;
  };

  static inline Box4f merge(const Box4f& a, const Box4f& b) { return Box4f(min(a.boxMin, b.boxMin), max(a.boxMax, b.boxMax)); }

  struct InnerNode : public Node
  {
    Box4f bounds[2];
    Node* children[2];
  
    InnerNode() {
      bounds  [0] = bounds[1]   = Box4f();
      children[0] = children[1] = nullptr;
      primsInNode = 0;
    }

    float sah() const override { 
      return 1.0f + (bounds[0].surfaceArea()*children[0]->sah() + bounds[1].surfaceArea()*children[1]->sah())/merge(bounds[0],bounds[1]).surfaceArea();
    }

    bool  isLeaf() const override { return false; }

    static void* create (RTCThreadLocalAllocator alloc, unsigned int numChildren, void* userPtr)
    {
      assert(numChildren == 2);
      void* ptr = rtcThreadLocalAlloc(alloc,sizeof(InnerNode),16);
      return (void*) new (ptr) InnerNode;
    }

    static void  setChildren (void* nodePtr, void** childPtr, unsigned int numChildren, void* userPtr)
    {
      assert(numChildren == 2);
      for (size_t i=0; i<2; i++)
        ((InnerNode*)nodePtr)->children[i] = (Node*) childPtr[i];
    }

    static void  setBounds (void* nodePtr, const RTCBounds** bounds, unsigned int numChildren, void* userPtr)
    {
      assert(numChildren == 2);
      for (size_t i=0; i<2; i++)
        ((InnerNode*)nodePtr)->bounds[i] = *(const Box4f*) bounds[i];
    }
  };

  struct LeafNode : public Node
  {
    unsigned id;
    Box4f bounds;

    LeafNode (unsigned a_id, const Box4f& a_bounds) : id(a_id), bounds(a_bounds) {}

    float sah() const override {
      return 1.0f;
    }

    bool  isLeaf() const override { return true; }

    static void* create (RTCThreadLocalAllocator alloc, const RTCBuildPrimitive* prims, size_t numPrims, void* userPtr)
    {
      assert(numPrims == 1);
      void* ptr = rtcThreadLocalAlloc(alloc,sizeof(LeafNode),16);
      return (void*) new (ptr) LeafNode(prims->primID,*(Box4f*)prims);
    }
  };

  static uint32_t EvaluatePrimsInNode(Node* node)
  {
    if(node == nullptr)
      return 0;

    if(node->isLeaf())
    {
      node->primsInNode = 1;
      return 1;
    }
    
    InnerNode* pInnderNode = (InnerNode*)node;
    node->primsInNode      = EvaluatePrimsInNode(pInnderNode->children[0]) + EvaluatePrimsInNode(pInnderNode->children[1]); 
    return node->primsInNode;
  }

  static std::vector<uint32_t> g_prims;

  static void GetPimListToGPrims(Node* node)
  {
    if(node == nullptr)
      return;
    
    if(node->isLeaf())
    {
      LeafNode* leaf = (LeafNode*)node;
      g_prims.push_back(leaf->id);
    }  
    else
    {
      InnerNode* pInnderNode = (InnerNode*)node;
      GetPimListToGPrims(pInnderNode->children[0]);
      GetPimListToGPrims(pInnderNode->children[1]);
    }
  }  

  struct NextNode
  {
    uint32_t       leftOffset = uint32_t(-1);
    cbvh::Interval interval   = cbvh::Interval(0,0);
  };

  static NextNode ConvertTreeRec(Node* node, const Triangle4f* a_tris, cbvh::BVHTree& a_res, uint32_t a_parentEscapeIndex)
  {
    if(node == nullptr)
      return NextNode();
    
    NextNode nodeData;
    if(node->isLeaf())
    {
      LeafNode* leaf = (LeafNode*)node;
      const uint32_t oldTris = a_res.indicesReordered.size();
      a_res.indicesReordered.push_back(leaf->id);
      
      cbvh::BVHNode currNode;
      currNode.boxMin = to_float3(leaf->bounds.boxMin);
      currNode.boxMax = to_float3(leaf->bounds.boxMax);
      currNode.leftOffset  = cbvh::LEAF_NORMAL;
      currNode.escapeIndex = cbvh::LEAF_NORMAL;
      nodeData.interval    = cbvh::Interval(oldTris, 1);
      a_res.intervals.push_back(nodeData.interval);
      a_res.nodes.push_back(currNode);
      return nodeData;      
    }

    if(node->primsInNode <= g_recommendedPrimsInLeaf) // #TODO: estimate relative bounding box overlap, probably we can do better than 4 triangle in leaf
    {
      InnerNode* pInnderNode = (InnerNode*)node;
      const float childrenOverlap = cbvh_internal::RelativeBoxOverlap(pInnderNode->bounds[0], pInnderNode->bounds[1]);
      //std::cout << "childrenOverlap = " << childrenOverlap << std::endl;
      if(childrenOverlap > g_overSplitThreshold)  // 0.25--0.5f
      {
        g_prims.resize(0);
        GetPimListToGPrims(node);
        
        // put triangles in side reordered buffer
        //
        const uint32_t oldTris = a_res.indicesReordered.size();
        a_res.indicesReordered.resize(a_res.indicesReordered.size() + g_prims.size());
        for(size_t triId = 0; triId < g_prims.size(); triId++)
          a_res.indicesReordered[oldTris + triId] = g_prims[triId];
        
        cbvh::BVHNode currNode;
        {
          currNode.boxMin = to_float3(min(pInnderNode->bounds[0].boxMin, pInnderNode->bounds[1].boxMin));
          currNode.boxMax = to_float3(max(pInnderNode->bounds[0].boxMax, pInnderNode->bounds[1].boxMax));
          currNode.leftOffset  = cbvh::LEAF_NORMAL;
          currNode.escapeIndex = cbvh::LEAF_NORMAL;
        }
        nodeData.interval = cbvh::Interval(oldTris, uint32_t(g_prims.size()));
        a_res.intervals.push_back(nodeData.interval);
        a_res.nodes.push_back(currNode);
        return nodeData;
      }
    }
    
    
    InnerNode* pInnderNode = (InnerNode*)node;
    cbvh::BVHNode leftNode, rightNode;
    {
      leftNode.boxMin  = to_float3(pInnderNode->bounds[0].boxMin);
      leftNode.boxMax  = to_float3(pInnderNode->bounds[0].boxMax);
      rightNode.boxMin = to_float3(pInnderNode->bounds[1].boxMin);
      rightNode.boxMax = to_float3(pInnderNode->bounds[1].boxMax);
      
      leftNode.leftOffset   = cbvh::LEAF_NORMAL;
      leftNode.escapeIndex  = cbvh::LEAF_NORMAL;
      rightNode.leftOffset  = cbvh::LEAF_NORMAL;
      rightNode.escapeIndex = cbvh::LEAF_NORMAL;
    }
      
    nodeData.leftOffset = uint32_t(a_res.nodes.size());
    a_res.intervals.push_back(cbvh::Interval(0,0));
    a_res.intervals.push_back(cbvh::Interval(0,0));
    a_res.nodes.push_back(leftNode);
    a_res.nodes.push_back(rightNode);
    
    a_res.nodes[nodeData.leftOffset + 0].escapeIndex = nodeData.leftOffset + 1;
    a_res.nodes[nodeData.leftOffset + 1].escapeIndex = a_parentEscapeIndex;

    NextNode dataLeft  = ConvertTreeRec(pInnderNode->children[0], a_tris, a_res, a_res.nodes[nodeData.leftOffset + 0].escapeIndex);
    NextNode dataRight = ConvertTreeRec(pInnderNode->children[1], a_tris, a_res, a_res.nodes[nodeData.leftOffset + 1].escapeIndex);

    a_res.nodes[nodeData.leftOffset + 0].leftOffset = dataLeft.leftOffset;  // fix leftOffset after we done with children
    a_res.nodes[nodeData.leftOffset + 1].leftOffset = dataRight.leftOffset; // fix leftOffset after we done with children
    a_res.intervals[nodeData.leftOffset + 0] = dataLeft.interval;           // same for intervals
    a_res.intervals[nodeData.leftOffset + 1] = dataRight.interval;          // same for intervals

    assert(dataLeft.interval.start + dataLeft.interval.count == dataRight.interval.start);
    nodeData.interval = cbvh::Interval(dataLeft.interval.start, dataLeft.interval.count + dataRight.interval.count);
    
    return nodeData;
  }

  static void ConvertTree(Node* node, const Triangle4f* a_tris, cbvh::BVHTree& a_res)
  {
    uint32_t totalPrims = EvaluatePrimsInNode(node);
    
    a_res.nodes.resize(0);            a_res.nodes.reserve(totalPrims/2);
    a_res.intervals.resize(0);        a_res.intervals.reserve(totalPrims/2);
    a_res.indicesReordered.resize(0); a_res.indicesReordered.reserve(totalPrims);
    
    g_prims.reserve(32);
    
    // (1) put root node first
    //
    assert(!node->isLeaf());
    InnerNode* pInnderNode = (InnerNode*)node;
    cbvh::BVHNode root;
    {
      root.boxMin = to_float3(min(pInnderNode->bounds[0].boxMin, pInnderNode->bounds[1].boxMin));
      root.boxMax = to_float3(max(pInnderNode->bounds[0].boxMax, pInnderNode->bounds[1].boxMax));  
      root.leftOffset  = cbvh::LEAF_NORMAL;
      root.escapeIndex = cbvh::ESCAPE_ROOT;
    }

    a_res.intervals.push_back(cbvh::Interval(0,pInnderNode->primsInNode));
    a_res.nodes.push_back(root);

    // (2) then convert tree format to out BVH2_Static format
    //
    NextNode rootData         = ConvertTreeRec(node, a_tris, a_res, 0xFFFFFFFE);
    a_res.nodes[0].leftOffset = rootData.leftOffset; 
    a_res.intervals[0]        = rootData.interval;
  }

  void build2(RTCBuildQuality quality, std::vector<RTCBuildPrimitive>& prims_i, const std::vector<Triangle4f>& a_tris, size_t extraSpace, cbvh::BVHTree& a_res)
  {
    auto a_device = rtcNewDevice("");
    rtcSetDeviceMemoryMonitorFunction(a_device,memoryMonitor,nullptr);

    RTCBVH bvh = rtcNewBVH(a_device);

    auto& prims = prims_i;
    prims.reserve(prims_i.size()+extraSpace);
    prims.resize(prims_i.size());

    UserData user;
    user.tris = a_tris.data();

    /* settings for BVH build */
    RTCBuildArguments arguments = rtcDefaultBuildArguments();
    arguments.byteSize = sizeof(arguments);
    arguments.buildFlags = RTC_BUILD_FLAG_NONE; //RTC_BUILD_FLAG_DYNAMIC; // ????????????????????
    arguments.buildQuality = quality;
    arguments.maxBranchingFactor = 2;
    arguments.maxDepth = 1024;
    arguments.sahBlockSize = 1;
    arguments.minLeafSize = 1;
    arguments.maxLeafSize = 1;
    arguments.traversalCost = 1.0f;
    arguments.intersectionCost = 1.0f;
    arguments.bvh = bvh;
    arguments.primitives     = prims.data();
    arguments.primitiveCount = prims.size();
    arguments.primitiveArrayCapacity = prims.capacity();
    arguments.createNode = InnerNode::create;
    arguments.setNodeChildren = InnerNode::setChildren;
    arguments.setNodeBounds = InnerNode::setBounds;
    arguments.createLeaf = LeafNode::create;
    arguments.splitPrimitive = splitPrimitive;
    arguments.buildProgress = buildProgress;
    arguments.userPtr = &user;
    
    auto before = std::chrono::high_resolution_clock::now();
    Node* root = (Node*) rtcBuildBVH(&arguments);
    float time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - before).count()/1000.f;
    g_buildTime += double(time);

    ConvertTree(root, a_tris.data(), a_res);

    rtcReleaseBVH(bvh);
    rtcReleaseDevice(a_device); 
  }
}

namespace cbvh_embree
{
  cbvh::BVHTree BuildBVH(const float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, cbvh::BVHPresets a_settings)
  {
    RTCBuildQuality quality = RTC_BUILD_QUALITY_HIGH;    
    if(a_settings.btype == cbvh::BVH_BUILDER_TYPE::BVH_CONSTRUCT_EMBREE_FAST)
      quality = RTC_BUILD_QUALITY_LOW;

    std::vector<RTCBuildPrimitive> primBoxData(a_indexNum/3);
    std::vector<Triangle4f>        primTriData(a_indexNum/3);

    for(size_t i=0;i<a_indexNum/3;i++) 
    {
      const uint iA = a_indices[i*3+0];
      const uint iB = a_indices[i*3+1];
      const uint iC = a_indices[i*3+2];

      Triangle4f tri;
      tri.A = a_vertices[iA];
      tri.B = a_vertices[iB];
      tri.C = a_vertices[iC];
      tri.A.w = LiteMath::as_float(int(i));
      primTriData[i] = tri;

      Box4f box;
      box.boxMin = min(tri.A, min(tri.B, tri.C));
      box.boxMax = max(tri.A, max(tri.B, tri.C));
      
      RTCBuildPrimitive prim;
      prim.lower_x = box.boxMin.x;
      prim.lower_y = box.boxMin.y;
      prim.lower_z = box.boxMin.z;
      prim.geomID = 0;
      prim.upper_x = box.boxMax.x;
      prim.upper_y = box.boxMax.y;
      prim.upper_z = box.boxMax.z;
      prim.primID = (unsigned) i;
      primBoxData[i] = prim;
    }

    embree::g_recommendedPrimsInLeaf = a_settings.primsInLeaf;
    embree::g_overSplitThreshold     = 0.5f;
    if(a_settings.primsInLeaf >= 4)
      embree::g_overSplitThreshold = 0.25f;

    cbvh::BVHTree res;
    embree::build2(quality, primBoxData, primTriData, (a_indexNum/6), res);  
    res.format = cbvh::CBVH_FORMATS::FMT_BVH2Node32_Interval32_Static;
    return res;
  }
};
