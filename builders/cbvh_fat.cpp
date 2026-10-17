#include <fstream>
#include <chrono>

#include "cbvh.h"
#include "cbvh_core.h"

using cbvh2::BVHNode;
using LiteMath::float4;

//#include "nanort/nanort.h"
#include "FatBVH.h"
//#include "NonunifNodeStorage.h"

int cbvh2::CalcInputGroupSize(const LayoutPresets& a_layout)
{
  int res = 0;
  if(a_layout.grSzXX <= 0)
  {
    if(a_layout.layout == cbvh2::LAYOUT_COLBVH_YM06)
      res = std::pow(2, a_layout.grSzId+3)-1; // std::log2(int(cluster_size+1))-3;
    else
      res = std::pow(2, a_layout.grSzId+1); // std::log2(int(leavesNumberInTreelet))-1;
  }
  else
    res = a_layout.grSzXX;
  return res;
}

int cbvh2::CalcActualGroupSize(const LayoutPresets& a_layout)
{
  int res = 0;
  if(a_layout.grSzXX <= 0)
  {
    if(a_layout.layout == cbvh2::LAYOUT_COLBVH_YM06)
      res = std::pow(2, a_layout.grSzId+3)-1; 
    else
      res = std::pow(2, a_layout.grSzId+2)+1;  
  }
  else
  {
    if(a_layout.layout == cbvh2::LAYOUT_COLBVH_YM06)
      res = a_layout.grSzXX;
    else
      res = 2*a_layout.grSzXX+1;
  }
  return res;
}

cbvh2::LayoutPresets cbvh2::LayoutPresetsFromString(const char* a_str)
{
  LayoutPresets layout;
  std::string m_layoutName(a_str);

  if(m_layoutName == "DepthFirst" || m_layoutName == "DFS" || m_layoutName == "DFL")
  {
    layout.layout = cbvh2::LAYOUT_DFS;
    layout.grSzId = 0;
    layout.grSzXX = 0;
  }
  else if(m_layoutName == "OrderedDepthFirst" || m_layoutName == "ODFS" || m_layoutName == "ODFL")
  {
    layout.layout = cbvh2::LAYOUT_ODFS;
    layout.grSzId = 0;
    layout.grSzXX = 0;
  }
  else if(m_layoutName == "BreadthFirst" || m_layoutName == "BFS" || m_layoutName == "BFL")
  {
    layout.layout = cbvh2::LAYOUT_BFS;
    layout.grSzId = 0;
    layout.grSzXX = 0;
  }
  else if (m_layoutName.find("Treelet") != std::string::npos || m_layoutName.find("TRB") != std::string::npos)
  {
    size_t leavesNumberInTreelet = 4;
    std::string str = m_layoutName.substr(m_layoutName.size() - 2);

    if (str.find_first_not_of("0123456789") == std::string::npos) 
      leavesNumberInTreelet = stoi(str);
    else
    {
      str = m_layoutName.substr(m_layoutName.size() - 1);
      if (str.find_first_not_of("0123456789") == std::string::npos)
        leavesNumberInTreelet = stoi(str);
    }

    layout.grSzId = std::log2(int(leavesNumberInTreelet))-1;
    layout.grSzXX = int(leavesNumberInTreelet); // actual size is 2*int(leavesNumberInTreelet)-1;

    if (m_layoutName.find("SuperSuperTreelet") != std::string::npos)
      layout.layout = cbvh2::LAYOUT_COLBVH_TRB3;
    else if (m_layoutName.find("SuperTreelet") != std::string::npos)
      layout.layout = cbvh2::LAYOUT_COLBVH_TRB2;
    else
      layout.layout = cbvh2::LAYOUT_COLBVH_TRB1;
   
    const bool a_aligned = m_layoutName.find("Aligned") != std::string::npos;
    const bool a_merge   = m_layoutName.find("Merged") != std::string::npos;

    if (m_layoutName.find("SuperTreelet") != std::string::npos && a_aligned && a_merge)
      layout.layout = cbvh2::LAYOUT_CALBVH;
  }
  else if (m_layoutName.find("Clusterized") != std::string::npos || m_layoutName.find("TRB") != std::string::npos || m_layoutName == "opt")
  {
    size_t cluster_size = 31; // default of MAX size of the lowest-level (smallest) cluster
    std::string str = m_layoutName.substr(m_layoutName.size() - 2);
    if (str.find_first_not_of("0123456789") == std::string::npos)
      cluster_size = stoi(str);
    else
    {
      str = m_layoutName.substr(m_layoutName.size() - 1);
      if (str.find_first_not_of("0123456789") == std::string::npos)
        cluster_size = stoi(str);
    }
    
    layout.layout = cbvh2::LAYOUT_COLBVH_YM06;
    layout.grSzId = std::log2(int(cluster_size+1))-3;
    layout.grSzXX = int(cluster_size);          // actual size is cluster_size+1;
  }

  return layout;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<BVHNode>  g_dummy1;
std::vector<Interval> g_dummy2;

struct TreeConverter
{
  TreeConverter(const std::vector<BVHNode>&  a_input,
                const std::vector<Interval>& a_intervals,
                std::vector<BVHNodeFat>&     a_out) : m_input(a_input),
                                                      m_intervals(a_intervals),
                                                      m_out(a_out) {}
                                                      //m_tree(nullptr) {}

  //TreeConverter(const nanort::BVHAccel<float>& a_tree,
  //              std::vector<BVHNodeFat>&       a_out) : m_input(g_dummy1),
  //                                                      m_intervals(g_dummy2),
  //                                                      m_out(a_out), 
  //                                                      m_tree(&a_tree) {}
    
  const std::vector<BVHNode>&    m_input;
  const std::vector<Interval>&   m_intervals;

  std::vector<BVHNodeFat>&       m_out;
  //const nanort::BVHAccel<float>* m_tree;

  uint32_t ProcessBVHNode(uint32_t currNodeId);
  uint32_t ProcessNRTNode(uint32_t currNodeId);
};

uint32_t TreeConverter::ProcessBVHNode(uint32_t currNodeId)
{
  BVHNode currNode = m_input[currNodeId];

  if(currNode.leftOffset & LEAF_BIT) // should not happen in general, except for single node tree
  {
    BVHNodeFat fatNode;
    fatNode.lmin_xyz_rmax_x = to_float4(currNode.boxMin, currNode.boxMax.x);
    fatNode.lmax_xyz_rmax_y = to_float4(currNode.boxMax, currNode.boxMax.y);
    fatNode.rmin_xyz_rmax_z = to_float4(currNode.boxMin, currNode.boxMax.z);
    fatNode.offs_left       = PackOffsetAndSize(m_intervals[currNodeId].start, m_intervals[currNodeId].count); 
    fatNode.offs_right      = PackOffsetAndSize(0, 0); 
    m_out.push_back(fatNode);
    return uint32_t(m_out.size()-1);
  }
  
  const uint32_t leftOffset  = currNode.leftOffset + 0;
  const uint32_t rightOffset = currNode.leftOffset + 1;

  const BVHNode leftNode     = m_input[leftOffset];
  const BVHNode rightNode    = m_input[rightOffset];
  
  assert(m_out.capacity() > m_out.size()+1);
  m_out.push_back(BVHNodeFat());
  uint32_t currNodeIndex  = uint32_t(m_out.size()-1);
  BVHNodeFat& fatNode     = m_out[currNodeIndex];
  fatNode.lmin_xyz_rmax_x = to_float4(leftNode.boxMin,  rightNode.boxMax.x);
  fatNode.lmax_xyz_rmax_y = to_float4(leftNode.boxMax,  rightNode.boxMax.y);
  fatNode.rmin_xyz_rmax_z = to_float4(rightNode.boxMin, rightNode.boxMax.z);

  if(leftNode.leftOffset & LEAF_BIT)
    fatNode.offs_left = PackOffsetAndSize(m_intervals[leftOffset].start, m_intervals[leftOffset].count); 
  else
    fatNode.offs_left = ProcessBVHNode(leftOffset);

  if(rightNode.leftOffset & LEAF_BIT)
    fatNode.offs_right = PackOffsetAndSize(m_intervals[rightOffset].start, m_intervals[rightOffset].count);
  else
    fatNode.offs_right = ProcessBVHNode(rightOffset);
  
  return currNodeIndex;
}

//uint32_t TreeConverter::ProcessNRTNode(uint32_t currNodeId)
//{
//  nanort::BVHNode currNode = m_tree->GetNodes()[currNodeId];
//
//  if(currNode.flag == 1) // should not happen in general, except for single node tree
//  {
//    BVHNodeFat fatNode;
//    fatNode.lmin_xyz_rmax_x = to_float4(float3(currNode.bmin), currNode.bmax[0]);
//    fatNode.lmax_xyz_rmax_y = to_float4(float3(currNode.bmax), currNode.bmax[1]);
//    fatNode.rmin_xyz_rmax_z = to_float4(float3(currNode.bmin), currNode.bmax[2]);
//    fatNode.offs_left       = PackOffsetAndSize(currNode.data[1], currNode.data[0]); 
//    fatNode.offs_right      = PackOffsetAndSize(0, 0); 
//    m_out.push_back(fatNode);
//    return uint32_t(m_out.size()-1);
//  }
//  
//  const uint32_t leftOffset  = currNode.data[0];
//  const uint32_t rightOffset = currNode.data[1];
//
//  const nanort::BVHNode leftNode     = m_tree->GetNodes()[leftOffset];
//  const nanort::BVHNode rightNode    = m_tree->GetNodes()[rightOffset];
//  
//  m_out.push_back(BVHNodeFat());
//  uint32_t currNodeIndex  = uint32_t(m_out.size()-1);
//  BVHNodeFat& fatNode     = m_out[currNodeIndex];
//  fatNode.lmin_xyz_rmax_x = to_float4(float3(leftNode.bmin),  rightNode.bmax[0]);
//  fatNode.lmax_xyz_rmax_y = to_float4(float3(leftNode.bmax),  rightNode.bmax[1]);
//  fatNode.rmin_xyz_rmax_z = to_float4(float3(rightNode.bmin), rightNode.bmax[2]);
//
//  if(leftNode.flag == 1)
//    fatNode.offs_left = PackOffsetAndSize(leftNode.data[1], leftNode.data[0]); 
//  else
//    fatNode.offs_left = ProcessNRTNode(leftOffset);
//
//  if(rightNode.flag == 1)
//    fatNode.offs_right = PackOffsetAndSize(rightNode.data[1], rightNode.data[0]);
//  else
//    fatNode.offs_right = ProcessNRTNode(rightOffset);
//  
//  return currNodeIndex;
//}

std::vector<BVHNodeFat> CreateFatTreeArray(const std::vector<BVHNode>& a_input, const std::vector<Interval>& a_intervals)
{
  std::vector<BVHNodeFat> result;
  result.reserve(a_input.size()); // this is important, array should not be reallocated!
  result.resize(0);                       

  TreeConverter tc(a_input, a_intervals, result);
  tc.ProcessBVHNode(0);
  
  result.shrink_to_fit();         
  return result;
}


struct TreeConverter2
{
  TreeConverter2(const std::vector<BVHNode>&  a_input,
                 std::vector<BVHNodeFat>&     a_out) : m_input(a_input),
                                                       m_out(a_out) {}
                                                       //m_tree(nullptr) {}

  //TreeConverter2(const nanort::BVHAccel<float>& a_tree,
  //               std::vector<BVHNodeFat>&       a_out) : m_input(g_dummy1),
  //                                                       m_out(a_out), 
  //                                                       m_tree(&a_tree) {}
    
  const std::vector<BVHNode>&    m_input;
  std::vector<BVHNodeFat>&       m_out;
  //const nanort::BVHAccel<float>* m_tree;

  uint32_t ProcessBVHNode(uint32_t currNodeId);
  uint32_t ProcessNRTNode(uint32_t currNodeId);
};

uint32_t TreeConverter2::ProcessBVHNode(uint32_t leftOffset)
{
  assert(leftOffset + 1 < m_input.size());
  assert(m_out.capacity() > m_out.size()+1);

  const BVHNode leftNode  = m_input[leftOffset + 0];
  const BVHNode rightNode = m_input[leftOffset + 1];

  m_out.push_back(BVHNodeFat());
  uint32_t currNodeIndex  = uint32_t(m_out.size()-1);
  BVHNodeFat& fatNode     = m_out[currNodeIndex];
  fatNode.lmin_xyz_rmax_x = to_float4(leftNode.boxMin,  rightNode.boxMax.x);
  fatNode.lmax_xyz_rmax_y = to_float4(leftNode.boxMax,  rightNode.boxMax.y);
  fatNode.rmin_xyz_rmax_z = to_float4(rightNode.boxMin, rightNode.boxMax.z);

  if(leftNode.leftOffset & LEAF_BIT)
    fatNode.offs_left = leftNode.leftOffset;
  else
    fatNode.offs_left = ProcessBVHNode(leftNode.leftOffset);

  if(rightNode.leftOffset & LEAF_BIT)
    fatNode.offs_right = rightNode.leftOffset;
  else
    fatNode.offs_right = ProcessBVHNode(rightNode.leftOffset);
  
  return currNodeIndex;
}


std::vector<BVHNodeFat> CreateFatTreeArray(const std::vector<BVHNode>& a_input)
{
  std::vector<BVHNodeFat> result;
  result.reserve(a_input.size()); // this is important, array should not be reallocated!
  result.resize(0);

  //for(size_t i=0;i<a_input.size();i++)
  //  std::cout << i << ": (" << a_input[i].leftOffset << "," << a_input[i].escapeIndex << ")" << std::endl;

  TreeConverter2 tc(a_input, result);
  tc.ProcessBVHNode(0);
  
  result.shrink_to_fit();         
  return result;
}

//std::vector<BVHNodeFat> CreateFatTreeArray(const nanort::BVHAccel<float>& a_tree, std::vector<uint32_t>& a_indicesReordered)
//{
//  const size_t n_trgs = a_tree.GetIndices().size();
//  a_indicesReordered.resize(n_trgs);
//  for (size_t i = 0; i < n_trgs; i++)
//    a_indicesReordered[i] = a_tree.GetIndices()[i];
//
//  std::vector<BVHNodeFat> result;
//  result.reserve(a_tree.GetNodes().size()); // this is important, array should not be reallocated!
//  result.resize(0);
//
//  TreeConverter tc(a_tree, result);
//  tc.ProcessNRTNode(0);
//  
//  result.shrink_to_fit();         
//  return result;
//}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern double g_buildTime;
extern uint64_t g_buildTris;

cbvh2::BVHTreeFat cbvh2::BuildBVHFat(const float* a_vpos3f,     size_t a_vertNum, size_t a_vByteStride, 
                                     const uint32_t* a_indices, size_t a_indexNum, BuilderPresets a_presets, LayoutPresets a_layout,
                                     std::vector<int>* pTreelet_root, std::vector<int>* pSuper_treelet_root, std::vector<int>* pGrStart)
{
  std::vector<BVHNodeFat> bvhFat;
  std::vector<uint32_t>   objIndicesReordered;
  
  // (1) build
  //
  //if (a_presets.quality == cbvh2::BVH_CONSTRUCT_NANORT)
  //{
  //  const size_t vStride = a_vByteStride / 4;
  //  const size_t n_trgs  = a_indexNum / 3;
  //  g_buildTris += n_trgs;
//
  //  std::vector<float> verts;
  //  verts.resize(a_vertNum * 3);
  //  for (size_t y = 0; y < a_vertNum; y++)
  //  {
  //    verts[y * 3 + 0] = a_vpos3f[y * vStride + 0];
  //    verts[y * 3 + 1] = a_vpos3f[y * vStride + 1];
  //    verts[y * 3 + 2] = a_vpos3f[y * vStride + 2];
  //  }
//
  //  nanort::BVHAccel<float> accel;
  //  nanort::BVHBuildOptions<float> build_options;
  //  build_options.cache_bbox = false;
  //  build_options.min_leaf_primitives = a_presets.primsInLeaf;
  //  nanort::TriangleMesh<float>    triangle_mesh(verts.data(), a_indices, sizeof(float) * 3);
  //  nanort::TriangleSAHPred<float> triangle_pred(verts.data(), a_indices, sizeof(float) * 3);
  //  
  //  auto before = std::chrono::high_resolution_clock::now();
  //  bool res = accel.Build(n_trgs, triangle_mesh, triangle_pred, build_options);
  //  if(!res)
  //  {
  //    std::cout << "nanort::Build is failed!" << std::endl; 
  //    exit(0);
  //  }
  //  float time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - before).count()/1000.f;
  //  g_buildTime += double(time); //
  //  
  //  bvhFat = CreateFatTreeArray(accel, objIndicesReordered);
  //}
  //else
  {
    a_presets.format    = cbvh2::BVH2_LEFT_OFFSET;
    auto bvhData        = cbvh2::BuildBVH(a_vpos3f, a_vertNum, a_vByteStride, a_indices, a_indexNum, a_presets);
    bvhFat              = CreateFatTreeArray(bvhData.nodes);
    objIndicesReordered = bvhData.indices;
  }
  
  // (2) Reorder according to required layout
  //
  std::vector<int> treelet_root, super_treelet_root, grStart;
  std::vector<uint2> depthRanges; 

  if(a_layout.layout == cbvh2::LAYOUT_DFS)
    FatBVH::ReorderDFL(bvhFat);
  else if(a_layout.layout == cbvh2::LAYOUT_ODFS)
    FatBVH::ReorderODFL(bvhFat);
  else if(a_layout.layout == cbvh2::LAYOUT_BFS)
  {
    FatBVH::ReorderBFL(bvhFat);
    depthRanges = FatBVH::ComputeDepthRanges(bvhFat);
  }
  else if (a_layout.layout >= LAYOUT_COLBVH_TRB1 &&  a_layout.layout <= LAYOUT_CALBVH)
  {
    size_t leavesNumberInTreelet = cbvh2::CalcInputGroupSize(a_layout);

    FatBVH::TreeletizLevel level;
    if (a_layout.layout == cbvh2::LAYOUT_COLBVH_TRB3)
      level = FatBVH::TreeletizLevel::SUPERSUPERTREELETS;
    else if (a_layout.layout == cbvh2::LAYOUT_COLBVH_TRB2 || a_layout.layout == cbvh2::LAYOUT_CALBVH)
      level = FatBVH::TreeletizLevel::SUPERTREELETS;
    else
      level = FatBVH::TreeletizLevel::TREELETS; 

    if (a_layout.layout != cbvh2::LAYOUT_CALBVH)
      FatBVH::ReorderTRB(bvhFat, treelet_root, super_treelet_root, leavesNumberInTreelet, false, level);
    else
      FatBVH::ReorderTRBNew(bvhFat, treelet_root, super_treelet_root, grStart, leavesNumberInTreelet, true, true);
  }
  else if (a_layout.layout == cbvh2::LAYOUT_COLBVH_YM06)
  {
    size_t cluster_size = cbvh2::CalcInputGroupSize(a_layout); // default of MAX size of the lowest-level (smallest) cluster
    FatBVH::ReorderByClusters(bvhFat, treelet_root, 0, cluster_size, false);
  }
  
  // output these data if we need it
  //
  if(pTreelet_root != nullptr)
    (*pTreelet_root) = treelet_root;
  if(pSuper_treelet_root != nullptr)
    (*pSuper_treelet_root) = super_treelet_root;
  if(pGrStart != nullptr)
    (*pGrStart) = grStart;

  //if(bvhFat.size() > 800) {
  //  auto& treelets = (m_layoutName.find("Aligned") != std::string::npos) ? grStart : treelet_root; 
  //  FatBVH::PrintForGraphViz(bvhFat, treelets, super_treelet_root, "z_nodes.txt"); 
  //}

  return BVHTreeFat(bvhFat, objIndicesReordered, depthRanges);
}

#ifdef HALFFLOAT 
static inline half_float::half roundToNeg(float a_val) { return half_float::nextafter(half_float::half(a_val), -std::numeric_limits<half_float::half>::max()); }
static inline half_float::half roundToPos(float a_val) { return half_float::nextafter(half_float::half(a_val), +std::numeric_limits<half_float::half>::max()); }

cbvh2::BVHTreeFat16 cbvh2::BuildBVHFat16(const float* a_vpos3f, size_t a_vertNum, size_t a_vByteStride, const uint32_t* a_indices, size_t a_indexNum, 
                                         BuilderPresets a_presets, LayoutPresets a_layout)
{
  BVHTreeFat originalBvh = BuildBVHFat(a_vpos3f, a_vertNum, a_vByteStride, a_indices, a_indexNum, a_presets, a_layout);
  
  std::vector<BVHNodeFat16>  nodes16(originalBvh.nodes.size());
  for(size_t i=0;i<nodes16.size();i++) 
  {
    auto node     = originalBvh.nodes[i];
    auto leftBox  = GetChildBoxLeft(node);
    auto rightBox = GetChildBoxRight(node);

    BVHNodeFat16 node16;
    {
      node16.lmin_xyz_rmax_x = half4(roundToNeg(leftBox.boxMin.x),  roundToNeg(leftBox.boxMin.y),  roundToNeg(leftBox.boxMin.z),  roundToPos(rightBox.boxMax.x));
      node16.lmax_xyz_rmax_y = half4(roundToPos(leftBox.boxMax.x),  roundToPos(leftBox.boxMax.y),  roundToPos(leftBox.boxMax.z),  roundToPos(rightBox.boxMax.y));
      node16.rmin_xyz_rmax_z = half4(roundToNeg(rightBox.boxMin.x), roundToNeg(rightBox.boxMin.y), roundToNeg(rightBox.boxMin.z), roundToPos(rightBox.boxMax.z));
      node16.offs_left       = node.offs_left;
      node16.offs_right      = node.offs_right;
    }

    nodes16[i] = node16;
  }

  return BVHTreeFat16(nodes16, originalBvh.indices);
}

static inline uint4 PackNode(const BVHNode& node)
{
  LiteMath::half XY[2];
  LiteMath::half ZX[2];
  LiteMath::half YZ[2];

  XY[0] = roundToNeg(node.boxMin.x);
  XY[1] = roundToNeg(node.boxMin.y);
  ZX[0] = roundToNeg(node.boxMin.z);

  ZX[1] = roundToPos(node.boxMax.x);
  YZ[0] = roundToPos(node.boxMax.y);
  YZ[1] = roundToPos(node.boxMax.z);
  
  uint4 res;
  memcpy(&res.x, XY, sizeof(uint32_t));
  memcpy(&res.y, ZX, sizeof(uint32_t));
  memcpy(&res.z, YZ, sizeof(uint32_t));
  res.w = node.leftOffset;
  return res;
}

cbvh2::BVHTreeFatCompressed cbvh2::BuildBVHHalf(const float* a_vpos3f,     size_t a_vertNum, size_t a_vByteStride, 
                                                const uint32_t* a_indices, size_t a_indexNum, BuilderPresets a_presets)
{

  BVHTreeCommon originalBvh = cbvh2::BuildBVH(a_vpos3f, a_vertNum, a_vByteStride, a_indices, a_indexNum, a_presets);
  std::vector<uint4> nodesHalf(originalBvh.nodes.size());
  for(size_t i=0;i<nodesHalf.size();i++)
    nodesHalf[i] = PackNode(originalBvh.nodes[i]);
  return BVHTreeFatCompressed(nodesHalf, originalBvh.indices, 0);
}
#endif


/*cbvh2::BVHTreeFatCompressed cbvh2::BuildBVHFatCompressed(const float* a_vpos3f,     size_t a_vertNum, size_t a_vByteStride, 
                                                         const uint32_t* a_indices, size_t a_indexNum, BuilderPresets a_presets, LayoutPresets a_layout)
{
  a_layout.layout = cbvh2::LAYOUT_CALBVH;
  std::vector<int> treelet_root, super_treelet_root, grStart;
  auto fatTree = BuildBVHFat(a_vpos3f, a_vertNum, a_vByteStride, a_indices, a_indexNum, a_presets, a_layout,
                             &treelet_root, &super_treelet_root, &grStart);

  NonunifNodeStorage nns;

  std::vector<uint32_t> a_bvhOffsets;
  a_bvhOffsets.emplace_back(0); // for the new SINGLE tree

  int res = -1;
  if (grStart.empty() && !treelet_root.empty())
    res = nns.Create(fatTree.nodes, a_bvhOffsets, treelet_root, treelet_root); // groups = treelets
  else
    res = nns.Create(fatTree.nodes, a_bvhOffsets, treelet_root, grStart);

  if (res != 0)
  {
    assert(false);
    return BVHTreeFatCompressed();
  }

  while(nns.data.size() % (2*a_layout.grSzXX) != 0)
    nns.data.push_back(uint4(0,0,0,0));

  return BVHTreeFatCompressed(nns.data, fatTree.indices, nns.m_log2_group_size);
}*/
