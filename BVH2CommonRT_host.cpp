#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>

#include "BVH2CommonRT.h"

#include "builders/cbvh_core.h"
//#include "aligned_alloc.h"
using cbvh::BVHNode;
using cbvh::Interval;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

constexpr size_t reserveSize = 1000;

void BVH2CommonRT::ClearGeom()
{
  m_vertPos.reserve(std::max<size_t>(100000, m_vertPos.capacity()));
  m_indices.reserve(std::max<size_t>(100000 * 3, m_indices.capacity()));
  m_primIndices.reserve(std::max<size_t>(100000, m_primIndices.capacity()));

  m_vertPos.resize(0);
  m_indices.resize(0);
  m_primIndices.resize(0);

  m_allNodes.reserve(std::max<size_t>(100000, m_allNodes.capacity()));
  m_allNodes.resize(0);

  m_geomOffsets.reserve(std::max(reserveSize, m_geomOffsets.capacity()));
  m_geomOffsets.resize(0);

  m_geomBoxes.reserve(std::max<size_t>(reserveSize, m_geomBoxes.capacity()));
  m_geomBoxes.resize(0);

  m_bvhOffsets.reserve(m_geomBoxes.capacity() + 1);
  m_bvhOffsets.resize(0);
  
  totalTrisMem = 0;
  totalTrisMem = 0;
  totalinstNumInst  = 0;
  m_tlasSize        = 0;
  m_tlasOffset      = 0;
  ClearScene();
}

size_t BVH2CommonRT::AppendTreeData(const std::vector<cbvh::BVHNode>& a_nodes, const std::vector<uint32_t>& a_indices, 
                                        const uint32_t *a_triIndices, size_t a_indNumber)
{
  const size_t oldSize = m_allNodes.size();
  size_t oldIndexSize  = m_indices.size();

  m_allNodes.insert(m_allNodes.end(), a_nodes.begin(), a_nodes.end());
  m_primIndices.insert(m_primIndices.end(), a_indices.begin(), a_indices.end());

  m_indices.resize(oldIndexSize + a_indices.size()*3);
  for(size_t i=0;i<a_indices.size();i++)
  {
    const uint32_t triId = a_indices[i];
    m_indices[oldIndexSize + 3*i+0] = a_triIndices[triId*3+0];
    m_indices[oldIndexSize + 3*i+1] = a_triIndices[triId*3+1];
    m_indices[oldIndexSize + 3*i+2] = a_triIndices[triId*3+2];
  }
  
  return oldSize;
}

uint32_t BVH2CommonRT::AddGeom_Triangles3f(const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, uint32_t a_qualityLevel, size_t vByteStride)
{
  const size_t vStride = vByteStride / 4;
  assert(vByteStride % 4 == 0);

  const uint32_t currGeomId = uint32_t(m_geomOffsets.size());
  const size_t oldSizeVert  = m_vertPos.size();
  const size_t oldSizeInd   = m_indices.size();

  m_geomOffsets.push_back(uint2(oldSizeInd, oldSizeVert));
  m_vertPos.resize(oldSizeVert + a_vertNumber);

  Box4f bbox;
  for (size_t i = 0; i < a_vertNumber; i++)
  {
    const float4 v = float4(a_vpos3f[i * vStride + 0], a_vpos3f[i * vStride + 1], a_vpos3f[i * vStride + 2], 1.0f);
    m_vertPos[oldSizeVert + i] = v;
    bbox.include(v);
  }

  // Build BVH for each geom and append it to big buffer
  //
  size_t oldBvhSize = 0;
  if((a_qualityLevel & BUILD_SKIP) == 0)
  {
    auto presets = cbvh2::BuilderPresetsFromString(m_builderName.c_str());
    auto bvhData = cbvh2::BuildBVH((const float*)(m_vertPos.data() + oldSizeVert), a_vertNumber, 16, a_triIndices, a_indNumber, presets);
  
    oldBvhSize = AppendTreeData(bvhData.nodes, bvhData.indices, a_triIndices, a_indNumber);
  }
  m_bvhOffsets.push_back(uint32_t(oldBvhSize));

  totalTrisMem += (a_indNumber/3);
  //m_geomSize.push_back(uint32_t(a_indNumber/3));

  bbox.boxMin.w = LiteMath::as_float(uint32_t(a_indNumber/3)); // store 'm_geomSize[geomId]';
  bbox.boxMax.w = LiteMath::as_float(GEOM_TYPE_TRIANGLE);      // store 'm_geomTags[geomId]'; Triangles are always have zero tag
  m_geomBoxes.push_back(bbox);

  return currGeomId;
}

void BVH2CommonRT::UpdateGeom_Triangles3f(uint32_t a_geomId, const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, uint32_t a_qualityLevel, size_t vByteStride)
{
  std::cout << "[BVH2CommonRT::UpdateGeom_Triangles3f]: " << "not implemeted!" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

uint32_t BVH2CommonRT::AddGeom_AABB(uint32_t a_typeId, const CRT_AABB* boxMinMaxF8, size_t a_boxNumber, uint32_t a_buildFlags, void** a_customPrimPtrs, size_t a_customPrimCount)
{
  Box4f bbox;
  for (size_t i = 0; i < a_boxNumber; i++) // TODO: may omit this loop, take two first bvh nodes
  {
    Box4f currBox(boxMinMaxF8[i].boxMin, boxMinMaxF8[i].boxMax);
    bbox.include(currBox);
  }
  
  if((a_buildFlags & BUILD_SKIP) == 0)
  {
    auto presets = cbvh2::BuilderPresetsFromString(m_builderName.c_str());
    auto bvhData = cbvh2::BuildBVH( (const cbvh2::BVHNode*)boxMinMaxF8, a_boxNumber, presets);
    
    m_bvhOffsets.push_back(uint32_t(m_allNodes.size()));
    m_allNodes.insert(m_allNodes.end(), bvhData.begin(), bvhData.end());
  }

  bbox.boxMin.w = LiteMath::as_float(uint32_t(1)); // store 'm_geomSize[geomId]';
  bbox.boxMax.w = LiteMath::as_float(a_typeId);    // store 'm_geomTags[geomId]'; Triangles are always have zero tag
  m_geomBoxes.push_back(bbox);

  const uint32_t currGeomId = uint32_t(m_geomOffsets.size());
  m_geomOffsets.push_back(uint2(0, 0));

  return currGeomId;
}

uint32_t BVH2CommonRT::AddCustomGeom_FromFile(const char *geom_type_name, const char *filename, ISceneObject *fake_this)
{
  if(std::string(geom_type_name) == "sphere")
  {
    CRT_AABB box;
    box.boxMin = float4(-1,-1,-1,0);
    box.boxMax = float4(+1,+1,+1,0);
    return fake_this->AddGeom_AABB(GEOM_TYPE_SPHERE, &box, 1);
  }
  else if(std::string(geom_type_name) == "mandelbulb")
  {
    CRT_AABB box;
    box.boxMin = float4(-1,-1,-1,0);
    box.boxMax = float4(+1,+1,+1,0);
    return fake_this->AddGeom_AABB(GEOM_TYPE_MANDELBULB, &box, 1);
  }
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BVH2CommonRT::ClearScene()
{
  m_instBoxes.reserve(std::max(reserveSize, m_instBoxes.capacity()));
  m_instMatricesInv.reserve(std::max(reserveSize, m_instMatricesInv.capacity()));
  m_instMatricesFwd.reserve(std::max(reserveSize, m_instMatricesFwd.capacity()));
  m_geomIdByInstId.reserve(std::max(reserveSize, m_geomIdByInstId.capacity()));

  m_instBoxes.resize(0);
  m_instMatricesInv.resize(0); // 
  m_instMatricesFwd.resize(0);
  m_geomIdByInstId.resize(0);

  // TLAS is kept in 'm_allNodes' until the next CommitScene() which replaces it: BLAS of geometry added
  // after the previous CommitScene() are located after the TLAS and must not be discarded here

  totalTrisInst    = 0;
  totalinstNumInst = 0;
  m_sceneBBox = CRT_AABB();
}

void BVH2CommonRT::CommitScene(uint32_t a_qualityLevel)
{

  //cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_OFFSET, cbvh2::BVH_CONSTRUCT_FAST, 1};
  cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_OFFSET, cbvh2::BVH_CONSTRUCT_MEDIUM, 1};
  if(a_qualityLevel == BUILD_LOW && m_tlasSize == 0)
    presets = {cbvh2::BVH2_LEFT_OFFSET, cbvh2::BVH_CONSTRUCT_FAST, 1}; // BUGGY a bit for some scenes, may loose instances. So only use when we need fast update.

  auto nodesTLAS = cbvh2::BuildBVH((const cbvh::BVHNode*)m_instBoxes.data(), m_instBoxes.size(), presets);

  // TLAS is stored at the end of 'm_allNodes'. If the previous TLAS is still at the end (no geometry was added after it),
  // it is replaced in place; otherwise BLAS of new geometry follow it, the previous TLAS stays unused and the new one is appended.
  // So BLAS data is only appended, and only 'm_allNodes[m_tlasOffset, m_tlasOffset + m_tlasSize)' is rewritten (Update_m_allNodes).
  //
  if(m_tlasSize != 0 && size_t(m_tlasOffset) + m_tlasSize == m_allNodes.size())
    m_allNodes.resize(m_tlasOffset);

  m_allNodes.reserve(m_allNodes.size() + 4*nodesTLAS.size()); // reserve for future updates
  m_tlasOffset = uint32_t(m_allNodes.size());
  m_tlasSize   = nodesTLAS.size();
  m_allNodes.insert(m_allNodes.end(), nodesTLAS.begin(), nodesTLAS.end());
 
  this->Update_m_allNodes(m_tlasOffset, m_tlasSize);
  this->Update_m_instMatricesInv(0, m_instMatricesInv.size());

  //std::cout << "[BVH2CommonRT]: totalBVHNodes  = " << m_allNodes.size() << std::endl;
  //std::cout << "[BVH2CommonRT]: totalInstances = " << nodesTLAS.size()  << std::endl;
}

uint32_t BVH2CommonRT::AddInstance(uint32_t a_geomId, const float4x4 &a_matrix)
{
  if((a_geomId & CRT_GEOM_MASK_AABB_BIT) != 0)
    a_geomId = (a_geomId & CRT_GEOM_MASK_AABB_BIT_RM);
    
  const auto &box = m_geomBoxes[a_geomId];

  // (1) mult mesh bounding box vertices with matrix to form new bouding box for instance
  float4 boxVertices[8]{
      a_matrix * float4{box.boxMin.x, box.boxMin.y, box.boxMin.z, 1.0f},

      a_matrix * float4{box.boxMax.x, box.boxMin.y, box.boxMin.z, 1.0f},
      a_matrix * float4{box.boxMin.x, box.boxMax.y, box.boxMin.z, 1.0f},
      a_matrix * float4{box.boxMin.x, box.boxMin.y, box.boxMax.z, 1.0f},

      a_matrix * float4{box.boxMax.x, box.boxMax.y, box.boxMin.z, 1.0f},
      a_matrix * float4{box.boxMax.x, box.boxMin.y, box.boxMax.z, 1.0f},
      a_matrix * float4{box.boxMin.x, box.boxMax.y, box.boxMax.z, 1.0f},

      a_matrix * float4{box.boxMax.x, box.boxMax.y, box.boxMax.z, 1.0f},
  };

  Box4f newBox;
  for (size_t i = 0; i < 8; i++)
    newBox.include(boxVertices[i]);
  
  m_sceneBBox.include(newBox);

  // (2) append bounding box and matrices
  //
  const uint32_t oldSize  = uint32_t(m_instBoxes.size());
  const uint32_t geomSize = LiteMath::as_uint(m_geomBoxes[a_geomId].boxMin.w); // m_geomSize[a_geomId];
  const uint32_t geomTags = LiteMath::as_uint(m_geomBoxes[a_geomId].boxMax.w); // m_geomTags[a_geomId];

  m_instBoxes.push_back(newBox);
  m_instMatricesFwd.push_back(a_matrix);
  m_instMatricesInv.push_back(inverse4x4(a_matrix));
  
  const uint32_t geomRemap = (a_geomId & GEOM_ID_MASK) | (geomTags << GEOM_ID_SHFT);
  m_geomIdByInstId.push_back(geomRemap);
  
  totalinstNumInst++;
  totalTrisInst += geomSize; 
  return oldSize;
}

void BVH2CommonRT::UpdateInstance(uint32_t a_instanceId, const float4x4 &a_matrix)
{
  std::cout << "[BVH2CommonRT::UpdateInstance]: " << "not implemeted!" << std::endl;
}

ISceneObject *MakeBVH2CommonRT(const char *a_implName, const char* a_buildName) { return new BVH2CommonRT(a_buildName); }