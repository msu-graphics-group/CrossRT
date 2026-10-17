#include "ProxyRT.h"

#include <cassert>

void RTX_Proxy::ClearGeom() 
{ 
  m_totalPrimsIS2 = 0;
  for(auto impl : m_imps)
    if(impl != nullptr) 
      impl->ClearGeom(); 

  m_remapTable.clear();
  m_offsetByTag.clear();
  m_geomTags.clear();
  m_remapGeom.clear();

  m_remapTable.reserve(1000);
  m_remapGeom.reserve(m_remapTable.capacity());
  m_geomTags.reserve(m_remapTable.capacity());
} 

uint32_t RTX_Proxy::AddGeom_Triangles3f(const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, uint32_t a_flags, size_t vByteStride) 
{
  uint32_t flagsFirst = (m_imps[1] == nullptr) ? a_flags : a_flags | BUILD_SKIP; // when both impl used, actual impl is HW, so skip build for SW impl, only use it's input data
  uint32_t res0 = m_imps[0]->AddGeom_Triangles3f(a_vpos3f, a_vertNumber, a_triIndices, a_indNumber, flagsFirst, vByteStride); 
  uint32_t res1 = uint32_t(-1);
  if(m_imps[1] != nullptr)
    res1 = m_imps[1]->AddGeom_Triangles3f(a_vpos3f, a_vertNumber, a_triIndices, a_indNumber, a_flags, vByteStride);
  
  m_totalPrimsIS2 += a_indNumber/3;
  m_remapGeom.push_back(LiteMath::uint2(res0, res1));

  return res0;
}
                               
void RTX_Proxy::UpdateGeom_Triangles3f(uint32_t a_geomId, const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, uint32_t a_flags, size_t vByteStride)
{
  uint32_t flagsFirst = (m_imps[1] == nullptr) ? a_flags : a_flags | BUILD_SKIP;
  m_imps[0]->UpdateGeom_Triangles3f(a_geomId, a_vpos3f, a_vertNumber, a_triIndices, a_indNumber, flagsFirst, vByteStride);
  if(m_imps[1] != nullptr)
    m_imps[1]->UpdateGeom_Triangles3f(a_geomId, a_vpos3f, a_vertNumber, a_triIndices, a_indNumber, a_flags,              vByteStride);
}
  
uint32_t RTX_Proxy::AddGeom_AABB(uint32_t a_typeId, const CRT_AABB* boxMinMaxF8, size_t a_boxNumber, uint32_t a_buildFlags, void** a_customPrimPtrs, size_t a_customPrimCount) 
{
  uint32_t flagsFirst = (m_imps[1] == nullptr) ? a_buildFlags : a_buildFlags | BUILD_SKIP;
  const uint32_t geomIdOrigin = m_imps[0]->AddGeom_AABB(a_typeId, boxMinMaxF8, a_boxNumber, flagsFirst, a_customPrimPtrs, a_customPrimCount);
  uint32_t geomId = uint32_t(-1);
  if(m_imps[1] != nullptr)
    geomId = m_imps[1]->AddGeom_AABB(a_typeId, boxMinMaxF8, a_boxNumber, a_buildFlags, a_customPrimPtrs, a_customPrimCount);
  else
    assert(false);
  
  //const uint32_t geomIdClean = geomId & (CRT_GEOM_MASK_AABB_BIT_RM);
  //if(m_imps[1] != nullptr)
  //{
  //  assert(m_geomTags.size() == geomIdClean); // m_geomTags[geomIdClean] = a_typeId;
  //}

  m_geomTags.push_back(a_typeId);           //
  m_remapGeom.push_back(LiteMath::uint2(geomIdOrigin, geomId));

  {
    size_t actualPrimsCount = (a_customPrimCount == 0) ? a_boxNumber : a_customPrimCount;
    uint32_t div            = uint32_t(a_boxNumber/actualPrimsCount);

    auto pRemapOffset = m_offsetByTag.find(a_typeId);
    if(pRemapOffset == m_offsetByTag.end()) 
      pRemapOffset = m_offsetByTag.insert(std::make_pair(a_typeId, LiteMath::uint2(0,div))).first;
    
    pRemapOffset->second.y = div;
    if(m_IS2Mode)
      pRemapOffset->second.x = m_totalPrimsIS2;
    
    m_remapTable.push_back(pRemapOffset->second);
    if(!m_IS2Mode)
      pRemapOffset->second.x += uint32_t(actualPrimsCount);

    m_totalPrimsIS2 += uint32_t(actualPrimsCount);
  }

  return geomIdOrigin;
}
  
void RTX_Proxy::UpdateGeom_AABB(uint32_t a_geomId, uint32_t a_typeId, const CRT_AABB* boxMinMaxF8, size_t a_boxNumber, uint32_t a_buildFlags, void** a_customPrimPtrs, size_t a_customPrimCount) 
{
  uint32_t flagsFirst = (m_imps[1] == nullptr) ? a_buildFlags : a_buildFlags | BUILD_SKIP;
  m_imps[0]->UpdateGeom_AABB(a_geomId, a_typeId, boxMinMaxF8, a_boxNumber, flagsFirst, a_customPrimPtrs, a_customPrimCount);
  if(m_imps[1] != nullptr)
    m_imps[1]->UpdateGeom_AABB(a_geomId, a_typeId, boxMinMaxF8, a_boxNumber, a_buildFlags, a_customPrimPtrs, a_customPrimCount);
}

RTX_Proxy::PrimitiveRemapTable RTX_Proxy::GetAABBToPrimTable() const
{
  PrimitiveRemapTable res;
  res.table     = m_remapTable.data();
  res.tableSize = uint32_t(m_remapTable.size());
  res.geomTable = m_geomTags.data();
  res.geomSize  = uint32_t(m_geomTags.size());
  return res;
}

void RTX_Proxy::ClearScene()  
{ 
  for(auto impl : m_imps) 
    if(impl != nullptr)
      impl->ClearScene(); 
} 

void RTX_Proxy::CommitScene(uint32_t options)  
{ 
  for(auto impl : m_imps) 
    if(impl != nullptr)
      impl->CommitScene(options); 
}

uint32_t RTX_Proxy::AddInstanceMotion(uint32_t a_geomId, const LiteMath::float4x4* a_matrices, uint32_t a_matrixNumber)  
{ 
  uint32_t res = 0;
  for(auto impl : m_imps) 
    if(impl != nullptr)
      res = impl->AddInstanceMotion(a_geomId, a_matrices, a_matrixNumber);
  return res; 
}

uint32_t RTX_Proxy::AddInstance(uint32_t a_geomId, const LiteMath::float4x4& a_matrix)  
{ 
  a_geomId = a_geomId & 0x7FFFFFFFu;
  assert(a_geomId < m_remapGeom.size());
  auto remap = m_remapGeom[a_geomId];
  assert(a_geomId == remap.x); // 

  uint32_t res0 = m_imps[0]->AddInstance(remap.x, a_matrix);
  if(m_imps[1] != nullptr)
    m_imps[1]->AddInstance(remap.y, a_matrix);
  return res0; 
}

void RTX_Proxy::UpdateInstance(uint32_t a_instanceId, const LiteMath::float4x4& a_matrix) 
{ 
  for(auto impl : m_imps)
    if(impl != nullptr) 
      impl->UpdateInstance(a_instanceId, a_matrix); 
}

CRT_Hit RTX_Proxy::RayQuery_NearestHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar)  
{ 
  return m_imps[0]->RayQuery_NearestHit(posAndNear, dirAndFar); 
}

bool RTX_Proxy::RayQuery_AnyHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar)  
{ 
  return m_imps[0]->RayQuery_AnyHit(posAndNear, dirAndFar);
}  

CRT_Hit RTX_Proxy::RayQuery_NearestHitMotion(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, float time)  
{ 
  return m_imps[0]->RayQuery_NearestHit(posAndNear, dirAndFar); 
}

bool RTX_Proxy::RayQuery_AnyHitMotion(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar, float time) 
{ 
  return m_imps[0]->RayQuery_AnyHit(posAndNear, dirAndFar); 
}

uint32_t RTX_Proxy::AddCustomGeom_FromFile(const char *geom_type_name, const char *filename, ISceneObject *fake_this)
{
  return m_imps[0]->AddCustomGeom_FromFile(geom_type_name, filename, fake_this);
}
