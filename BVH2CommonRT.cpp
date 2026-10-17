#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>

#include "BVH2CommonRT.h"

using std::acos;
using std::atan2;
using std::pow;
using std::log;
using std::min;
using std::max;
using std::abs;

//Octahedral Normal Vectors (ONV) encoding https://jcgt.org/published/0003/02/01/
static float2 encode_normal(float3 v)
{
  float2 p = float2(v.x, v.y) * (1.0f / (abs(v.x) + abs(v.y) + abs(v.z)));
  float2 signNotZero = float2((p.x >= 0.0f) ? +1.0f : -1.0f, (p.y >= 0.0f) ? +1.0f : -1.0f);
  return (v.z <= 0.0f) ? ((1.0f - abs(float2(p.y, p.x))) * signNotZero) : p;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline bool Quadratic2(float A, float B, float C, float *t0, float *t1) 
{
  float discrim = B * B - 4.0f * A * C;
  if (discrim < 0.f) 
    return false;
  float rootDiscrim = std::sqrt(discrim);
  float floatRootDiscrim   = float(rootDiscrim);
  // Compute quadratic _t_ values
  float q;
  if ((float)B < 0.0f)
      q = -.5f * (B - floatRootDiscrim);
  else
      q = -.5f * (B + floatRootDiscrim);
  *t0 = q / A;
  *t1 = C / q;
  if ((float)*t0 > (float)*t1) 
  {
    // std::swap(*t0, *t1);
    float temp = *t0;
    *t0 = *t1;
    *t1 = temp;
  }
  return true;
}

static inline float3 myfaceforward(const float3 n, const float3 v) { return (dot(n, v) < 0.f) ? (-1.0f)*n : n; }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline float mandelbulb_sdf(float3 pos) 
{
  const float mandelbulb_power    = 8.0f;
  const int   mandelbulb_iter_num = 16;
  
	float3 z = pos;
	float dr = 1.0f;
	float r  = 0.0f;
	for (int i = 0; i < mandelbulb_iter_num ; i++)
	{
		r = length(z);
		if (r > 1.5f) 
      break;
		
		// convert to polar coordinates
		float theta = acos(z.z / r);
		float phi   = atan2(z.y, z.x);

		dr = pow( r, mandelbulb_power-1.0f)*mandelbulb_power*dr + 1.0f;
		
		// scale and rotate the point
		float zr = pow(r, mandelbulb_power);
		theta = theta*mandelbulb_power;
		phi   = phi*mandelbulb_power;
		
		// convert back to cartesian coordinates
		z = pos + zr*float3(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta));
	}
	return 0.5f*log(r)*r/dr;
}

static inline float ray_marching_sdf(const float3 ray_pos, const float3 ray_dir)
{
  const float epsilon = 0.0002f;
  const float2 boxHit = RayBoxIntersection2(ray_pos, SafeInverse(ray_dir), float3(-1,-1,-1), float3(+1,+1,+1));

	float depth = boxHit.x;
	int steps   = 0;
  float dist  = epsilon*2.0f;

	while(depth <= boxHit.y && dist > epsilon && steps < 100)
	{
		dist   = mandelbulb_sdf(ray_pos + depth*ray_dir);
		depth += dist;
		steps++;
	}

  if(depth > boxHit.y || steps > 100)
    depth = -1.0f;
   
	return depth;
}

static inline float3 estimate_normal(float3 z, float eps)
{
  const float3 z1 = z + float3(eps, 0, 0);
  const float3 z2 = z - float3(eps, 0, 0);
  const float3 z3 = z + float3(0, eps, 0);
  const float3 z4 = z - float3(0, eps, 0);
  const float3 z5 = z + float3(0, 0, eps);
  const float3 z6 = z - float3(0, 0, eps);

  const float dx = mandelbulb_sdf(z1) - mandelbulb_sdf(z2);
  const float dy = mandelbulb_sdf(z3) - mandelbulb_sdf(z4);
  const float dz = mandelbulb_sdf(z5) - mandelbulb_sdf(z6);
 
  return normalize(float3(dx, dy, dz) / (2.0f*eps));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

uint32_t BVH2CommonRT::IntersectAllPrimitivesInLeaf(float4 rayPosAndNear, float4 rayDirAndFar, CRT_LeafInfo info, CRT_Hit *pHit)
{
  const uint32_t geomIdType = m_geomIdByInstId[info.instId];
  const uint32_t geomId     = (geomIdType & GEOM_ID_MASK);
  const uint32_t geomType   = (geomIdType & GEOM_TP_MASK) >> GEOM_ID_SHFT;

  uint32_t hitTag = 0;

  if(geomType == GEOM_TYPE_TRIANGLE)
  {  
    const uint2 a_geomOffsets = m_geomOffsets[info.geomId]; // or geomId
    const uint32_t a_start    = info.aabbId; // pass through
    const uint32_t a_count    = info.primId; // pass through

    for (uint32_t triId = a_start; triId < a_start + a_count; triId++)
    {
      const uint32_t A = m_indices[a_geomOffsets.x + triId*3 + 0];
      const uint32_t B = m_indices[a_geomOffsets.x + triId*3 + 1];
      const uint32_t C = m_indices[a_geomOffsets.x + triId*3 + 2];
  
      const float3 A_pos = to_float3(m_vertPos[a_geomOffsets.y + A]);
      const float3 B_pos = to_float3(m_vertPos[a_geomOffsets.y + B]);
      const float3 C_pos = to_float3(m_vertPos[a_geomOffsets.y + C]);
  
      const float3 edge1 = B_pos - A_pos;
      const float3 edge2 = C_pos - A_pos;
      const float3 pvec = cross(to_float3(rayDirAndFar), edge2);
      const float3 tvec = to_float3(rayPosAndNear) - A_pos;
      const float3 qvec = cross(tvec, edge1);
  
      const float invDet = 1.0f / dot(edge1, pvec);
      const float v = dot(tvec, pvec) * invDet;
      const float u = dot(qvec, to_float3(rayDirAndFar)) * invDet;
      const float t = dot(edge2, qvec) * invDet;
  
      if (v >= -1e-6f && u >= -1e-6f && (u + v <= 1.0f + 1e-6f) && t > rayPosAndNear.w && t < pHit->t) // if (v > -1e-6f && u > -1e-6f && (u + v < 1.0f+1e-6f) && t > tMin && t < hit.t)
      {
        pHit->t = t;
        pHit->primId = triId;
        pHit->instId = info.instId;
        pHit->geomId = info.geomId | (GEOM_TYPE_TRIANGLE << GEOM_ID_BITS);
        pHit->coords[0] = u;
        pHit->coords[1] = v;
        hitTag = GEOM_TYPE_TRIANGLE;
      }
    }

  } // triangles
  #ifdef ENABLE_SPHERES
  else if(geomType == GEOM_TYPE_SPHERE)
  {
    constexpr float  radius = 1.0f;
    //constexpr float3 center = float3(0,0,0);
  
    // Compute _t0_ and _t1_ for ray--element intersection
    const float3 o = to_float3(rayPosAndNear); // - center;
    const float  A = rayDirAndFar.x * rayDirAndFar.x + rayDirAndFar.y * rayDirAndFar.y + rayDirAndFar.z * rayDirAndFar.z;
    const float  B = 2 * (rayDirAndFar.x * o.x + rayDirAndFar.y * o.y + rayDirAndFar.z * o.z);
    const float  C = o.x * o.x + o.y * o.y + o.z * o.z - radius * radius;
    float  t0, t1;
    if (!Quadratic2(A, B, C, &t0, &t1)) 
      return 0;
    
    const float tHit = std::min(t0, t1);
    
    if(tHit > rayPosAndNear.w && tHit < pHit->t)
    {
      pHit->t         = tHit;
      pHit->primId    = info.instId;
      pHit->instId    = info.instId;
      pHit->geomId    = info.geomId | (GEOM_TYPE_SPHERE << GEOM_ID_BITS);
    
      hitTag = GEOM_TYPE_SPHERE;
    }
  }
  #endif
  #ifdef ENABLE_MANDELBULB
  else if(geomType == GEOM_TYPE_MANDELBULB)
  {
    float tHit = ray_marching_sdf(to_float3(rayPosAndNear), to_float3(rayDirAndFar));
    if(tHit > rayPosAndNear.w && tHit < pHit->t)
    {
      pHit->t         = tHit;
      pHit->primId    = info.instId;
      pHit->instId    = info.instId;
      pHit->geomId    = info.geomId | (GEOM_TYPE_MANDELBULB << GEOM_ID_BITS);

      hitTag = GEOM_TYPE_MANDELBULB;
    }
  }
  #endif
  return hitTag;
}      

void BVH2CommonRT::IntersectionComplete(float4 rayPosAndNear, float4 rayDirAndFar, CRT_Hit *pHit)
{
  if(pHit->instId != uint32_t(-1) && pHit->geomId != uint32_t(-1))
  {
    const uint32_t geomIdType = m_geomIdByInstId[pHit->instId];
    const uint32_t geomId     = (geomIdType & GEOM_ID_MASK);
    const uint32_t geomType   = (geomIdType & GEOM_TP_MASK) >> GEOM_ID_SHFT;
  
    const float3 ray_pos = matmul4x3(m_instMatricesInv[pHit->instId], to_float3(rayPosAndNear));
    const float3 ray_dir = matmul3x3(m_instMatricesInv[pHit->instId], to_float3(rayDirAndFar));

    if(geomType == GEOM_TYPE_TRIANGLE) // remap primitive id only for triangles
    {
      const uint2 geomOffsets = m_geomOffsets[pHit->geomId];
      pHit->primId = m_primIndices[geomOffsets.x/3 + pHit->primId];
    }
    #ifdef ENABLE_SPHERES
    else if(geomType == GEOM_TYPE_SPHERE)
    {
      // Compute surface normal of element at ray intersection point
      float3 norm = normalize(ray_pos + ray_dir*pHit->t);
             norm = myfaceforward(norm, -1.0f*to_float3(rayDirAndFar));
       
      float2 normEncoded = encode_normal(norm);
      pHit->coords[0] = normEncoded.x;
      pHit->coords[1] = normEncoded.y;
    }
    #endif
    #ifdef ENABLE_MANDELBULB
    else if(geomType == GEOM_TYPE_MANDELBULB)
    {
      float3 norm = estimate_normal(ray_pos + ray_dir*pHit->t, 1e-4f);
             norm = myfaceforward(norm, -1.0f*to_float3(rayDirAndFar));

      float2 normEncoded = encode_normal(norm);
      pHit->coords[0] = normEncoded.x;
      pHit->coords[1] = normEncoded.y;
    }
    #endif
  }
}

//extern bool g_debugPrint;

CRT_Hit BVH2CommonRT::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
{
  bool stopOnFirstHit = (dirAndFar.w <= 0.0f);
  if(stopOnFirstHit)
    dirAndFar.w *= -1.0f;

  #ifdef ENABLE_METRICS
  ResetVarLC();
  m_stats.raysNumber++;
  #endif

  [[threadlocal]] uint32_t stack[STACK_SIZE];

  CRT_Hit hit;
  hit.t      = dirAndFar.w;
  hit.primId = uint32_t(-1);
  hit.instId = uint32_t(-1);
  hit.geomId = uint32_t(-1);
  
  ////////////////////////////////////////////////////////////////////////////// instancing variables are a bit more than common
  int instTop  = 0;
  uint32_t instId = uint32_t(-1);
  
  float3 ray_pos = to_float3(posAndNear);
  float3 ray_dir = to_float3(dirAndFar);
  float3 inv_dir = SafeInverse(ray_dir);

  uint32_t bvhOffset = m_tlasOffset;
  ////////////////////////////////////////////////////////////////////////////// instancing variables are a bit more than common

  int top = 0;
  uint32_t leftNodeOffset = 0; 

  while (top >= 0 && !(stopOnFirstHit && hit.primId != uint32_t(-1)))
  {
    while (top >= 0 && ((leftNodeOffset & LEAF_BIT) == 0))
    {
      //if(g_debugPrint)
      //  std::cout << "leftNodeOffset = " << leftNodeOffset << std::endl;

      #ifdef ENABLE_METRICS
      m_stats.NC  += 2;
      m_stats.BLB += 2*sizeof(BVHNode);
      #endif
      
      //auto nodeAddrB = reinterpret_cast<uintptr_t>(m_allNodes.data() + bvhOffset);
      //auto nodeAddr  = reinterpret_cast<uintptr_t>(m_allNodes.data() + bvhOffset + leftNodeOffset);
      //assert(nodeAddrB % 64 == 0);
      //assert(nodeAddr  % 64 == 0);

      const BVHNode node0 = m_allNodes[bvhOffset + leftNodeOffset + 0];
      const BVHNode node1 = m_allNodes[bvhOffset + leftNodeOffset + 1];

      const float2 tm0 = RayBoxIntersection2(ray_pos, inv_dir, node0.boxMin, node0.boxMax);
      const float2 tm1 = RayBoxIntersection2(ray_pos, inv_dir, node1.boxMin, node1.boxMax);

      const bool hitChild0 = (tm0.x <= tm0.y) && (tm0.y >= posAndNear.w) && (tm0.x <= hit.t);
      const bool hitChild1 = (tm1.x <= tm1.y) && (tm1.y >= posAndNear.w) && (tm1.x <= hit.t);

      //if(g_debugPrint)
      //{
      //  std::cout << "node0.boxMin = (" << node0.boxMin.x << ", " << node0.boxMin.y << ", " << node0.boxMin.z << ")" << std::endl;
      //  std::cout << "node0.boxMax = (" << node0.boxMax.x << ", " << node0.boxMax.y << ", " << node0.boxMax.z << ")" << std::endl;
      //}

      // traversal decision
      //
      leftNodeOffset = hitChild0 ? node0.leftOffset : node1.leftOffset;

      if (hitChild0 && hitChild1)
      {
        leftNodeOffset = (tm0.x <= tm1.x) ? node0.leftOffset : node1.leftOffset; // GPU style branch
        stack[top]     = (tm0.x <= tm1.x) ? node1.leftOffset : node0.leftOffset; // GPU style branch
        top++;
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
      }

      if (!hitChild0 && !hitChild1) // both miss, stack.pop()
      {
        top--;
        leftNodeOffset = stack[std::max(top,0)];
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
      }

      if (top >= 0 && top < instTop && bvhOffset != m_tlasOffset)
      {
        ray_pos   = to_float3(posAndNear);
        ray_dir   = to_float3(dirAndFar);
        inv_dir   = SafeInverse(ray_dir);
        bvhOffset = m_tlasOffset;
      }

    } // end while (searchingForLeaf)

    if (top >= 0 && leftNodeOffset != 0xFFFFFFFF && bvhOffset != m_tlasOffset)  // leaf node of BLAS, intersect triangles
    {
      const uint32_t start = EXTRACT_START(leftNodeOffset);
      const uint32_t count = EXTRACT_COUNT(leftNodeOffset);
     
      CRT_LeafInfo leafInfo;
      leafInfo.aabbId = start; // pass-through, bad code, wrong usage!
      leafInfo.primId = count; // pass-through, bad code, wrong usage!
      leafInfo.instId = instId;
      leafInfo.geomId = m_geomIdByInstId[instId] & GEOM_ID_MASK;
      leafInfo.rayxId = 0; 
      leafInfo.rayyId = 0; 
      
      const float4 rayPosAndNear2 = to_float4(ray_pos, posAndNear.w);
      const float4 rayDirAndFar2  = to_float4(ray_dir, dirAndFar.w);
      IntersectAllPrimitivesInLeaf(rayPosAndNear2, rayDirAndFar2, leafInfo, &hit); 

      #ifdef ENABLE_METRICS
      m_stats.LC++;
      m_stats.LC2++;
      m_stats.TC+=count;
      m_stats.BLB+=count*(3*sizeof(uint32_t) + 3*sizeof(float4));
      #endif
    }
    else if (top >= 0 && bvhOffset == m_tlasOffset)                            // leaf node of BLAS, intersect BLAS next
    {
      instId    = EXTRACT_START(leftNodeOffset);
      bvhOffset = m_bvhOffsets[m_geomIdByInstId[instId] & GEOM_ID_MASK];
      
      //if(g_debugPrint)
      //  std::cout << "TLAS -> BLAS for inst(" << instId << ") at " << leftNodeOffset << std::endl;
        
      leftNodeOffset = 0;

      ray_pos = matmul4x3(m_instMatricesInv[instId], to_float3(posAndNear));
      ray_dir = matmul3x3(m_instMatricesInv[instId], to_float3(dirAndFar));
      inv_dir = SafeInverse(ray_dir);
      instTop = top;
    }
    
    if(top >= 0 && (leftNodeOffset != 0 || leftNodeOffset == 0xFFFFFFFF)) // continue BVH traversal, except the cases when we have just enter from TLAS to BLAS (see upper code)
    {                                                                     // (leftNodeOffset == 0xFFFFFFFF); // for empty leaves
      top--;
      leftNodeOffset = stack[std::max(top,0)];

      if (top < instTop && bvhOffset != m_tlasOffset)
      {
        ray_pos   = to_float3(posAndNear);
        ray_dir   = to_float3(dirAndFar);
        inv_dir   = SafeInverse(ray_dir);
        bvhOffset = m_tlasOffset;
      }
      
      #ifdef ENABLE_METRICS
      m_stats.SOC++;
      m_stats.SBL+=sizeof(uint32_t); 
      #endif
    }
    
  } // end while (top >= 0)
  
  IntersectionComplete(posAndNear, dirAndFar, &hit);

  return hit;
}

bool BVH2CommonRT::RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar)
{
  dirAndFar.w *= -1.0f;
  CRT_Hit hit = RayQuery_NearestHit(posAndNear, dirAndFar);
  return (hit.geomId != uint32_t(-1));
}
