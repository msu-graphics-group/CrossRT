#include <cstdint>

// Build statistics, accumulated by BVH builders.
// Compiled only when CrossRT is used standalone: inside LiteRT they are defined in utils/raytracing/cbvh_embree2.cpp
//
double   g_buildTime = 0.0; // used to get build time from deep underground of code
uint64_t g_buildTris = 0;   // used to get total tris processed by builder
