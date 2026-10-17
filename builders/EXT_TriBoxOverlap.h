#define EXT_TRIBOXOVERLAP_GUARDIAN

namespace EXTERNAL_TOOLS
{
  int triBoxOverlap(float *boxcenter,float *boxhalfsize,float triverts[3][3]);
  int planeBoxOverlap(float normal[3],float d, float maxbox[3]);
};
