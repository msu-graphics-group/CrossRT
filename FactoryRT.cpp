#include "CrossRT.h"
#include "ProxyRT.h"
#include <string>
#include <iostream>
#include <memory>

ISceneObject* CreateEmbreeRT(const char* a_implName);
ISceneObject* MakeBruteForceRT(const char* a_implName);
ISceneObject *MakeBVH2CommonRT(const char *a_implName, const char* a_buildName);

ISceneObject* CreateSceneRT(const char* a_implName)
{
  const std::string className(a_implName); //
  if (className.find("BruteForce") != std::string::npos)
    return MakeBruteForceRT(a_implName);
  else if (className.find("BVH2Common") != std::string::npos)
    return MakeBVH2CommonRT(a_implName, "cbvh_embree2");
  else
  {
    //std::cout << "[CreateSceneRT]: Use Embree for '" << className.c_str() << "' " << std::endl;
    return CreateEmbreeRT(a_implName);
  }
}

void DeleteSceneRT(ISceneObject* a_impl) { delete a_impl; }

//std::shared_ptr<ISceneObject> MakeSceneRT(const char* a_implName) // create proxy with a single implementation
//{
//  ISceneObject* pImplNaked = CreateSceneRT(a_implName);
//  std::shared_ptr<ISceneObject> pImpl(pImplNaked);
//  return std::make_shared<RTX_Proxy>(pImpl, nullptr);
//}
