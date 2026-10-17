include(../litert.cmake)
litert_import(../LiteScene)
include(dependencies/embree)

add_library(CrossRT STATIC)

set(SOURCES_SWRT
    BVH2CommonRT.cpp
    BVH2CommonRT_host.cpp
    BruteForceRT.cpp
    builders/cbvh_core.cpp
    builders/cbvh.cpp
    builders/cbvh_embree2.cpp
    builders/cbvh_esc2.cpp
    builders/cbvh_esc.cpp
    builders/cbvh_fat.cpp
    builders/cbvh_hlbvh.cpp
    builders/cbvh_split.cpp
    builders/EXT_TriBoxOverlap.cpp
    builders/FatBVHClusterize.cpp
    builders/FatBVHPrint.cpp
    builders/FatBVHReorder.cpp
    builders/lbvh.cpp
    builders/lbvh_host.cpp
)

if (USE_EMBREE4)
    set(SOURCES_EMBREE EmbreeRT4.cpp)
else()
    set(SOURCES_EMBREE EmbreeRT.cpp)
endif()

if (USE_VULKAN)
    target_sources(CrossRT PRIVATE
        ${SOURCES_EMBREE}
        ProxyRT.cpp
        ${SOURCES_SWRT}
        FactoryRT.cpp
        VulkanRTX.cpp
    )

    target_include_directories(CrossRT PUBLIC .)
    target_link_libraries(CrossRT PUBLIC LiteScene LiteScene_SceneManager litert_embree)
else()
    target_sources(CrossRT PRIVATE
        ${SOURCES_EMBREE}
        ProxyRT.cpp
        ${SOURCES_SWRT}
        FactoryRT.cpp
    )
    target_include_directories(CrossRT PUBLIC .)
    target_link_libraries(CrossRT PUBLIC LiteScene litert_embree)
endif()
