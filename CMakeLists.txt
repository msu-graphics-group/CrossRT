#
# CrossRT can be used in two ways:
#
# 1) Inside LiteRT via litert_import(CrossRT): LiteScene and Embree are taken from LiteRT.
#
# 2) Standalone via add_subdirectory(<path to CrossRT>), then target_link_libraries(<app> CrossRT).
#    Embree 3 is taken from binaries stored near this file (lib, lib_win64, bin_win64, lib_macos),
#    or Embree 4 is searched with find_package(embree 4) if USE_EMBREE4=ON.
#    Headers of LiteMath (and LiteScene, vkutils, Vulkan for USE_VULKAN=ON) must be visible from parent
#    project (include_directories before add_subdirectory). Sources of LiteScene and vkutils are compiled
#    by parent project, because they are linked to the final executable.
#

if(COMMAND litert_import)
    set(CROSSRT_INSIDE_LITERT ON)
else()
    set(CROSSRT_INSIDE_LITERT OFF)
endif()

if(NOT DEFINED USE_VULKAN)
    set(USE_VULKAN ON)
endif()

if(CROSSRT_INSIDE_LITERT)
    litert_import(../LiteScene)
    include(dependencies/embree)
    set(CROSSRT_EMBREE_TARGET litert_embree)
else()
    if(USE_EMBREE4)
        find_package(embree 4 CONFIG REQUIRED)
        add_library(CrossRT_embree INTERFACE)
        target_link_libraries(CrossRT_embree INTERFACE embree)
        target_compile_definitions(CrossRT_embree INTERFACE USE_EMBREE4=1)
    else()
        add_library(CrossRT_embree SHARED IMPORTED)
        if(WIN32)
            set_target_properties(CrossRT_embree PROPERTIES
                IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/bin_win64/embree3.dll
                IMPORTED_IMPLIB   ${CMAKE_CURRENT_SOURCE_DIR}/lib_win64/embree3.lib)
        elseif(APPLE)
            set_target_properties(CrossRT_embree PROPERTIES
                IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/lib_macos/libembree3.dylib
                INTERFACE_LINK_LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR}/lib_macos/libtbb.dylib)
        else()
            set_target_properties(CrossRT_embree PROPERTIES
                IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/libembree3.so
                INTERFACE_LINK_LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR}/lib/libtbb.so)
        endif()
    endif()
    set(CROSSRT_EMBREE_TARGET CrossRT_embree)
endif()

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

target_sources(CrossRT PRIVATE
    ${SOURCES_EMBREE}
    ProxyRT.cpp
    ${SOURCES_SWRT}
    FactoryRT.cpp
)

if (USE_VULKAN)
    target_sources(CrossRT PRIVATE VulkanRTX.cpp)
endif()

target_include_directories(CrossRT PUBLIC .)
target_link_libraries(CrossRT PUBLIC ${CROSSRT_EMBREE_TARGET})

if(CROSSRT_INSIDE_LITERT)
    target_link_libraries(CrossRT PUBLIC LiteScene)
    if (USE_VULKAN)
        target_link_libraries(CrossRT PUBLIC LiteScene_SceneManager)
    endif()
else()
    target_sources(CrossRT PRIVATE builders/cbvh_stats.cpp)

    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(CrossRT PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()
