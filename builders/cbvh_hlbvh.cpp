#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "cbvh.h"
#include "cbvh_core.h"
#include "cbvh_internal.h"

#include "lbvh.h"

using LiteMath::float4;
using LiteMath::float3;
using LiteMath::float2;
using LiteMath::uint2;
using LiteMath::int4;
using LiteMath::uint3;
using LiteMath::sign;
using LiteMath::dot;
using cbvh::BVHNode;
using cbvh::Interval;

extern double g_buildTime;

namespace cbvh
{
  static constexpr uint32_t START_MASK = 0x00FFFFFF;
  static constexpr uint32_t END_MASK   = 0xFF000000;
  static constexpr uint32_t SIZE_MASK  = 0x7F000000;
  static constexpr uint32_t LEAF_BIT   = 0x80000000;
  
  static inline uint32_t PackOffsetAndSize(uint32_t start, uint32_t size)
  {
    return LEAF_BIT | ((size << 24) & SIZE_MASK) | (start & START_MASK);
  }
  
  static inline uint32_t EXTRACT_START(uint32_t a_leftOffset)  { return  a_leftOffset & START_MASK; }
  static inline uint32_t EXTRACT_COUNT(uint32_t a_leftOffset)  { return (a_leftOffset & SIZE_MASK) >> 24; }

  struct PrintDFS_LR
  {
    PrintDFS_LR(const std::vector<BVHNode>& a_nodes, std::ofstream& a_out) : nodes(a_nodes), out(a_out) {}
    const std::vector<BVHNode>& nodes; 
    std::ofstream& out;
  
    void VisitNode(uint32_t a_nodeOffset)
    {
      const auto& currNode = nodes[a_nodeOffset];
  
      // print node here
      out << a_nodeOffset << std::endl;
  
      if((currNode.leftOffset & LEAF_BIT) == 0) 
      {
        out << a_nodeOffset << " -> " << currNode.leftOffset << std::endl;
        out << a_nodeOffset << " -> " << currNode.escapeIndex << std::endl;
        VisitNode(currNode.leftOffset);
        VisitNode(currNode.escapeIndex);
      }
    }
  };

  struct PrintDFS_LO
  {
    PrintDFS_LO(const std::vector<BVHNode>& a_nodes, std::ofstream& a_out) : nodes(a_nodes), out(a_out) {}
    const std::vector<BVHNode>& nodes; 
    std::ofstream& out;
  
    void VisitNode(uint32_t a_leftOffset)
    {
      if((a_leftOffset & LEAF_BIT) != 0) 
        return;

      const auto& leftNode  = nodes[a_leftOffset + 0];
      const auto& rightNode = nodes[a_leftOffset + 1];
  
      // print node here
      out << a_leftOffset + 0 << std::endl;
      out << a_leftOffset + 1 << std::endl;
  
      if((leftNode.leftOffset & LEAF_BIT) == 0) 
      {
        out << a_leftOffset+0 << " -> " << leftNode.leftOffset + 0 << std::endl;
        out << a_leftOffset+0 << " -> " << leftNode.leftOffset + 1 << std::endl;
        VisitNode(leftNode.leftOffset);
      }

      if((rightNode.leftOffset & LEAF_BIT) == 0) 
      {
        out << a_leftOffset+1 << " -> " << rightNode.leftOffset + 0 << std::endl;
        out << a_leftOffset+1 << " -> " << rightNode.leftOffset + 1 << std::endl;
        VisitNode(rightNode.leftOffset);
      }
    }
  };

  void PrintBVH2ForGraphViz_LR(const std::vector<BVHNode>& a_nodes, const char* a_fileName)
  {
    std::ofstream fout(a_fileName);
    PrintDFS_LR trav(a_nodes, fout);
    fout << "digraph D {" << std::endl;
    fout << "rankdir=\"LR\"" << std::endl;
    trav.VisitNode(0);
    fout << "}" << std::endl;
  }

  void PrintBVH2ForGraphViz_LO(const std::vector<BVHNode>& a_nodes, const char* a_fileName)
  {
    std::ofstream fout(a_fileName);
    PrintDFS_LO trav(a_nodes, fout);
    fout << "digraph D {" << std::endl;
    fout << "rankdir=\"LR\"" << std::endl;
    trav.VisitNode(0);
    fout << "}" << std::endl;
  }

  template<class T>
  bool CompareBuffers(const std::vector<T>& a, const std::vector<T>& b)
  {
    if(a.size() != b.size())
      return false;

    for(size_t i=0;i<a.size();i++)
    {
      if(a[i] != b[i])
      {
        std::cout << "diff at " << i << ", " << a[i] << " != " << b[i] << std::endl;
        return false;
      }
    }

    return true;  
  }

  template<>
  bool CompareBuffers(const std::vector<uint2>& a, const std::vector<uint2>& b)
  {
    if(a.size() != b.size())
      return false;

    for(size_t i=0;i<a.size();i++)
    {
      if(a[i].x != b[i].x || a[i].y != b[i].y)
      {
        std::cout << "diff at " << i << ", (" << a[i].x << "," << a[i].y << ") != (" << b[i].x << "," << b[i].y << ")" << std::endl;
        return false;
      }
    }

    return true;  
  }

  bool CompareBuffersXOnly(const std::vector<uint2>& a, const std::vector<uint2>& b)
  {
    if(a.size() != b.size())
      return false;

    for(size_t i=0;i<a.size();i++)
    {
      if(a[i].x != b[i].x)
      {
        std::cout << "diff at " << i << ", (" << a[i].x << "," << a[i].y << ") != (" << b[i].x << "," << b[i].y << ")" << std::endl;
        return false;
      }
    }

    return true;
  }

  bool CompareLBVHNodesAfterKarras(const BVHNode* a, const BVHNode* b, uint32_t leavesNum, bool testLeavesOnly = false)
  {
    if(!testLeavesOnly){
      for(uint32_t i=0; i<leavesNum-1;i++){
        if(a[i].leftOffset != b[i].leftOffset || a[i].escapeIndex != b[i].escapeIndex){
          std::cout << "diff at " << i << " (" << a[i].leftOffset << "," << a[i].leftOffset << ") != ";
          std::cout << " (" << b[i].leftOffset << "," << b[i].leftOffset << ")" << std::endl;
          return false;
        }
      }
    }

    for(uint32_t i=leavesNum-1; i<2*leavesNum-1;i++){
      const float diffMin = LiteMath::length(a[i].boxMin - b[i].boxMin);
      const float diffMax = LiteMath::length(a[i].boxMax - b[i].boxMax);
      if(a[i].leftOffset != b[i].leftOffset || a[i].escapeIndex != b[i].escapeIndex || diffMin > 1e-5f || diffMax > 1e-5f){
        std::cout << "diff at " << i << " (" << a[i].leftOffset << "," << a[i].leftOffset << ") != ";
        std::cout << " (" << b[i].leftOffset << "," << b[i].leftOffset << ") | (" << diffMin << ", " << diffMax << ") " << std::endl;
        return false;
      }
    }
    return true;
  }

  const bool USE_PEACH_TREES = true;
  const bool DEBUG_BUILDER   = false;

  BVHTree BuildLBVH(const float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum, CBVH_FORMATS a_desiredFormat, bool onGPU)
  {
    const bool   useBoxes  = (a_indices == nullptr);
    const size_t boxNumber = useBoxes ? a_vertNum : a_indexNum/3;
    if(useBoxes)
      a_indexNum = 0;

    const cbvh2::FMT format = (USE_PEACH_TREES && a_desiredFormat == FMT_BVH2Node32_Interval32_Static) ? cbvh2::BVH2_LEFT_OFFSET : cbvh2::BVH2_LEFT_RIGHT;

    std::shared_ptr<LBVHBuilder> pBuilder = std::make_shared<LBVHBuilder>();

    size_t reservedNodes = pBuilder->Reserve(boxNumber, format); ///< normally you have to do it only once!
    pBuilder->CommitDeviceData();                                ///< normally you have to do it only once!

    BVHTree result; 
    result.nodes.resize(reservedNodes);
    result.intervals.resize(reservedNodes);
    result.indicesReordered.resize(boxNumber);
    
    if(useBoxes)
      pBuilder->BuildFromBoxes(a_vertices, boxNumber, 
                               result.nodes.data(), uint32_t(reservedNodes), result.indicesReordered.data());
    else
      pBuilder->BuildFromTriangles(a_vertices, a_vertNum, a_indices, a_indexNum, 
                                   result.nodes.data(), uint32_t(reservedNodes), result.indicesReordered.data());

    // convert to old format if needed
    //
    result.format       = cbvh::FMT_BVH2Node32_Interval32_Dynamic;
    result.leavesNumber = uint32_t(pBuilder->GetLeavesNumber());
    result.depthRanges.resize(0);
    //result.nodes.resize(pBuilder->GetNodesNumber());      // DON'T !!! (because of peach convert)
    //result.intervals.resize(pBuilder->GetNodesNumber());  // DON'T !!! (because of peach convert)

    if(result.nodes.size()%2 != 0) // should be always true for LBVH, but not for peach convert
    {
      result.nodes.push_back(cbvh2::DummyNode());
      result.intervals.push_back(cbvh::Interval(0,0));
    }

    const char* funcName = useBoxes ? "BuildFromBoxes" : "BuildFromTriangles";

    float timings[8] = {};
    pBuilder->GetExecutionTime(funcName,timings);

    std::cout << "[lbvh]: BoxNum     = " << boxNumber << std::endl;
    std::cout << "[lbvh]: 2*BoxNum-1 = " << 2*boxNumber-1 << std::endl;
    std::cout << "[lbvh]: NodesNum   = " << pBuilder->GetNodesNumber() << std::endl;
    
    if(false)
    {
      const float ev_codes_ms    = timings[0];
      const float sort_ms        = timings[1];
      const float form_leaves_ms = timings[2];
      const float karras_ms      = timings[3];
      const float refit_ms       = timings[4];
      const float convert_ms     = timings[5];
      const float total_ms       = ev_codes_ms + sort_ms + form_leaves_ms + karras_ms + refit_ms + convert_ms;

      std::cout << "[lbvh_cpu]: ev_codes     = " << std::fixed << std::setprecision(3) << ev_codes_ms    << " ms" << std::endl;
      std::cout << "[lbvh_cpu]: sort         = " << std::fixed << std::setprecision(3) << sort_ms        << " ms" << std::endl;
      std::cout << "[lbvh_cpu]: form_leaves  = " << std::fixed << std::setprecision(3) << form_leaves_ms << " ms" << std::endl;
      std::cout << "[lbvh_cpu]: karras_12    = " << std::fixed << std::setprecision(3) << karras_ms      << " ms" << std::endl;
      std::cout << "[lbvh_cpu]: refit        = " << std::fixed << std::setprecision(3) << refit_ms       << " ms" << std::endl;
      std::cout << "[lbvh_cpu]: convert      = " << std::fixed << std::setprecision(3) << convert_ms     << " ms" << std::endl;
      std::cout << "[lbvh_cpu]: total        = " << std::fixed << std::setprecision(3) << total_ms << " ms" << std::endl;

      g_buildTime += total_ms;
    }
    
    if(useBoxes)
    {
      if(USE_PEACH_TREES && a_desiredFormat == FMT_BVH2Node32_Interval32_Static)
      {
        result.format = cbvh::FMT_BVH2Node32_Interval32_Static;
        result.invalidIntervals = true;
        return result;
      }
      else if(a_desiredFormat == FMT_BVH2Node32_Interval32_Dynamic)
      {
        result.format = cbvh::FMT_BVH2Node32_Interval32_Dynamic;
        result.invalidIntervals = true;
        return result;
      }
    }   

    // we drop support for intervals in new builder, therefore now we have to extract them from leaf nodes
    //
    result.invalidIntervals = false;
    for(uint32_t i=0;i<result.intervals.size();i++) 
    {
      const uint32_t offset = result.nodes[i].leftOffset;
      if((offset & LEAF_BIT) != 0)
      {
        const uint32_t start = EXTRACT_START(offset);
        const uint32_t count = EXTRACT_COUNT(offset);
        result.intervals[i] = Interval(start, count);
      }
      else
        result.intervals[i] = Interval(0, 0);
    }   

    if(USE_PEACH_TREES && a_desiredFormat == FMT_BVH2Node32_Interval32_Static)
    {
      result.format = cbvh::FMT_BVH2Node32_Interval32_Static;
      return result;
    }
    else
    {
      if(a_desiredFormat == FMT_BVH2Node32_Interval32_Dynamic)
        return result;
      else if(a_desiredFormat == FMT_BVH2Node32_Interval32_Static)
      {
        return ConvertBVH2DynamicToBVH2Flat(result);
        //auto tmp = ConvertBVH2DynamicToBVH2Flat(result);
        //if(useBoxes)
        //  PrintBVH2ForGraphViz_LO(tmp.nodes, "z_boxes2.txt");
        //return tmp;
      }
      else
        return cbvh::ConvertBVH2DynamicToBVH4Flat(result);
    }
  }
}
