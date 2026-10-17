#include "FatBVH.h"
#include <queue>
#include <stack>
#include <iostream>
#include <fstream>
#include <cassert>
#include <unordered_set>

struct PrintDFS
{
  PrintDFS(const std::vector<BVHNodeFat>& a_nodes, std::ofstream& a_out) : nodes(a_nodes), out(a_out) {}
  
  const std::vector<BVHNodeFat>& nodes; 
  std::ofstream& out;

  void VisitNode(uint32_t a_nodeOffset)
  {
    const auto& currNode = nodes[a_nodeOffset];

    // print node here
    out << a_nodeOffset << std::endl;

    if(!(currNode.offs_left & LEAF_BIT)) {
      VisitNode(currNode.offs_left);
      out << a_nodeOffset << " -> " << currNode.offs_left << std::endl;
    }

    if(!(currNode.offs_right & LEAF_BIT)) {
      VisitNode(currNode.offs_right);
      out << a_nodeOffset << " -> " << currNode.offs_right << std::endl;
    }
  }

};

void FatBVH::PrintForGraphViz(const std::vector<BVHNodeFat>& a_nodes, const std::vector<int>& a_treeletRoots, const std::vector<int>& a_treeletRootsSuper, const char* a_fileName)
{
  std::ofstream fout(a_fileName);
  PrintDFS trav(a_nodes, fout);
  fout << "digraph D {" << std::endl;
  fout << "rankdir=\"LR\"" << std::endl;
  trav.VisitNode(0);
  
  if(a_treeletRoots.size() == 0) {
    fout << "}" << std::endl;
    return;
  }

  if(a_treeletRootsSuper.size() == 0) {
    for(size_t trId=0;trId<a_treeletRoots.size()-1;trId++) {
      int start = a_treeletRoots[trId];
      int end   = a_treeletRoots[trId+1];
      fout << "subgraph cluster_p" << trId <<  "{";
      for(int i=start; i<end;i++) 
        fout << i << ";";
      fout << "}" << std::endl;
    }
  }
  else {
    std::vector< std::vector<int> > superTreelsts(a_treeletRootsSuper.size());
    std::unordered_set<int> processedTreelets;
    for(size_t sId=0;sId<a_treeletRootsSuper.size()-1;sId++) {
      int start_s = a_treeletRootsSuper[sId+0];
      int end_s   = a_treeletRootsSuper[sId+1];
      for(size_t trId=0;trId<a_treeletRoots.size()-1;trId++) {
        int start = a_treeletRoots[trId+0];
        int end   = a_treeletRoots[trId+1];
        if(start >= start_s && end <= end_s) {
          superTreelsts[sId].push_back(trId);
          processedTreelets.insert(trId);
        }
      }
    }
    
    size_t clustId = 0;
    for(size_t j=0;j<superTreelsts.size();j++) {
      fout << "subgraph cluster_g" << j <<  "{" << std::endl;
      for(size_t trId=0;trId<superTreelsts[j].size();trId++,clustId++) {
        int trId2 = superTreelsts[j][trId];
        int start = a_treeletRoots[trId2+0];
        int end   = a_treeletRoots[trId2+1];
        fout << "subgraph cluster_p" << clustId <<  "{";
        for(int i=start; i<end;i++) 
          fout << i << ";";
        fout << "}" << std::endl;
      }
      fout << "}" << std::endl;
    }
    
    for(size_t trId=0;trId<a_treeletRoots.size()-1;trId++,clustId++) { // print treeelets that were not in 'superTreelsts' due to the last element process bug
      int start = a_treeletRoots[trId];
      int end   = a_treeletRoots[trId+1];
      if(processedTreelets.find(trId) == processedTreelets.end()) {
        fout << "subgraph cluster_p" << clustId <<  "{";
        for(int i=start; i<end;i++) 
          fout << i << ";";
        fout << "}" << std::endl;
      }
    }
      
  }

  fout << "}" << std::endl;
}