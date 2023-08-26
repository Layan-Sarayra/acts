#include "ActsExamples/Vertexing/KDEgenerator.hpp"
#include "VertexingHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TTree.h>

namespace ActsExamples {

KDEAlgorithm::KDEAlgorithm(const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("KDEAlgorithm", lvl), m_cfg(cfg){

    outFile = new TFile("OUTFILE.root", "RECREATE");

    // Open the ROOT file and access the TTree and TBranches
    const char* rootFilePath = "/eos/user/l/lalsaray/odd_output/tracksummary_ambi.root";

    TFile* inputFile = new TFile(rootFilePath);
    if (!inputFile || inputFile->IsZombie()) {
        ACTS_ERROR("Error opening ROOT file");
    }
    
    inputTree = dynamic_cast<TTree*>(inputFile->Get("tracksummary"));
    if (!inputTree) {
        ACTS_ERROR("Error getting TTree from ROOT file");
        inputFile->Close();
    }
    inputTree->SetMakeClass(1);

    //Set Object Pointer
    err_eLOC0_fit = 0;
    inputTree->SetBranchAddress("err_eLOC0_fit", &err_eLOC0_fit);
  }


ProcessCode KDEAlgorithm::execute(const AlgorithmContext& ctx) const {

  Long64_t nentries = inputTree->GetEntries();

  std::cout << "Number of events = " << nentries << std::endl; // this should give you 100 events. CHECK IT 

  for(Long64_t entry = 0; entry<nentries; entry++){
    Long64_t ientry = inputTree->LoadTree(entry);
    if(ientry < 0) break;
    inputTree->GetEntry(ientry);
  
    std::cout << "size of err_eLOC0_fit = " << err_eLOC0_fit->size() << std::endl;
    
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
 
} //namespace ActsExamples 