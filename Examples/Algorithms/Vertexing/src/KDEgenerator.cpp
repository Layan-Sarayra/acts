#include "ActsExamples/Vertexing/KDEgenerator.hpp"
#include "VertexingHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <TFile.h>
#include <TTree.h>

using namespace std;

namespace ActsExamples {

KDEAlgorithm::KDEAlgorithm(const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("KDEAlgorithm", lvl), m_cfg(cfg) {}

ProcessCode KDEAlgorithm::execute(const AlgorithmContext& ctx) const {

    ///ACTS_INFO(m_cfg.message);
  const char* rootFilePath = "/eos/user/l/lalsaray/odd_output/tracksummary_ambi.root";
  TFile* inputFile = TFile::Open(rootFilePath);
  // if (!inputFile || inputFile->IsZombie()) {
  //     ACTS_ERROR("Error opening ROOT file");
  //     return ActsExamples::ProcessCode::ABORT;
  // }

  TTree* inputTree = dynamic_cast<TTree*>(inputFile->Get("tracksummary"));
  // if (!inputTree) {
  //     ACTS_ERROR("Error getting TTree from ROOT file");
  //     inputFile->Close();
  //     return ActsExamples::ProcessCode::ABORT;
  // }

  std::vector<float> err_eLOC0_fit;
  TBranch* br_err_eLOC0_fit = inputTree->GetBranch("err_eLOC0_fit");
  br_err_eLOC0_fit->SetAddress(&err_eLOC0_fit);

  Long64_t entry = 0;  // Change this to the desired event number
  inputTree->GetEntry(entry);

  ACTS_INFO("Entry " << entry << ": err_eLOC0_fit = " << err_eLOC0_fit[0]);

  inputFile->Close();

  return ProcessCode::SUCCESS;
}

}

std::array<Kernel,2> const ActsExamples::KDEAlgorithm :: getKDE_from_PDF() const{
  std::array<Kernel,2> kernels;



  return kernels;    
}
