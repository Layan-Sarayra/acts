#include "ActsExamples/Vertexing/KDEgenerator.hpp"
#include "VertexingHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <TFile.h>
#include <TTree.h>

ActsExamples::KDEAlgorithm::KDEAlgorithm(const Config& cfg,
					 Acts::Logging::Level lvl)
  : ActsExamples::IAlgorithm("KDEAlgorithm", lvl), m_cfg(cfg) {}

ActsExamples::ProcessCode ActsExamples::KDEAlgorithm::execute(const ActsExamples::AlgorithmContext& ctx) const {
    
  ///ACTS_INFO(m_cfg.message);

  const char* rootFilePath = "/eos/user/l/lalsaray/odd_output/tracksummary_ambi.root";
  TFile* inputFile = TFile::Open(rootFilePath);
  if (!inputFile || inputFile->IsZombie()) {
    ACTS_ERROR("Error opening ROOT file");
    return ActsExamples::ProcessCode::ABORT;
  }

  TTree* inputTree = dynamic_cast<TTree*>(inputFile->Get("tracksummary"));
  if (!inputTree) {
    ACTS_ERROR("Error getting TTree from ROOT file");
    inputFile->Close();
    return ActsExamples::ProcessCode::ABORT;
  }

  std::vector<float> err_eLOC0_fit;
  TBranch* br_err_eLOC0_fit = inputTree->GetBranch("err_eLOC0_fit");
  br_err_eLOC0_fit->SetAddress(&err_eLOC0_fit);

  Long64_t nEntries = inputTree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    br_err_eLOC0_fit->GetEntry(i);
    ACTS_INFO("Entry " << i << ": err_eLOC0_fit = " << err_eLOC0_fit[0]);
  }

  inputFile->Close();

  return ActsExamples::ProcessCode::SUCCESS;
}
