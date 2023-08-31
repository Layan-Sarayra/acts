#pragma once // Ensure this file is included only once during compilation

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <vector>
#include <memory>
#include <string>

#include <TTree.h>
#include <TFile.h>

namespace ActsExamples {

class KDEAlgorithm final : public IAlgorithm {
  public:
    struct Config {
      // Define configuration options here if needed
      std::string message;
    };

    // Constructor
    KDEAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

    // This function will be called on each event by the sequencer
    ProcessCode execute(const AlgorithmContext& ctx) const final;

    // Accessor for configuration
    const Config& config() const { return m_cfg; }

  
  private:
   Config m_cfg; // Configuration data for the algorithm
  

   //Define all member vars here
   TTree *inputTree;
   TFile *inputFile;
   TFile *outFile;
   std::vector<float> *eLOC0_fit;
   std::vector<float> *eLOC1_fit;
   std::vector<float> *err_eLOC0_fit;
   std::vector<float> *err_eLOC1_fit;
   std::vector<float> *err_eLOC0LOC1_fit;
};

} // namespace ActsExamples