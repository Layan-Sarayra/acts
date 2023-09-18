#pragma once // Ensure this file is included only once during compilation

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <vector>
#include <memory>
#include <string>

#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

namespace ActsExamples {
class KDEAlgorithm final : public IAlgorithm {
  public:
    struct Config {
      // Define configuration options here if needed
      std::string message;
    };

    // Constructor
    KDEAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

    // Destructor:
    ~KDEAlgorithm();

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
   std::vector<float> *z_0;
   std::vector<float> *err_eLOC0_fit;
   std::vector<float> *sigma_z_0;
   std::vector<float> *err_eLOC0LOC1_fit;
   mutable std::vector<double> sortedTracks; //mutable allows you to modify a member variable within the 'const' execute function.
   mutable std::vector<double> filteredTracks;

   double bandwidth;
   TH1F* kdeHistogram;

   mutable Long64_t ientry;
   mutable Long64_t eventNumber;

};

} // namespace ActsExamples