#pragma once // Ensure this file is included only once during compilation

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <vector>
#include <memory>
#include <string>
#include <algorithm>

#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

namespace ActsExamples {

struct KDEData {
    double z0_candidate;
    double kdeValue;
};


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
    ProcessCode finalize() override;

    // Accessor for configuration
    const Config& config() const { return m_cfg; }

  
  private:
   Config m_cfg; // Configuration data for the algorithm
   mutable Long64_t ientry;
   mutable Long64_t eventNumber;

   //Define all member vars here

   TFile* inputFile = nullptr;
   TTree* inputTree = nullptr;
   TFile* outFile = nullptr;
   TH1F* kdeHistogram = nullptr; 
   
   std::vector<float> *d_0;
   std::vector<float> *z_0;
   std::vector<float> *sigma_d0;
   std::vector<float> *sigma_z0;
   std::vector<float> *sigma_d0_z0;

   mutable std::vector<double> sortedTracks; //mutable allows you to modify a member variable within the 'const' execute function.
   mutable std::vector<double> filteredTracks;
   mutable std::vector<KDEData> accumulatedData;

   mutable double z_min;
   mutable double z_max;


    //Define grid search parameters
   double bandwidth;
   int nbins;
    
};

} // namespace ActsExamples