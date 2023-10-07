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

    // This function will run after all events are processed
    ProcessCode finalize() override;

    void copyBranches();

    // Accessor for configuration
    const Config& config() const { return m_cfg; }

  
  private:
   Config m_cfg; // Configuration data for the algorithm

   mutable Long64_t ientry;
   mutable Long64_t performance_entry;
   mutable Long64_t eventNumber;

   //Define all member vars here

   // first input ROOT file, aka tracksummary_ambi.root
   TFile* inputFile = nullptr;
   TTree* inputTree = nullptr;

   std::vector<float> *d_0;
   std::vector<float> *z_0;
   std::vector<float> *sigma_d0;
   std::vector<float> *sigma_z0;
   std::vector<float> *sigma_d0_z0;


   // second input ROOT file, aka performance_vertexing.root
   TFile* performanceFile = nullptr;
   TTree* performanceTree = nullptr;

   std::vector<float> *truthX;
   std::vector<float> *truthY;
   std::vector<float> *truthZ;
   std::vector<float> *recoX;
   std::vector<float> *recoY;
   std::vector<float> *recoZ;


   // defining output ROOT file, aka KDE_output_file.root
   TFile* outFile = nullptr;
   TTree* outputTree = nullptr;
   // TH1F* kdeHistogram = nullptr;
   // TH1F* recoZHistogram = nullptr;

   
   //from the 1st input root file
    std::vector<float> *m_recoTrack_z0;
    std::vector<float> *m_recoTrack_d0;
    std::vector<float> *m_recoTrack_ErrD0;
    std::vector<float> *m_recoTrack_ErrZ0;
    std::vector<float> *m_recoTrack_ErrD0Z0;

   //from the 2nd input root file
    std::vector<float> *m_truthVtx_x;
    std::vector<float> *m_truthVtx_y;
    std::vector<float> *m_truthVtx_z;
    std::vector<float> *m_recoVtx_x;
    std::vector<float> *m_recoVtx_y;
    std::vector<float> *m_recoVtx_z;   

   //for the histogram
   std::vector<float> *m_kernelA_zdata;      


   mutable std::vector<std::pair<double, int>> sortedTracks; //mutable allows you to modify a member variable within the 'const' execute function.
   mutable std::vector<std::pair<double, int>> filteredTracks;
   mutable std::vector<KDEData> accumulatedData;

   int bins = 12000;
   double z_min = -240.0;
   double z_max = 240.0;
    
};

} // namespace ActsExamples