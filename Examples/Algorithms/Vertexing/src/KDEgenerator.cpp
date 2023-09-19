#include "ActsExamples/Vertexing/KDEgenerator.hpp"
#include "VertexingHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>

namespace ActsExamples {
KDEAlgorithm::KDEAlgorithm(const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("KDEAlgorithm", lvl), m_cfg(cfg){

    // Open the ROOT file and access the TTree and set the TBranches

    const char* rootFilePath = "/eos/user/l/lalsaray/odd_output/tracksummary_ambi.root";
    inputFile = new TFile(rootFilePath);
    if (!inputFile || inputFile->IsZombie()) {
        ACTS_ERROR("Error opening ROOT file");
    }
    
    inputTree = dynamic_cast<TTree*>(inputFile->Get("tracksummary"));
    if (!inputTree) {
        ACTS_ERROR("Error getting TTree from ROOT file");
        inputFile->Close();
    }

    // call on the inputTree and enable automatic creation of classes to store structured data 
    inputTree->SetMakeClass(1);


    //Set Objects Pointer and set branches addresses
    bandwidth = 1.0;
    nbins = 60;
    z_min = -180.0;
    z_max = 180.0;

    d_0 = 0;
    z_0 = 0;
    sigma_d0 = 0;
    sigma_z0 = 0;
    sigma_d0_z0 = 0;

    inputTree->SetBranchAddress("eLOC0_fit", &d_0);
    inputTree->SetBranchAddress("eLOC1_fit", &z_0);
    inputTree->SetBranchAddress("err_eLOC0_fit", &sigma_d0);
    inputTree->SetBranchAddress("err_eLOC1_fit", &sigma_z0);
    inputTree->SetBranchAddress("err_eLOC0LOC1_fit", &sigma_d0_z0);

    // Create a histogram to store the KDE results
    kdeHistogram = new TH1F("kdeHistogram", "Kernel Density Estimation", 100, z_min, z_max);

    outFile = new TFile("/eos/user/l/lalsaray/KDE_output/KDE_output_file.root", "RECREATE");
          
}


ProcessCode KDEAlgorithm::execute(const AlgorithmContext&) const{

    eventNumber++;

    ACTS_INFO("Entering Execute for event " << eventNumber);

    // Load only for the current event
    ientry = inputTree->LoadTree(eventNumber);
    if (ientry < 0) return ActsExamples::ProcessCode::ABORT;
    inputTree->GetEntry(ientry);

    // // to find the minimum and maximum values of z_0
    // double min_z0 = *std::min_element(z_0->begin(), z_0->end());
    // double max_z_0 = *std::max_element(z_0->begin(), z_0->end());
    // std::cout << "Minimum value of z_0: " << min_z0 << std::endl;
    // std::cout << "Maximum value of z_0: " << max_z_0 << std::endl;

    // Clear previous event data
    sortedTracks.clear();
    filteredTracks.clear();

    // Process the current event
    for (size_t i = 0; i < z_0->size(); ++i) {
        double value = (*z_0)[i] - 3.0 * (*sigma_z0)[i];
        sortedTracks.push_back(value);
        // std::cout << "size of value (z0-3dz0) = " << value << std::endl;
    }

    // Sort the tracks in increasing order of "z_0 - 3 * sigma_z0"
    std::sort(sortedTracks.begin(), sortedTracks.end());

    // Find the minimum and maximum value of z_0
    double min_z = *std::min_element(z_0->begin(), z_0->end());
    double max_z = *std::max_element(z_0->begin(), z_0->end());

    // Iterate through sortedTracks and add tracks that meet the condition of 3sigmas window
    for (size_t i = 0; i < sortedTracks.size(); ++i) {
        double current_z_0 = (*z_0)[i];
        double current_sigma_z0 = (*sigma_z0)[i];
        
        if ((current_z_0 - 3 * current_sigma_z0) < max_z && (current_z_0 + 3 * current_sigma_z0) > min_z) {
            filteredTracks.push_back(sortedTracks[i]);
        }
    }

    // Perform grid search
    double step = (z_max - z_min) / nbins;
    double bestKDEValue = -1.0;
    double bestZ0 = 0.0;

    for (int i = 0; i < nbins; ++i) {
        double z0_candidate = z_min + i * step;
        double kdeValue = 0.0;

        // Calculate KDE value for this z0_candidate
        for (size_t j = 0; j < filteredTracks.size(); ++j) {
            double dataPoint = filteredTracks[j];
            double diff = (z0_candidate - dataPoint) / bandwidth;
            double kernel = exp(-0.5 * diff * diff) / (sqrt(2 * M_PI) * bandwidth);
            kdeValue += kernel;
        }

        // Update bestKDEValue and bestZ0 if this kdeValue is greater
        if (kdeValue > bestKDEValue) {
            bestKDEValue = kdeValue;
            bestZ0 = z0_candidate;
        }

        kdeHistogram->Fill(z0_candidate, kdeValue);
    }

    ACTS_INFO("Best KDE Value = " << bestKDEValue << " at z_0 = " << bestZ0);

    kdeHistogram->Write();

    return ActsExamples::ProcessCode::SUCCESS;

}


KDEAlgorithm::~KDEAlgorithm() {
    
    //we will perfom memory clean-up here in the destructor after all events are run

    //histogram
    if (kdeHistogram) {
        delete kdeHistogram;
        kdeHistogram = nullptr;
    }

    //output file
    if (outFile) {
        outFile->Close();
        delete outFile;
        outFile = nullptr;
    }

    //input file
    if (inputFile) {
        inputFile->Close();
        delete inputFile;
        inputFile = nullptr;
    }

    ACTS_INFO("Memory Clean-up Performed in Destructor");
}


} //namespace ActsExamples 