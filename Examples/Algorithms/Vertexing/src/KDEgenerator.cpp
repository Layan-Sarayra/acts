#include "ActsExamples/Vertexing/KDEgenerator.hpp"
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
    : IAlgorithm("KDEAlgorithm", lvl), m_cfg(cfg) {
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

    eventNumber = 0;

    d_0 = 0;
    z_0 = 0;
    sigma_d0 = 0;
    sigma_z0 = 0;
    sigma_d0_z0 = 0;

    z_min = std::numeric_limits<double>::max();
    z_max = std::numeric_limits<double>::min();

    inputTree->SetBranchAddress("eLOC0_fit", &d_0);
    inputTree->SetBranchAddress("eLOC1_fit", &z_0);
    inputTree->SetBranchAddress("err_eLOC0_fit", &sigma_d0);
    inputTree->SetBranchAddress("err_eLOC1_fit", &sigma_z0);
    inputTree->SetBranchAddress("err_eLOC0LOC1_fit", &sigma_d0_z0);

    outFile = new TFile("/eos/user/l/lalsaray/KDE_output/KDE_output_file.root", "RECREATE");
          
}


ProcessCode KDEAlgorithm::execute(const AlgorithmContext&) const {

    eventNumber++;

    // Load only for the current event
    ientry = inputTree->LoadTree(eventNumber);
    if (ientry < 0) return ActsExamples::ProcessCode::ABORT;
    inputTree->GetEntry(ientry);

    // Clear previous event data
    sortedTracks.clear();
    filteredTracks.clear();

    // Process the current event
    for (size_t i = 0; i < z_0->size(); ++i) {
        double value = (*z_0)[i] - 3.0 * (*sigma_z0)[i];
        sortedTracks.push_back(value);
    }

    // Sort the tracks
    std::sort(sortedTracks.begin(), sortedTracks.end());

    // Find the minimum and maximum value of z_0
    double event_z_min = *std::min_element(z_0->begin(), z_0->end());
    double event_z_max = *std::max_element(z_0->begin(), z_0->end());

    // Update the overall z_min and z_max
    if (event_z_min < z_min) {z_min = event_z_min;}
    if (event_z_max > z_max) {z_max = event_z_max;}

    // Filter tracks
    for (size_t i = 0; i < sortedTracks.size(); ++i) {
        double current_z_0 = (*z_0)[i];
        double current_sigma_z0 = (*sigma_z0)[i];
        
        if ((current_z_0 - 3 * current_sigma_z0) < z_max && (current_z_0 + 3 * current_sigma_z0) > z_min) {
            filteredTracks.push_back(sortedTracks[i]);
        }
    }

    // Lambda functions for binning
    auto bin_center = [this](int const& i) -> double {
        return ((i + 0.5) / this->nbins) * (this->z_max - this->z_min) + this->z_min;
    };

    auto bin_min = [this](int const& i) -> double {
        return ((i + 0.0) / this->nbins) * (this->z_max - this->z_min) + this->z_min;
    };

    auto bin_max = [this](int const& i) -> double {
        return ((i + 1.0) / this->nbins) * (this->z_max - this->z_min) + this->z_min;
    };


    // Perform grid search
    double bestKDEValue = -1.0;
    double bestZ0 = 0.0;

    for (int i = 0; i < nbins; ++i) {
        double z0_candidate = bin_center(i);
        z_min = bin_min(i);
        z_max = bin_max(i);        
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

        // Store the KDE data for this event
        KDEData dataPoint;
        dataPoint.z0_candidate = z0_candidate;
        dataPoint.kdeValue = kdeValue;
        accumulatedData.push_back(dataPoint);
    }

    ACTS_INFO("Best KDE Value = " << bestKDEValue << " at z_0 = " << bestZ0);

    return ActsExamples::ProcessCode::SUCCESS;

}

ProcessCode KDEAlgorithm::finalize() {
    // Create a histogram to store the KDE results
    kdeHistogram = new TH1F("kdeHistogram", "Kernel Density Estimation", 60, z_min, z_max);

    // Fill the histogram using the accumulated data
    for (const auto& dataPoint : accumulatedData) {
        kdeHistogram->Fill(dataPoint.z0_candidate, dataPoint.kdeValue);
    }

    // Write the histogram to file
    kdeHistogram->Write();

    //clean-up

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

    ACTS_INFO("Memory Clean-up Performed While Finalizing");

    return ActsExamples::ProcessCode::SUCCESS;
}

}  // namespace ActsExamples