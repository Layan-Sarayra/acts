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
    eLOC0_fit = 0;
    z_0 = 0;
    err_eLOC0_fit = 0;
    sigma_z_0 = 0;
    err_eLOC0LOC1_fit = 0;

    inputTree->SetBranchAddress("eLOC0_fit", &eLOC0_fit);
    inputTree->SetBranchAddress("eLOC1_fit", &z_0);
    inputTree->SetBranchAddress("err_eLOC0_fit", &err_eLOC0_fit);
    inputTree->SetBranchAddress("err_eLOC1_fit", &sigma_z_0);
    inputTree->SetBranchAddress("err_eLOC0LOC1_fit", &err_eLOC0LOC1_fit);

    // Define KDE parameters
    bandwidth = 1.0; // You can adjust this value as needed

    // Create a histogram to store the KDE results
    kdeHistogram = new TH1F("kdeHistogram", "Kernel Density Estimation", 100, -7.0, 7.0);
    kdeHistogram->SetMaximum(5000);

    outFile = new TFile("/eos/user/l/lalsaray/KDE_output/KDE_output_file.root", "RECREATE");
          
}


ProcessCode KDEAlgorithm::execute(const AlgorithmContext&) const{

    eventNumber++;

    std::cout << "Entering Execute for event " << eventNumber << std::endl;

    // Load only for the current event
    ientry = inputTree->LoadTree(eventNumber);
    if (ientry < 0) return ActsExamples::ProcessCode::ABORT;
    inputTree->GetEntry(ientry);

    // Clear previous event data
    sortedTracks.clear();
    filteredTracks.clear();

    // Process the current event
    for (size_t i = 0; i < z_0->size(); ++i) {
        double value = (*z_0)[i] - 3.0 * (*sigma_z_0)[i];
        sortedTracks.push_back(value);
        // std::cout << "size of value (z0-3dz0) = " << value << std::endl;
    }

    // Sort the tracks in increasing order of "z_0 - 3 * sigma_z_0"
    std::sort(sortedTracks.begin(), sortedTracks.end());

    // Find the minimum value of z_0
    double min_z_0 = *std::min_element(z_0->begin(), z_0->end());

    // Iterate through sortedTracks and add tracks that meet the condition to the filteredTracks vector defined in the .hpp file
    for (size_t i = 0; i < sortedTracks.size(); ++i) {
        if (sortedTracks[i] >= min_z_0 - 3.0 * (*sigma_z_0)[i]) {
            filteredTracks.push_back(sortedTracks[i]);
        }
        // std::cout << "size of sortedTracks =" << sortedTracks.size() << std::endl;
    }

    // Perform KDE on the filtered tracks
    for (size_t i = 0; i < filteredTracks.size(); ++i) {
        double dataPoint1 = filteredTracks[i];
        double kdeValue = 0.0;

        for (size_t j = 0; j < sortedTracks.size(); ++j) {
            double dataPoint2 = filteredTracks[j];
            double diff = (dataPoint1 - dataPoint2) / bandwidth;
            double kernel = exp(-0.5 * diff * diff) / (sqrt(2 * M_PI) * bandwidth);
            kdeValue += kernel;
        }
        kdeHistogram->Fill(dataPoint1, kdeValue);
        // std::cout << "size of filteredTracks =" << filteredTracks.size() << std::endl;

    }

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

    std::cout << "Memory Clean-up Performed in Destructor" << std::endl;
}


} //namespace ActsExamples 