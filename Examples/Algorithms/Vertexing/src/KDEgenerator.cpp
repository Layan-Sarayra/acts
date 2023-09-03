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
    eLOC1_fit = 0;
    err_eLOC0_fit = 0;
    err_eLOC1_fit = 0;
    err_eLOC0LOC1_fit = 0;
    inputTree->SetBranchAddress("eLOC0_fit", &eLOC0_fit);
    inputTree->SetBranchAddress("eLOC1_fit", &eLOC1_fit);
    inputTree->SetBranchAddress("err_eLOC0_fit", &err_eLOC0_fit);
    inputTree->SetBranchAddress("err_eLOC1_fit", &err_eLOC1_fit);
    inputTree->SetBranchAddress("err_eLOC0LOC1_fit", &err_eLOC0LOC1_fit);

    outFile = new TFile("/eos/user/l/lalsaray/KDE_output/KDE_output_file.root", "RECREATE");

  }


ProcessCode KDEAlgorithm::execute(const AlgorithmContext&) const {

    Long64_t nentries = inputTree->GetEntries();

    // Define KDE parameters
    double bandwidth = 1.0; // You can adjust this value as needed

    // Create a histogram to store the KDE results
    TH1F* kdeHistogram = new TH1F("kdeHistogram", "Kernel Density Estimation", 100, -7.0, 7.0);
    kdeHistogram->SetMaximum(5000);
    outFile->cd(); //switches to write data to the latest directory; aka outFile instead of inputFile

    for (Long64_t entry = 0; entry < nentries; entry++) {
        // Must load the tree first to get the data; otherwise, 0 values would be printed
        Long64_t ientry = inputTree->LoadTree(entry);
        if (ientry < 0) break;
        inputTree->GetEntry(ientry);

        // Add the values of "eLOC0_fit - 3 * eLOC1_fit" to the sortedTracks vector defined in the .hpp file
        for (size_t i = 0; i < eLOC0_fit->size(); ++i) {
            double value = (*eLOC0_fit)[i] - 3.0 * (*eLOC1_fit)[i];
            sortedTracks.push_back(value);
        }
    }

    // Sort the tracks in increasing order of "eLOC0_fit - 3 * eLOC1_fit"
    std::sort(sortedTracks.begin(), sortedTracks.end());

    // Find the minimum value of eLOC0_fit
    double minELOC0_fit = *std::min_element(eLOC0_fit->begin(), eLOC0_fit->end());

    // Iterate through sortedTracks and add tracks that meet the condition to the filteredTracks vector defined in the .hpp file
    for (size_t i = 0; i < sortedTracks.size(); ++i) {
        if (sortedTracks[i] >= minELOC0_fit - 3.0 * (*err_eLOC0_fit)[i]) {
            filteredTracks.push_back(sortedTracks[i]);
        }
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
    }

    kdeHistogram->Write();

    // Clean up memory by closing the ROOT file
    delete kdeHistogram;
    outFile->Close();

    // Close and clear up the input file memory
    if (inputFile) {
        inputFile->Close();
        delete inputFile;
    }

    return ActsExamples::ProcessCode::SUCCESS;
}
 
} //namespace ActsExamples 