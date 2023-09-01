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

    outFile = new TFile("/eos/user/l/lalsaray/KDE_output/KDE_output_file.root", "RECREATE");

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

  }


ProcessCode KDEAlgorithm::execute(const AlgorithmContext&) const {

    Long64_t nentries = inputTree->GetEntries();
    
    // Define KDE parameters
    double bandwidth = 1.0; // You can adjust this value as needed
    
    // Create a histogram to store the KDE results
    TH1F* kdeHistogram = new TH1F("kdeHistogram", "Kernel Density Estimation", 100, -150.0, 350.0);

    //std::cout << "Number of events = " << nentries << std::endl; // this should give you 100 events.

    for(Long64_t entry = 0; entry<nentries; entry++){
        //must load the tree first to get the data; otherwise 0 values would be printed
        Long64_t ientry = inputTree->LoadTree(entry);
        if(ientry < 0) break;
        //then we get the entries (events stored)
        inputTree->GetEntry(ientry);

        for (size_t i = 0; i < eLOC0_fit->size(); ++i) {
            double dataPoint1 = (*eLOC0_fit)[i];
            double kdeValue = 0.0;

            // for (size_t j = 0; j < eLOC0_fit->size(); ++j) {
            //     double dataPoint2 = (*eLOC0_fit)[j];
            //     double diff = (dataPoint1 - dataPoint2) / bandwidth;
            //     double kernel = exp(-0.5 * diff * diff) / (sqrt(2 * M_PI) * bandwidth);
            //     kdeValue += kernel;
            // }
            kdeHistogram->Fill(dataPoint1, kdeValue);
        }

        //std::cout << "size of eLOC0_fit = " << eLOC0_fit->size() << std::endl;
    
    }


    kdeHistogram->Write();

    //clean up memory by closing the ROOT file
    delete kdeHistogram;
    outFile->Close();
    //inputFile->Close();

    //close and clear up the input file memory
    if (inputFile) {
        inputFile->Close();
        delete inputFile;
    }

    return ActsExamples::ProcessCode::SUCCESS;
}
 
} //namespace ActsExamples 