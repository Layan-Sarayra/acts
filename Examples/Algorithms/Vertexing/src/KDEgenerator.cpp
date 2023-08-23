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

namespace ActsExamples {

KDEAlgorithm::KDEAlgorithm(const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("KDEAlgorithm", lvl), m_cfg(cfg) {}

ProcessCode KDEAlgorithm::execute(const AlgorithmContext& ctx) const {

    // Open the ROOT file and access the TTree and TBranches
    const char* rootFilePath = "/eos/user/l/lalsaray/odd_output/tracksummary_ambi.root";

    TFile* inputFile = TFile::Open(rootFilePath);
    std::cout << "I broke here :) 1" << std::endl;
    if (!inputFile || inputFile->IsZombie()) {
        ACTS_ERROR("Error opening ROOT file");
        return ActsExamples::ProcessCode::ABORT;
    }

    TTree* inputTree = dynamic_cast<TTree*>(inputFile->Get("tracksummary"));
    std::cout << "I broke here :) 2" << std::endl;
    if (!inputTree) {
        ACTS_ERROR("Error getting TTree from ROOT file");
        inputFile->Close();
        return ActsExamples::ProcessCode::ABORT;
    }

    // Define variables to store data from the branches of interest
    std::vector<float> err_eLOC0_fit;

    TBranch* br_err_eLOC0_fit = inputTree->GetBranch("err_eLOC0_fit");
    std::cout << "I broke here :) 3" << std::endl;
    if (!br_err_eLOC0_fit) {
        ACTS_ERROR("Branch 'err_eLOC0_fit' not found in the ROOT file");
        inputFile->Close();
        return ActsExamples::ProcessCode::ABORT;
    }
    br_err_eLOC0_fit->SetAddress(&err_eLOC0_fit);
    std::cout << "I broke here :) 4" << std::endl;

    // Define KDE parameters
    double bandwidth = 1.0;

    // text file to save the values
    std::ofstream outputFile("KDE_values.txt");
    std::cout << "I broke here :) 5" << std::endl;
    if (!outputFile.is_open()) {
        ACTS_ERROR("Error opening text file for output");
        inputFile->Close();
        return ActsExamples::ProcessCode::ABORT;
    }

    // Loop over the entries to perform calculations
    for (Long64_t entry = 0; entry < inputTree->GetEntries(); ++entry) {
        inputTree->GetEntry(entry);
        std::cout << "I broke here :) 6" << std::endl;

        // Check if err_eLOC0_fit vector is empty
        if (err_eLOC0_fit.empty()) {
            std::cout << "err_eLOC0_fit is empty for entry " << entry << ", skipping the loop." << std::endl;
            continue;  // Skip this iteration of the loop
        }

        for (size_t i = 0; i < err_eLOC0_fit.size(); ++i) {
            double dataPoint = err_eLOC0_fit[i];
            double kdeValue = 0.0;
            std::cout << "I broke here :) 7" << std::endl;

            for (size_t j = 0; j < err_eLOC0_fit.size(); ++j) {
                double diff = (dataPoint - err_eLOC0_fit[j]) / bandwidth;
                double kernel = exp(-0.5 * diff * diff) / (sqrt(2 * M_PI) * bandwidth);
                kdeValue += kernel;
            }
            std::cout << "I broke here :) 8" << std::endl;

            // Save the KDE value to the text file
            outputFile << kdeValue << "\n";
            std::cout << "I broke here :) 10" << std::endl;
        }
    }

    outputFile.close();
    inputFile->Close();
    return ActsExamples::ProcessCode::SUCCESS;
}

} // namespace ActsExamples