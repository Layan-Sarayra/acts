#include "ActsExamples/Vertexing/KDEgenerator.hpp"
#include "VertexingHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

using namespace std;

namespace ActsExamples {

KDEAlgorithm::KDEAlgorithm(const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("KDEAlgorithm", lvl), m_cfg(cfg) {}

ProcessCode KDEAlgorithm::execute(const AlgorithmContext& ctx) const {

    // Open the ROOT file and access the TTree and TBranches
    const char* rootFilePath = "/eos/user/l/lalsaray/odd_output/tracksummary_ambi.root";

    TFile* inputFile = TFile::Open(rootFilePath);
    if (!inputFile || inputFile->IsZombie()) {
        ACTS_ERROR("Error opening ROOT file");
        return ActsExamples::ProcessCode::ABORT;
    }

    TTree* inputTree = dynamic_cast<TTree*>(inputFile->Get("tracksummary"));
    if (!inputTree) {
        ACTS_ERROR("Error getting TTree from ROOT file");
        inputFile->Close();
        return ActsExamples::ProcessCode::ABORT;
    }

    std::vector<float> err_eLOC0_fit;  // chose data/branch needed
    
    TBranch* br_err_eLOC0_fit = inputTree->GetBranch("err_eLOC0_fit");
    br_err_eLOC0_fit->SetAddress(&err_eLOC0_fit);

    // Define KDE parameters
    double bandwidth = 1.0; // You can adjust this value based on your data

    // Create a histogram to store the KDE results
    TH1F* kdeHistogram = new TH1F("kdeHistogram", "Kernel Density Estimation", 100, -5.0, 5.0);

    // Loop over the entries and perform calculations
    for (Long64_t entry = 0; entry < inputTree->GetEntries(); ++entry) {
        inputTree->GetEntry(entry);

        for (size_t i = 0; i < err_eLOC0_fit.size(); ++i) {
            double dataPoint = err_eLOC0_fit[i];

            double kdeValue = 0.0;
            for (size_t j = 0; j < err_eLOC0_fit.size(); ++j) {
                double diff = (dataPoint - err_eLOC0_fit[j]) / bandwidth;
                double kernel = exp(-0.5 * diff * diff) / (sqrt(2 * M_PI) * bandwidth);
                kdeValue += kernel;
            }
            kdeHistogram->Fill(dataPoint, kdeValue);
        }
    }

    // Now, `kdeHistogram` contains the KDE estimate of the data

    // Save the histogram
    TFile* outputFile = TFile::Open("/eos/user/l/lalsaray/KDE_output.root", "RECREATE");
    kdeHistogram->Write();
    outputFile->Close();

    //clean up memory by closing the ROOT file
    delete kdeHistogram;
    inputFile->Close();

    return ActsExamples::ProcessCode::SUCCESS;
}

} // namespace ActsExamples