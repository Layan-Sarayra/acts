#include "ActsExamples/Vertexing/KDEgenerator.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>

namespace ActsExamples {


// x: the point in the 4D space (d_0, z_0, phi, theta) where we want to evaluate the Gaussian PDF. 
// and the covariance matrix describes how the variables are spread out and correlated. 
double GetGaussianPDF(const Eigen::VectorXd& x, const Eigen::VectorXd& mean, const Eigen::MatrixXd& covMat) {
    Eigen::MatrixXd covMat_inverse = covMat.inverse();
    double covMat_det = covMat.determinant();
    Eigen::VectorXd diff = x - mean;
    double chisq = (diff.transpose() * covMat_inverse * diff)(0, 0);
    return std::exp(-0.5 * chisq) / (2 * M_PI * std::sqrt(covMat_det));
}


KDEAlgorithm::KDEAlgorithm(const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("KDEAlgorithm", lvl), m_cfg(cfg) {

    ACTS_INFO("********** IN CONSTRUCTOR **********");


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

    eventNumber = -1;

    d_0 = 0;
    z_0 = 0;
    phi = 0;
    theta = 0;
    sigma_d0 = 0;
    sigma_z0 = 0;
    sigma_phi = 0;
    sigma_theta = 0;
    sigma_d0_z0 = 0;
    sigma_d0_phi = 0;
    sigma_d0_theta = 0;
    sigma_z0_phi = 0;
    sigma_z0_theta = 0;
    sigma_theta_phi = 0;


    inputTree->SetBranchAddress("eLOC0_fit", &d_0);
    inputTree->SetBranchAddress("eLOC1_fit", &z_0);
    inputTree->SetBranchAddress("ePHI_fit", &phi);
    inputTree->SetBranchAddress("eTHETA_fit", &theta);    

    inputTree->SetBranchAddress("err_eLOC0_fit", &sigma_d0);
    inputTree->SetBranchAddress("err_eLOC1_fit", &sigma_z0);
    inputTree->SetBranchAddress("err_ePHI_fit", &sigma_phi);
    inputTree->SetBranchAddress("err_eTHETA_fit", &sigma_theta);    

    inputTree->SetBranchAddress("err_eLOC0LOC1_fit", &sigma_d0_z0);
    inputTree->SetBranchAddress("err_eLOC0PHI_fit", &sigma_d0_phi);
    inputTree->SetBranchAddress("err_eLOC0THETA_fit", &sigma_d0_theta);
    inputTree->SetBranchAddress("err_eLOC1PHI_fit", &sigma_z0_phi);
    inputTree->SetBranchAddress("err_eLOC1THETA_fit", &sigma_z0_theta);
    inputTree->SetBranchAddress("err_eTHETAPHI_fit", &sigma_theta_phi);



    outFile = new TFile("/eos/user/l/lalsaray/KDE_output/KDE_output_file.root", "RECREATE");
          
}


ProcessCode KDEAlgorithm::execute(const AlgorithmContext&) const {

    ACTS_INFO("********** IN EXECUTE **********");

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
        return ((i + 0.5) / nbins) * (z_max - z_min) + z_min;};

    auto bin_min = [this](int const& i) -> double {
        return ((i + 0.0) / nbins) * (z_max - z_min) + z_min;};

    auto bin_max = [this](int const& i) -> double {
        return ((i + 1.0) / nbins) * (z_max - z_min) + z_min;};


    // Perform grid search
    double bestKDEValue = -1.0;
    double bestZ0 = 0.0;

    Eigen::MatrixXd covMat(2, 2);  // 2x2 covariance matrix for d_0 and z_0
    Eigen::VectorXd x(2);  // 2D vector for d_0 and z_0
    Eigen::VectorXd mean(2);  // 2D vector for the mean of d_0 and z_0

    mean << 0.0, 0.0;

    // Calculate the mean for d_0 and z_0
    for (size_t j = 0; j < sortedTracks.size(); ++j) {
        mean(0) += (*d_0)[j];
        mean(1) += (*z_0)[j];}

    //make sure we're not dividing by zero
    if (sortedTracks.size() > 0) {
        mean /= static_cast<double>(sortedTracks.size());}


    for (int i = 0; i < nbins; ++i) {
        double z0_candidate = bin_center(i);
        double local_z_min = bin_min(i);
        double local_z_max = bin_max(i);  
        double kdeValue = 0.0;

        std::cout << "Processing bin " << i << " with z0_candidate: " << z0_candidate << std::endl;

        // Filter tracks for this bin
        filteredTracks.clear();
        for (size_t j = 0; j < sortedTracks.size(); ++j) {
            double current_z_0 = (*z_0)[j];
            double current_sigma_z0 = (*sigma_z0)[j];
        
            if ((current_z_0 - 3 * current_sigma_z0) <= local_z_max && (current_z_0 + 3 * current_sigma_z0) >= local_z_min) {
                filteredTracks.push_back(sortedTracks[j]);
            }
        }

        // std::cout << "Number of filtered tracks: " << filteredTracks.size() << std::endl;

        // Calculate KDE value for this z0_candidate
        for (size_t j = 0; j < filteredTracks.size(); ++j) {
            // Construct 2x2 covariance matrix using d_0 and z_0
            covMat << (*sigma_d0)[j] * (*sigma_d0)[j], (*sigma_d0_z0)[j],
                     (*sigma_d0_z0)[j], (*sigma_z0)[j] * (*sigma_z0)[j];

            x << std::abs((*d_0)[j]), (*z_0)[j];

            // Calculate PDF value
            double pdfValue = GetGaussianPDF(x, mean, covMat);
            kdeValue += pdfValue;
        }

        // std::cout << "kdeValue for bin " << i << ": " << kdeValue << std::endl;

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

    ACTS_INFO("********** IN FINALIZE **********");

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

    return ActsExamples::ProcessCode::SUCCESS;
}

}  // namespace ActsExamples