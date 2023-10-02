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


// x: the point in the 4D space (d_0, z_0) where we want to evaluate the Gaussian PDF. 
// and the covariance matrix describes how the variables are spread out and correlated. 
double GetGaussianPDF(const Eigen::Vector2d& x, const Eigen::Vector2d& mean, const Eigen::MatrixXd& covMat) {
    Eigen::MatrixXd covMat_inverse = covMat.inverse();
    double covMat_det = covMat.determinant();
    Eigen::Vector2d diff = x - mean;
    double chisq = (diff.transpose() * covMat_inverse * diff)(0, 0);
    return std::exp(-0.5 * chisq) / (2 * M_PI * std::sqrt(covMat_det));
}


KDEAlgorithm::KDEAlgorithm(const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("KDEAlgorithm", lvl), m_cfg(cfg) {

    ACTS_INFO("********** IN CONSTRUCTOR **********");


    // Open the ROOT file and access the TTree and set the TBranches

    const char* rootFilePath = "/eos/user/r/rgarg/Rocky/ACTS_Project/PVFinder/Data_Layan/odd_output/tracksummary_ambi.root";
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

    // Clear data of the previous event
    sortedTracks.clear();
    filteredTracks.clear();

    // Process the current event
    for (size_t i = 0; i < z_0->size(); ++i) {
        double value = (*z_0)[i] - 3.0 * (*sigma_z0)[i];
        sortedTracks.emplace_back(value, i); //storing the value while keeping track of the index number associated with the sortedTrack
    }

    // Validating sorting tracks function is working properly
    // First:
    // Print pre-sorted tracks
    // std::cout << "pre-sorted tracks = ";
    // for (const auto& track : sortedTracks) {
    //     std::cout << "(" << track.first << ", " << track.second << ") ";
    // }
    // std::cout << std::endl;

    // Sort the tracks
    std::sort(sortedTracks.begin(), sortedTracks.end(),
        [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
            return a.first < b.first;
        }
    );

    // Second:
    // Print sorted tracks
    // std::cout << "sorted tracks = ";
    // for (const auto& track : sortedTracks) {
    //     std::cout << "(" << track.first << ", " << track.second << ") ";
    // }
    // std::cout << std::endl;


    // Define Lambda functions for binning, will be called later
    auto bin_center = [] (int const& nbins, float const& min, float const& max, int const& i) -> float {
        return ((i + 0.5f) / nbins) * (max - min) + min;
    };
    auto bin_min = [] (int const& nbins, float const& min, float const& max, int const& i) -> float {
        return ((i + 0.0f) / nbins) * (max - min) + min;
    };
    auto bin_max = [] (int const& nbins, float const& min, float const& max, int const& i) -> float {
        return ((i + 1.0f) / nbins) * (max - min) + min;
    };



    // Initiate variables used to perform grid search later
    double bestKDEValue = -1.0;
    double bestZ0 = 0.0;

    int bins = 1200;

    for (int i = 0; i < bins; ++i) {

        double z0_candidate = bin_center(bins, z_min, z_max, i);
        double local_z_min = bin_min(bins, z_min, z_max, i);
        double local_z_max = bin_max(bins, z_min, z_max, i);  
   
        std::vector<double> kernel_value = {-1.0, -1.0};
        Eigen::Vector3d best;         

        // std::cout << "Processing bin " << i << " with z0_candidate: " << z0_candidate << std::endl;

        // Filter tracks for this bin
        for (size_t j = 0; j < sortedTracks.size(); ++j) {
            int index = sortedTracks[j].second;
            double current_z_0 = (*z_0)[index];
            double current_sigma_z0 = (*sigma_z0)[index];
        
            if ((current_z_0 - 3 * current_sigma_z0) <= local_z_max && (current_z_0 + 3 * current_sigma_z0) >= local_z_min) {
                filteredTracks.push_back(sortedTracks[j]);
            }
            // add a break statement to make code more efficient
            if ((current_z_0 - 3 * current_sigma_z0) > local_z_max) {
                break;  // Add break statement
            }        
        }    

        // std::cout << "Number of filtered tracks: " << filteredTracks.size() << std::endl;

        // here we define a lambda variable to evaluate the pdf, set kernel_value maximum and it's position (just defined here, used afterwards) 
        auto eval_pdf_max_and_position = [&kernel_value, &best, this, &bin_center, &z0_candidate](Eigen::Vector3d& p) {

            double this_kernel = 0.0, this_kernel_sq = 0.0;

            Eigen::VectorXd x(2);  // 2D vector for d_0 and z_0
            // Eigen::VectorXd mean(2);  // 2D vector for the mean of d_0 and z_0

            // Iterate over the filtered tracks
            for (const auto& track : filteredTracks) {
                int index = track.second;
                // Construct 2x2 covariance matrix using d_0 and z_0
                Eigen::MatrixXd covMat(2, 2);
                covMat << (*sigma_d0)[index] * (*sigma_d0)[index], (*sigma_d0_z0)[index],
                          (*sigma_d0_z0)[index], (*sigma_z0)[index] * (*sigma_z0)[index];

                // Calculate the PDF for this track and point p, and accumulate the kernel value

                double pdf_for_this_track = GetGaussianPDF(Eigen::Vector2d(p.x(), p.z()), Eigen::Vector2d((*d_0)[index], (*z_0)[index]), covMat);
                this_kernel += pdf_for_this_track;
                this_kernel_sq += pdf_for_this_track * pdf_for_this_track;
            }

            // Check and update the best kernel values and positions
            if (this_kernel > kernel_value[0]) {
                kernel_value[0] = this_kernel;
                kernel_value[1] = this_kernel_sq;
                best = p;
            }
        };

        int nbxy = 60;
        float xymin = -0.6f, xymax = 0.6f;

        for (int bx = 0; bx < nbxy; bx++) {
            for (int by = 0; by < nbxy; by++) {
                Eigen::Vector3d p(bin_center(nbxy, xymin, xymax, bx), bin_center(nbxy, xymin, xymax, by), z0_candidate);
                eval_pdf_max_and_position(p);
            }
        }
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