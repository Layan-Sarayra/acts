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


// x: the point in the 2D space (d_0, z_0) where we want to evaluate the Gaussian PDF. 
// and the covariance matrix describes how the variables are spread out and correlated. 
double GetGaussianPDF(const Eigen::Vector2d& x, const Eigen::Vector2d& mean, const Eigen::MatrixXd& covMat) {
    Eigen::MatrixXd covMat_inverse = covMat.inverse();
    double covMat_det = covMat.determinant();

    //debugging lines... validate the covariance matrix elements
    if (covMat_det < 0) {
        std::cerr << "Determinant of covariance matrix is negative" << std::endl;
        return 0;
    }
    if (covMat_det == 0) {
        std::cerr << "Covariance matrix is singular" << std::endl;
        return 0;
    }

    Eigen::Vector2d diff = x - mean;
    double chisq = (diff.transpose() * covMat_inverse * diff)(0, 0);
    return std::exp(-0.5 * chisq) / (2 * M_PI * std::sqrt(covMat_det));
}


KDEAlgorithm::KDEAlgorithm(const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("KDEAlgorithm", lvl), m_cfg(cfg) {

    ACTS_INFO("********** IN CONSTRUCTOR **********");

    // Open the first ROOT file and access the TTree and set the TBranches
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


    // Open the second ROOT file and access the TTree and set the TBranches
    const char* performanceFilePath = "/eos/user/r/rgarg/Rocky/ACTS_Project/PVFinder/Data_Layan/odd_output/performance_vertexing.root";
    performanceFile = new TFile(performanceFilePath);
    if (!performanceFile || performanceFile->IsZombie()) {
        ACTS_ERROR("Error opening ROOT file");
    }
    
    performanceTree = dynamic_cast<TTree*>(performanceFile->Get("vertexing"));
    if (!performanceTree) {
        ACTS_ERROR("Error getting TTree from ROOT file");
        inputFile->Close();
    }

    performanceTree->SetBranchAddress("truthX", &truthX);
    performanceTree->SetBranchAddress("truthY", &truthY);   
    performanceTree->SetBranchAddress("truthZ", &truthZ);
    performanceTree->SetBranchAddress("recoX", &recoX);   
    performanceTree->SetBranchAddress("recoY", &recoY);
    performanceTree->SetBranchAddress("recoZ", &recoZ);


    // Defining output file and output tree & branches
    outFile = new TFile("/eos/user/l/lalsaray/KDE_output/KDE_output_file.root", "RECREATE");
    outFile->cd();
    outputTree = new TTree("PVFinderData", "PVFinderData");

    // write data from the first ROOT file; tracksummary_ambi.root into the output file KDE_output_file.root
    outputTree->Branch("RecoTrack_z0", &m_recoTrack_z0);
    outputTree->Branch("RecoTrack_d0", &m_recoTrack_d0);
    outputTree->Branch("RecoTrack_ErrD0", &m_recoTrack_ErrD0);
    outputTree->Branch("RecoTrack_ErrZ0", &m_recoTrack_ErrZ0);
    outputTree->Branch("RecoTrack_ErrD0Z0", &m_recoTrack_ErrD0Z0);

    // write data from the second ROOT file into the output file; performance_vertexing.root
    outputTree->Branch("TruthVertex_x", &m_truthVtx_x);    
    outputTree->Branch("TruthVertex_y", &m_truthVtx_y);
    outputTree->Branch("TruthVertex_z", &m_truthVtx_z);
    outputTree->Branch("RecoVertex_x", &m_recoVtx_x);
    outputTree->Branch("RecoVertex_y", &m_recoVtx_y);
    outputTree->Branch("RecoVertex_z", &m_recoVtx_z);

    //write the histogram into the output file
    outputTree->Branch("KernelA_zdata", &m_kernelA_zdata);

    eventNumber = -1;


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
    // // Print sorted tracks
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


    for (int i = 0; i < bins; ++i) {

        double z0_candidate = bin_center(bins, z_min, z_max, i);
        double local_z_min = bin_min(bins, z_min, z_max, i);
        double local_z_max = bin_max(bins, z_min, z_max, i);  
   
        std::vector<double> kernel_value = {-1.0, -1.0};
        Eigen::Vector3d best;         

        // std::cout << "Processing bin " << i << " with z0_candidate: " << z0_candidate << std::endl; 

        filteredTracks.clear();
        // Filter tracks for this bin
        for (size_t j = 0; j < sortedTracks.size(); ++j) {
            int index = sortedTracks[j].second;
            double current_z_0 = (*z_0)[index];
            double current_sigma_z0 = (*sigma_z0)[index];

            if ((current_z_0 - 3 * current_sigma_z0) <= local_z_max && (current_z_0 + 3 * current_sigma_z0) >= local_z_min) {
                filteredTracks.push_back(sortedTracks[j]);

                // std::cout << "filteredTracks value is = (" << filteredTracks.back().first << ", " << filteredTracks.back().second << ")" << std::endl;
                // // std::cout << "Number of filtered tracks: " << filteredTracks.size() << std::endl;
                // std::cout << "Current z_0: " << current_z_0 << ", Current sigma_z0: " << current_sigma_z0 << std::endl;
                // std::cout << "Local z_min: " << local_z_min << ", Local z_max: " << local_z_max << std::endl; 
            }
            // add a break statement to make code more efficient
            if ((current_z_0 - 3 * current_sigma_z0) > local_z_max) {
                break;  // Add break statement
            }
        } 


        KDEData dataPoint;
        dataPoint.z0_candidate = z0_candidate;
        dataPoint.kdeValue = -1;

        // here we define a lambda variable to evaluate the pdf, set kernel_value maximum and it's position (just defined here, used afterwards) 
        auto eval_pdf_max_and_position = [&kernel_value, &best, &dataPoint, this, &bin_center, &z0_candidate](Eigen::Vector3d& p) {
            // p is the xyz point in space

            double this_kernel = 0.0, this_kernel_sq = 0.0;
            Eigen::VectorXd x(2);  // 2D vector for d_0 and z_0

            // Iterate over the filtered tracks
            for (size_t k = 0; k < filteredTracks.size(); ++k) {
                int index = filteredTracks[k].second;
                // Construct 2x2 covariance matrix using d_0 and z_0
                Eigen::MatrixXd covMat(2, 2);
                covMat << (*sigma_d0)[index] , (*sigma_d0_z0)[index],
                          (*sigma_d0_z0)[index], (*sigma_z0)[index];

                // std::cout << "the covariance matrix = " << covMat << std::endl;                 
                // std::cout << "Determinant of covMat: " << covMat.determinant() << std::endl;


                // Calculate the PDF for this track and point p, and accumulate the kernel value
                double pdf_for_this_track = GetGaussianPDF(Eigen::Vector2d(std::sqrt(p.x() * p.x() + p.y() * p.y()), p.z()), Eigen::Vector2d((*d_0)[index], (*z_0)[index]), covMat);
                this_kernel += pdf_for_this_track;
                this_kernel_sq += pdf_for_this_track * pdf_for_this_track;

                // std::cout << "GetGaussianPDF return: " << pdf_for_this_track << std::endl;

            }


            // Check and update the best kernel values and positions
            if (this_kernel > kernel_value[0]) {
                kernel_value[0] = this_kernel;
                kernel_value[1] = this_kernel_sq;
                best = p;
            }


            dataPoint.kdeValue = kernel_value[0];

        };

        int nbxy = 60;
        float xymin = -0.6f, xymax = 0.6f;

        for (int bx = 0; bx < nbxy; bx++) {
            for (int by = 0; by < nbxy; by++) {
                Eigen::Vector3d p(bin_center(nbxy, xymin, xymax, bx), bin_center(nbxy, xymin, xymax, by), z0_candidate);
                eval_pdf_max_and_position(p);
            }
        }


        // add the data to the accumulatedData vector
        accumulatedData.push_back(dataPoint);


        ACTS_INFO("KDE Value = " << kernel_value[0] << " z0 = " << z0_candidate);

    }

    return ActsExamples::ProcessCode::SUCCESS;
}


void KDEAlgorithm::copyBranches() {
    ACTS_INFO("********** COPYING BRANCHES **********");

    // Loop over the entries and fill the output tree
    Long64_t nentries = inputTree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        inputTree->GetEntry(i);
        performanceTree->GetEntry(i);
        outputTree->Fill();
    }

    // Write the output tree and close the files
    outFile->cd();
    outputTree->Write();

}


ProcessCode KDEAlgorithm::finalize() {

    ACTS_INFO("********** IN FINALIZE **********");

    outFile->cd();  // Make sure we're in the correct directory in the output file

    outputTree->Write();

    // Create a histogram to store the KDE results
    kdeHistogram = new TH1F("kdeHistogram", "Kernel Density Estimation", bins, z_min, z_max);

    // Fill the histogram using the accumulated data
    for (size_t i = 0; i < accumulatedData.size(); ++i) {
        kdeHistogram->Fill(accumulatedData[i].z0_candidate, accumulatedData[i].kdeValue);
    }

    // // Extract data from the histogram and store it in m_kernelA_zdata
    // m_kernelA_zdata->clear();  // Clear any previous data
    // for (int i = 1; i <= kdeHistogram->GetNbinsX(); ++i) {
    //     m_kernelA_zdata->push_back(kdeHistogram->GetBinContent(i));
    // }
    // //Fill Again After extracting data from the histogram and before writing the histogram
    // outputTree->Fill();

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

    //input file; tracksummary_ambi.root
    if (inputFile) {
        inputFile->Close();
        delete inputFile;
        inputFile = nullptr;
    }


    //second input file; performance_vertexing.root
    if (performanceFile) {
        performanceFile->Close();
        delete performanceFile;
        performanceFile = nullptr;
    }

    return ActsExamples::ProcessCode::SUCCESS;
}

}  // namespace ActsExamples