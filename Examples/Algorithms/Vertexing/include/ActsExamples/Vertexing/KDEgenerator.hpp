#pragma once /* to ensure that this file is only included -read- once even if code executed itirates in a single compilation */

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <vector>
#include <memory>
#include <string>


//zmin and zmax corresponding to +-24cm (with 12000 bins, each bin is 40 micrometers)
constexpr unsigned int kde_points = 12000, n_features = 4;
constexpr float zmin = -240.f, zmax = 240.f; //in mm
constexpr float binWidth = 0.04f; //in mm

struct Kernel {
  std::array<float, kde_points> zdata{};
  std::array<float, kde_points> xmax{};
  std::array<float, kde_points> ymax{};
  void clear() {
    zdata.fill(0.f);
    xmax.fill(0.f);
    ymax.fill(0.f);
  }
};


namespace ActsExamples {

  class KDEAlgorithm final : public IAlgorithm {
  public:
    struct Config { /// defines a nested structure within the class to group related configuration options
      /// The message we are going to log.
      std::string message;
    };
    ///logging: record information and events during execution by sending (info, warning, error, etc) to console -output- 
    KDEAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

    /// This function will be called on each event by the sequencer. (coded in IAlgorithm -> AlgorithmContext)
    ProcessCode execute(const AlgorithmContext& ctx) const final;

    std::array<Kernel,2> const getKDE_from_PDF() const;


    const Config& config() const { return m_cfg; }

  private:
    Config m_cfg;
  };

}
