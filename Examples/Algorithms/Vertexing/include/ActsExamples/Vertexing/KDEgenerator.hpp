#pragma once /* to ensure that this file is only included -read- once even if code executed itirates in a single compilation */

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <vector>
#include <memory>
#include <string>

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

    const Config& config() const { return m_cfg; }

  private:
    Config m_cfg;
  };

}
