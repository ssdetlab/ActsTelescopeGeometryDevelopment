// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Plugins/Json/TrackParametersJsonConverter.hpp"
#include "ActsExamples/TrackFinding/ITrackLookupGridReader.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"

#include <fstream>
#include <memory>

#include <nlohmann/json.hpp>

namespace ActsExamples {

namespace Experimental {

class JsonTrackLookupGridReader final : public ITrackLookupGridReader {
 public:
    using AxisGen = Acts::GridAxisGenerators::EqOpenEqOpen;

  struct Config {
    /// Binning of the grid
    std::pair<std::size_t, std::size_t> bins;
    /// Grid bounds
    std::pair<Acts::ActsScalar, Acts::ActsScalar> xBounds;
    std::pair<Acts::ActsScalar, Acts::ActsScalar> yBounds;
  };

  JsonTrackLookupGridReader(const Config& config) : m_cfg(config){};

  ~JsonTrackLookupGridReader() override = default;

  TrackLookupGrid readLookup(const std::string& path) final {
    std::cout << "READING LOOKUP" << std::endl;

    // Read the json file
    std::ifstream ifj(path);
    nlohmann::json jLookup;
    ifj >> jLookup;

    std::cout << "JLOOKUP: " << jLookup.dump(4) << std::endl;

    auto axisGen = Acts::GridAxisGenerators::EqOpenEqOpen{
        {m_cfg.xBounds.first, m_cfg.xBounds.second}, m_cfg.bins.first,
        {m_cfg.yBounds.first, m_cfg.yBounds.second}, m_cfg.bins.second};

    // Convert the grid to json
    return Acts::GridJsonConverter::fromJson<
        AxisGen,
        std::unique_ptr<Acts::CurvilinearTrackParameters>>(
            jLookup, axisGen);
  };

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The config of the writer
  Config m_cfg;
};

}  // namespace Experimental

}  // namespace ActsExamples
