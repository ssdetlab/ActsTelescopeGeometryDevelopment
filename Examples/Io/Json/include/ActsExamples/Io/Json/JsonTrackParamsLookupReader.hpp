// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Plugins/Json/TrackParametersJsonConverter.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "ActsExamples/TrackFinding/ITrackParamsLookupReader.hpp"

#include <fstream>
#include <memory>

#include <nlohmann/json.hpp>

namespace ActsExamples {

/// @brief Json reader for track parameter lookup tables
///
/// This reader is used to read track parameter lookup tables
/// from a json file to be later used in track parameter estimation
/// for seeding
class JsonTrackParamsLookupReader final : public ITrackParamsLookupReader {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Reference tracking layers
    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
        refLayers;
    /// Binning of the grid to be emposed
    /// onto the reference layers
    std::pair<std::size_t, std::size_t> bins;
  };

  JsonTrackParamsLookupReader(const Config& config) : m_cfg(config){};

  ~JsonTrackParamsLookupReader() override = default;

  Lookup readLookup(const std::string& path) final {
    // Read the json file
    std::ifstream ifj(path);
    nlohmann::json jLookup;
    ifj >> jLookup;

    Lookup lookup;
    // Iterate over the json and deserialize the grids
    for (const auto& jGrid : jLookup) {
      Acts::GeometryIdentifier id(jGrid["geo_id"]);

      if (m_cfg.refLayers.find(id) == m_cfg.refLayers.end()) {
        throw std::invalid_argument("Geometry identifier not found");
      }

      const auto* refSurface = m_cfg.refLayers.at(id);

      // Get bounds to construct the lookup grid
      auto bounds =
          dynamic_cast<const Acts::RectangleBounds*>(&refSurface->bounds());

      if (bounds == nullptr) {
        throw std::invalid_argument("Only rectangle bounds supported");
      }

      // Axis is not deserilizable, so we need to recreate it
      auto halfX = bounds->halfLengthX();
      auto halfY = bounds->halfLengthY();

      LookupAxisGen axisGen{{-halfX, halfX},
                            m_cfg.bins.first,
                            {-halfY, halfY},
                            m_cfg.bins.second};

      // Deserialize the grid
      LookupGrid grid =
          Acts::GridJsonConverter::fromJson<LookupAxisGen, LookupPair>(
              jGrid["grid"], axisGen);

      lookup.insert({id, grid});
    }

    return lookup;
  };

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The config of the writer
  Config m_cfg;
};

}  // namespace ActsExamples
