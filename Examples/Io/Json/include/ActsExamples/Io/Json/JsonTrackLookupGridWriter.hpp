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
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/EnumBitwiseOperators.hpp"
#include "ActsExamples/TrackFinding/ITrackLookupGridWriter.hpp"

#include <fstream>
#include <memory>

#include <nlohmann/json.hpp>

namespace ActsExamples {

namespace Experimental {

enum class JsonFormat : std::uint8_t {
  NoOutput = 0,
  Json = 1,
  Cbor = 2,
  All = std::numeric_limits<std::uint8_t>::max()
};

ACTS_DEFINE_ENUM_BITWISE_OPERATORS(JsonFormat)

class JsonTrackLookupGridWriter final : public ITrackLookupGridWriter {
 public:
  struct Config {
    /// Output file name
    std::string path;
    /// Output format of the file
    JsonFormat writeFormat = JsonFormat::Json;
  };

  /// Constructor
  ///
  /// @param config The configuration struct of the writer
  /// @param level The log level
  JsonTrackLookupGridWriter(const Config& config) : m_cfg(config){};

  /// Virtual destructor
  ~JsonTrackLookupGridWriter() override = default;

  /// Write out the material map
  ///
  /// @param detMaterial is the SurfaceMaterial and VolumeMaterial maps
  void writeLookup(
    const TrackLookupGrid& ipGrid,
    const TrackLookupGrid& refGrid) final {
    // Convert the grid to json
    nlohmann::json jLookupIp = Acts::GridJsonConverter::toJson(ipGrid);
    nlohmann::json jLookupRef = Acts::GridJsonConverter::toJson(refGrid);

    // Write the json file
    if (ACTS_CHECK_BIT(m_cfg.writeFormat, JsonFormat::Json)) {
      std::ofstream ofjIp(m_cfg.path + "-IP" + ".json", std::ios::out);
      ofjIp << std::setw(4) << jLookupIp << std::endl;
    std::ofstream ofjRef(m_cfg.path + "-REF" + ".json");
        ofjRef << std::setw(4) << jLookupRef << std::endl;
    } else if (ACTS_CHECK_BIT(m_cfg.writeFormat, JsonFormat::Cbor)) {
      std::vector<std::uint8_t> cborOutIp = nlohmann::json::to_cbor(jLookupIp);
      std::ofstream ofjIp(m_cfg.path + "-IP" + ".cbor", std::ios::out | std::ios::binary);
      ofjIp.write(reinterpret_cast<char*>(cborOutIp.data()), cborOutIp.size());
        std::vector<std::uint8_t> cborOutRef = nlohmann::json::to_cbor(jLookupRef);
        std::ofstream ofjRef(m_cfg.path + "-REF" + ".cbor", std::ios::out | std::ios::binary);
    }
  };

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The config of the writer
  Config m_cfg;
};

}  // namespace Experimental

}  // namespace ActsExamples
