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
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/EnumBitwiseOperators.hpp"
#include "ActsExamples/TrackFinding/ITrackParamsLookupWriter.hpp"

#include <fstream>
#include <memory>

#include <nlohmann/json.hpp>

namespace ActsExamples {

class JsonTrackParamsLookupWriter final : public ITrackParamsLookupWriter {
    public:
        struct Config {
            /// Output file name
            std::string path;
        };

        /// Constructor
        ///
        /// @param config The configuration struct of the writer
        /// @param level The log level
        JsonTrackParamsLookupWriter(const Config& config) : m_cfg(config){};
        
        /// Virtual destructor
        ~JsonTrackParamsLookupWriter() override = default;
        
        /// Write out the material map
        ///
        /// @param detMaterial is the SurfaceMaterial and VolumeMaterial maps
        void writeLookup(const Lookup& lookup) final {
            nlohmann::json jLookup;
            for (const auto& [id, grid] : lookup) {
                nlohmann::json jGrid;
                jGrid["geo_id"] = id.value();
                jGrid["grid"] = Acts::GridJsonConverter::toJson(grid);
    
                jLookup.push_back(jGrid);
            }
    
            // Write the json file
            std::ofstream ofj(m_cfg.path + ".json", std::ios::out);
            ofj << std::setw(4) << jLookup << std::endl;
        };
        
        /// Readonly access to the config
        const Config& config() const { return m_cfg; }
        
        private:
        /// The config of the writer
        Config m_cfg;
};

}  // namespace ActsExamples
