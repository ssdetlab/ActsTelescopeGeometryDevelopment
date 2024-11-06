// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/TrackFinding/ITrackParamsLookupReader.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

namespace ActsExamples {

class TrackParamsLookupProvider {
    public:
        using SourceLinkCalibrator = 
            Acts::Delegate<Acts::Vector2(
                const Acts::GeometryContext&, 
                const ActsExamples::IndexSourceLink&)>;

        struct Extensions {
            /// Source link calibration
            SourceLinkCalibrator sourceLinkCalibrator;
        };

        struct Config {
            /// Lookup reader
            std::shared_ptr<ITrackParamsLookupReader> trackLookupReader;
            /// Lookup path
            std::string lookupPath;
        };
        
        TrackParamsLookupProvider(const Config& config)
            : m_cfg(std::move(config)),
            m_lookup(std::make_shared<Lookup>(
                m_cfg.trackLookupReader->readLookup(m_cfg.lookupPath))) {}

        /// Lookup the track parameters at a given position
        std::pair<
            Acts::CurvilinearTrackParameters, 
            Acts::CurvilinearTrackParameters> lookup(
                const Acts::GeometryContext& gctx, 
                const Acts::SourceLink& pivot) const {
                    auto idxSl = pivot.get<ActsExamples::IndexSourceLink>();
    
                    Acts::Vector2 localPos = m_extensions.sourceLinkCalibrator(
                        gctx, idxSl);
    
                    auto bin = m_lookup->at(
                        idxSl.geometryId()).localBinsFromPosition(localPos);
    
                    return {
                        *m_lookup->at(idxSl.geometryId()).atLocalBins(bin).first,  
                        *m_lookup->at(idxSl.geometryId()).atLocalBins(bin).second};
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

        Extensions& extensions() { return m_extensions; }

    private:
        Config m_cfg;

        Extensions m_extensions;

        std::shared_ptr<Lookup> m_lookup;
};

}  // namespace ActsExamples
