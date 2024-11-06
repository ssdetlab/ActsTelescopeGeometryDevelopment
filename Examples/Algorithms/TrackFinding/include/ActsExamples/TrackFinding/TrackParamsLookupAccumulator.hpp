// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/TrackFinding/TrackParamsLookupTable.hpp"

namespace ActsExamples {

class TrackParamsLookupAccumulator {
    public:
        struct Config {
            /// Axis generator
            LookupAxisGen axisGen;
        };

        TrackParamsLookupAccumulator(const Config& config)
            : m_cfg(std::move(config)),
            m_ipGrid(m_cfg.axisGen()),
            m_refGrid(m_cfg.axisGen()) {}

        void addTrack(
            const Acts::CurvilinearTrackParameters& ipTrackParameters,
            const Acts::CurvilinearTrackParameters& refTrackParameters,
            const Acts::Vector2& position);

        LookupGrid finalizeLookup();

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;

        std::mutex m_writeMutex;

        LookupAccumGrid m_ipGrid;
        LookupAccumGrid m_refGrid;
};

}  // namespace ActsExamples
