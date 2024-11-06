// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackParamsLookupAccumulator.hpp"

#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridIterator.hpp"
#include "Acts/EventData/TrackParameters.hpp"

void ActsExamples::TrackParamsLookupAccumulator::addTrack(
    const Acts::CurvilinearTrackParameters& ipTrackParameters,
    const Acts::CurvilinearTrackParameters& refTrackParameters,
    const Acts::Vector2& position) {
        std::lock_guard<std::mutex> lock(m_writeMutex);
        auto bin = m_ipGrid.localBinsFromPosition(position);
        
        m_ipGrid.atLocalBins(bin).push_back(ipTrackParameters);
        m_refGrid.atLocalBins(bin).push_back(refTrackParameters);
}

ActsExamples::LookupGrid 
ActsExamples::TrackParamsLookupAccumulator::finalizeLookup() {
    auto meanTrack = [](const std::vector<Acts::CurvilinearTrackParameters>& tracks) {
        Acts::Vector4 fourPosition = Acts::Vector4::Zero();
        Acts::Vector3 direction = Acts::Vector3::Zero();
        Acts::ActsScalar qOverP = 0;
        for (const auto& track : tracks) {
            fourPosition += track.fourPosition();
            direction += track.direction();
            qOverP += track.qOverP();
        }
        fourPosition /= tracks.size();
        direction /= tracks.size();
        qOverP /= tracks.size();
    
        return Acts::CurvilinearTrackParameters(
            fourPosition, 
            direction, 
            qOverP, 
            std::nullopt, 
            tracks.front().particleHypothesis());
    };

    ActsExamples::LookupGrid lookupGrid(m_cfg.axisGen());

    for (auto it = m_ipGrid.begin(); it != m_ipGrid.end(); ++it) {
        auto bin = it.localBinsIndices();
        if (m_ipGrid.atLocalBins(bin).empty() || 
            m_refGrid.atLocalBins(bin).empty()) {
                continue;
        }
        
        lookupGrid.atLocalBins(bin) = 
            std::make_pair(
                std::make_shared<Acts::CurvilinearTrackParameters>(
                    meanTrack(m_ipGrid.atLocalBins(bin))),
                std::make_shared<Acts::CurvilinearTrackParameters>(
                    meanTrack(m_refGrid.atLocalBins(bin))));
    }

    return lookupGrid;
}
