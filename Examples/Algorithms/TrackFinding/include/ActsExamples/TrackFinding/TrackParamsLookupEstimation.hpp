// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once


#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/TrackFinding/ITrackParamsLookupWriter.hpp"
#include "ActsExamples/TrackFinding/TrackParamsLookupAccumulator.hpp"

namespace ActsExamples {

class TrackParamsLookupEstimation : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Reference tracking layers
            std::unordered_map<
                Acts::GeometryIdentifier, 
                const Acts::Surface*> refLayers;
            /// Binning of the grid to be emposed
            /// onto the reference layers
            std::pair<std::size_t, std::size_t> bins;
            /// Input SourceLink collection
            std::string inputHits = "InputHits";
            /// Input Measurement collection
            std::string inputParticles = "InputParticles";
            /// Track lookup grid writer
            std::vector<std::shared_ptr<ITrackParamsLookupWriter>>
                trackLookupGridWriters{};
        };
        
        /// @brief Constructor
        TrackParamsLookupEstimation(const Config& config, Acts::Logging::Level level);

        /// @brief Destructor
        ~TrackParamsLookupEstimation();

        /// @brief The execute method
        ProcessCode execute(const AlgorithmContext& ctx) const override;

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;
        
        ReadDataHandle<SimParticleContainer> m_inputParticles{
            this,
            "InputSimParticles"};
        
        ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};

        std::unordered_map<Acts::GeometryIdentifier, 
           std::unique_ptr<TrackParamsLookupAccumulator>> m_accumulators;
};

}  // namespace ActsExamples
