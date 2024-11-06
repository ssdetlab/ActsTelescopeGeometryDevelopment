// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/TrackFinding/TrackParamsLookupProvider.hpp"

namespace ActsExamples {

class TrackParamsLookupValidation : public IAlgorithm {
    public:
        struct Config {
            /// Reference tracking layers
            std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> refLayers;
            /// Track lookup grid provider
            std::shared_ptr<TrackParamsLookupProvider> lookup;
            /// Input SourceLink collection
            std::string inputHits = "InputHits";
            /// Input Measurement collection
            std::string inputParticles = "InputParticles";
            /// Output IpPars collection
            std::string outputIpPars = "OutputIpPars";
            /// Output RefLayerPars collection
            std::string outputRefLayerPars = "OutputRefLayerPars";
            /// Output IpParsEst collection
            std::string outputIpParsEst = "OutputIpParsEst";
            /// Output RefLayerParsEst collection
            std::string outputRefLayerParsEst = "OutputRefLayerParsEst";
        };
    
        TrackParamsLookupValidation(const Config& config, Acts::Logging::Level level);
    
        ProcessCode execute(const AlgorithmContext& ctx) const final;

        // Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;

        ReadDataHandle<SimParticleContainer> m_inputParticles{
            this,
            "InputSimParticles"};
        
        ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};

        WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>> m_outputIpPars{
            this, "OutputIpPars"};
        WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>> m_outputRefLayerPars{
            this, "OutputRefLayerPars"};
        
        WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>> m_outputIpParsEst{
            this, "OutputIpParsEst"};
        WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>> m_outputRefLayerParsEst{
            this, "OutputRefLayerParsEst"};
};


}  // namespace ActsExamples
