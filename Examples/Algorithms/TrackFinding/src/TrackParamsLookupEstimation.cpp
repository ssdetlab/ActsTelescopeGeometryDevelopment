// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackParamsLookupEstimation.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

ActsExamples::TrackParamsLookupEstimation::TrackParamsLookupEstimation(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("TrackParamsLookupEstimation", level),
    m_cfg(std::move(config)) {
        for (const auto& [geoId, refSurface] : m_cfg.refLayers) {
            auto bounds = 
                dynamic_cast<const Acts::RectangleBounds*>(&refSurface->bounds());

            if (bounds == nullptr) {
                throw std::invalid_argument("Only rectangle bounds supported");
            }

            auto halfX = bounds->halfLengthX();
            auto halfY = bounds->halfLengthY();

            LookupAxisGen axisGen{
                {-halfX, halfX}, m_cfg.bins.first, 
                {-halfY, halfY}, m_cfg.bins.second};

            TrackParamsLookupAccumulator::Config accConfig{
                axisGen};

            m_accumulators[refSurface->geometryId()] =
                std::make_unique<TrackParamsLookupAccumulator>(accConfig);
        }
        
        m_inputParticles.initialize(m_cfg.inputParticles);
        m_inputSimHits.initialize(m_cfg.inputHits);
}

ActsExamples::TrackParamsLookupEstimation::~TrackParamsLookupEstimation() {
    ActsExamples::Lookup lookup;

    for (auto& [id, acc] : m_accumulators) {
        lookup.insert({id, acc->finalizeLookup()});
    }
    for (auto& writer : m_cfg.trackLookupGridWriters) {
        writer->writeLookup(lookup);
    }
};

ActsExamples::ProcessCode 
ActsExamples::TrackParamsLookupEstimation::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
        const auto& particles = m_inputParticles(ctx);
        const auto& hits = m_inputSimHits(ctx);
        
        for (const auto& [geoId, refSurface] : m_cfg.refLayers) {
            auto refLayerHits = hits.equal_range(
                geoId);
    
            for (auto hit = refLayerHits.first; hit != refLayerHits.second; ++hit) {
                const auto& id = hit->particleId();
                const auto& particle = particles.find(id);
            
                if (particle == particles.end()) {
                    throw std::invalid_argument("Particle not found");
                }
            
                auto refLayerPars = Acts::CurvilinearTrackParameters(
                    hit->fourPosition(), 
                    hit->direction(), 
                    particle->qOverP(), 
                    std::nullopt,
                    particle->hypothesis());
            
                auto ipPars = Acts::CurvilinearTrackParameters(
                    particle->fourPosition(), 
                    particle->direction(), 
                    particle->qOverP(),
                    std::nullopt, 
                    particle->hypothesis());
            
                auto localPos = refSurface->globalToLocal(
                    ctx.geoContext, 
                    hit->position(), 
                    Acts::Vector3{0, 1, 0}).value();
    
                m_accumulators.at(geoId)->addTrack(
                    ipPars, 
                    refLayerPars, 
                    localPos);
            }
        }
    
        return ActsExamples::ProcessCode::SUCCESS;
}
