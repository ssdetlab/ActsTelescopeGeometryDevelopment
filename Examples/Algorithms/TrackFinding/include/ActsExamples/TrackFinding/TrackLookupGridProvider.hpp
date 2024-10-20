#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/TrackFinding/ITrackLookupGridReader.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

namespace ActsExamples {

namespace Experimental {

class TrackLookupGridProvider {
    public:
        struct Config {
            /// Track lookup grid reader
            std::shared_ptr<ITrackLookupGridReader> trackLookupGridReader;
            /// Ip lookup grid path
            std::string ipLookupGridPath;
            /// Reference lookup grid path
            std::string refLookupGridPath;
        };
        
        TrackLookupGridProvider(const Config& config)
            : m_cfg(std::move(config)),
            m_ipGrid(std::make_shared<TrackLookupGrid>(
                m_cfg.trackLookupGridReader->readLookup(m_cfg.ipLookupGridPath))),
            m_refGrid(std::make_shared<TrackLookupGrid>(
                m_cfg.trackLookupGridReader->readLookup(m_cfg.refLookupGridPath))) {}

        void linkMeasurements(const MeasurementContainer& measurements) {
            m_measurements = std::make_shared<MeasurementContainer>(measurements);
        }

        /// Lookup the track parameters at a given position
        std::pair<Acts::CurvilinearTrackParameters, Acts::CurvilinearTrackParameters> lookup(
            const Acts::GeometryContext& gctx, const Acts::SourceLink& pivot) const {
                auto idxSl = pivot.get<ActsExamples::IndexSourceLink>();

                auto measurement = m_measurements->at(idxSl.index());

                assert(measurement.contains(Acts::eBoundLoc0) &&
                    "Measurement does not contain the required bound loc0");
                assert(measurement.contains(Acts::eBoundLoc1) &&
                    "Measurement does not contain the required bound loc1");

                Acts::Vector2 localPos{
                    measurement.parameters()[Acts::eBoundLoc0],
                    measurement.parameters()[Acts::eBoundLoc1]};

                auto bin = m_ipGrid->localBinsFromPosition(localPos);

                return {*m_ipGrid->atLocalBins(bin), *m_refGrid->atLocalBins(bin)};
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;

        std::shared_ptr<TrackLookupGrid> m_ipGrid;
        std::shared_ptr<TrackLookupGrid> m_refGrid;

        std::shared_ptr<MeasurementContainer> m_measurements;
};


class DummyGridTest : public IAlgorithm {
    public:
        struct Config {
            /// Reference tracking layers
            const Acts::Surface* refLayer;
            /// Track lookup grid provider
            std::shared_ptr<TrackLookupGridProvider> lookup;
            /// Input SourceLink collection
            std::string inputHits = "InputHits";
            /// Input Measurement collection
            std::string inputParticles = "InputParticles";
        };
    
        DummyGridTest(const Config& config, Acts::Logging::Level level)
            : IAlgorithm("DummyGridTest", level),
            m_cfg(std::move(config)) {
                m_inputSimHits.initialize(m_cfg.inputHits);
                m_inputParticles.initialize(m_cfg.inputParticles);
        }
    
        ProcessCode execute(const AlgorithmContext& ctx) const final {
            auto particles = m_inputParticles(ctx);
            auto hits = m_inputSimHits(ctx);

            MeasurementContainer measurements;
            std::vector<Acts::SourceLink> sourceLinks;
            std::vector<Acts::CurvilinearTrackParameters> ipTrackParameters;
            std::vector<Acts::CurvilinearTrackParameters> refLayerTrackParameters;
            for (const auto& hit : hits) {
                if (hit.geometryId() != m_cfg.refLayer->geometryId()) {
                    continue;
                }
                
                IndexSourceLink isl{hit.geometryId(), measurements.size()};
                Acts::SourceLink sl{isl};
                sourceLinks.push_back(sl);
                
                std::array<Acts::BoundIndices, 2u> indices = {
                    Acts::BoundIndices{Acts::eBoundLoc0},
                    Acts::BoundIndices{Acts::eBoundLoc1}
                };

                Acts::Vector2 localPos = m_cfg.refLayer->globalToLocal(
                    ctx.geoContext, hit.position(), Acts::Vector3{0, 1, 0}).value();
                measurements.emplaceMeasurement<2u>(sl, indices, localPos, Acts::SquareMatrix2::Identity());

                auto id = hit.particleId();
                auto particle = particles.find(id);
            
                if (particle == particles.end()) {
                    throw std::invalid_argument("Particle not found");
                }
            
                auto refLayerPars = Acts::CurvilinearTrackParameters(
                    hit.fourPosition(), hit.direction(), particle->qOverP(), std::nullopt,
                    particle->hypothesis());
            
                auto ipPars = Acts::CurvilinearTrackParameters(
                    particle->fourPosition(), particle->direction(), particle->qOverP(),
                    std::nullopt, particle->hypothesis());
            
                ipTrackParameters.push_back(std::move(ipPars));
                refLayerTrackParameters.push_back(std::move(refLayerPars));
            }
            m_cfg.lookup->linkMeasurements(measurements);

            Acts::Vector4 fourPositionErrIp{0., 0., 0., 0.};
            Acts::Vector4 fourPositionErrRef{0., 0., 0., 0.};

            Acts::Vector3 directionErrIp{0., 0., 0.};
            Acts::Vector3 directionErrRef{0., 0., 0.};

            double qOverPErrIp = 0;
            double qOverPErrRef = 0;
            for (int i = 0; i < sourceLinks.size(); ++i) {
                auto sl = sourceLinks.at(i);
                std::cout << "----------------------------------------------------" << std::endl;
                auto [ipParsEst, refLayerParsEst] = m_cfg.lookup->lookup(ctx.geoContext, sl);

                fourPositionErrIp += 
                    (ipParsEst.fourPosition() - ipTrackParameters.at(i).fourPosition()).cwiseAbs() / 
                        ipTrackParameters.at(i).fourPosition().norm();
                fourPositionErrRef +=
                    (refLayerParsEst.fourPosition() - refLayerTrackParameters.at(i).fourPosition()).cwiseAbs() / 
                        refLayerTrackParameters.at(i).fourPosition().norm();
                    
                directionErrIp +=
                    (ipParsEst.direction() - ipTrackParameters.at(i).direction()).cwiseAbs() / 
                        ipTrackParameters.at(i).direction().norm();
                directionErrRef +=
                    (refLayerParsEst.direction() - refLayerTrackParameters.at(i).direction()).cwiseAbs() / 
                        refLayerTrackParameters.at(i).direction().norm();

                qOverPErrIp += (ipParsEst.qOverP() - ipTrackParameters.at(i).qOverP()) / ipTrackParameters.at(i).qOverP();
                qOverPErrRef += (refLayerParsEst.qOverP() - refLayerTrackParameters.at(i).qOverP()) / refLayerTrackParameters.at(i).qOverP();
            }
            fourPositionErrIp /= sourceLinks.size();
            fourPositionErrRef /= sourceLinks.size();
            directionErrIp /= sourceLinks.size();
            directionErrRef /= sourceLinks.size();
            qOverPErrIp /= sourceLinks.size();
            qOverPErrRef /= sourceLinks.size();

            std::cout << "--- IP TRACK PARAMETERS ---" << std::endl;
            std::cout << "Four Position Error: " << fourPositionErrIp.transpose() << std::endl;
            std::cout << "Direction Error: " << directionErrIp.transpose() << std::endl;
            std::cout << "qOverP Error: " << qOverPErrIp << std::endl;

            std::cout << "--- REF LAYER TRACK PARAMETERS ---" << std::endl;
            std::cout << "Four Position Error: " << fourPositionErrRef.transpose() << std::endl;
            std::cout << "Direction Error: " << directionErrRef.transpose() << std::endl;
            std::cout << "qOverP Error: " << qOverPErrRef << std::endl;


            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;

        ReadDataHandle<SimParticleContainer> m_inputParticles{
            this,
            "InputSimParticles"};
        
        ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
};


}  // namespace Experimental

}  // namespace ActsExamples
