#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/TrackFinding/ITrackLookupGridWriter.hpp"
#include "Acts/Utilities/Interpolation.hpp"

#include <fstream>

namespace ActsExamples {

namespace Experimental {

class TrackLookupGridAccumulator {
    public:
        struct Config {
            /// Binning of the grid
            std::pair<std::size_t, std::size_t> bins;
            /// Grid bounds
            std::pair<Acts::ActsScalar, Acts::ActsScalar> xBounds;
            std::pair<Acts::ActsScalar, Acts::ActsScalar> yBounds;
        };

        // struct InterpolableTrackParameters {
            // std::shared_ptr<Acts::CurvilinearTrackParameters> trackParameters = nullptr;

            // InterpolableTrackParameters() = default;
            // InterpolableTrackParameters(
                // Acts::Vector4 fourPos,
                // Acts::Vector3 dir,
                // Acts::ActsScalar qop,
                // std::optional<Acts::BoundMatrix> cov = std::nullopt, 
                // Acts::ParticleHypothesis particle = Acts::ParticleHypothesis::electron())
                // : trackParameters(std::make_shared<Acts::CurvilinearTrackParameters>(
                    // fourPos, dir, qop, cov, particle)) {}

            // InterpolableTrackParameters(const InterpolableTrackParameters& other) {
                // trackParameters = other.trackParameters;
            // }
            // InterpolableTrackParameters(InterpolableTrackParameters&& other) {
                // trackParameters = std::move(other.trackParameters);
            // }
            // InterpolableTrackParameters& operator=(const InterpolableTrackParameters& other) {
                // trackParameters = other.trackParameters;
                // return *this;
            // }
            // InterpolableTrackParameters& operator=(InterpolableTrackParameters&& other) {
                // trackParameters = std::move(other.trackParameters);
                // return *this;
            // }

            // InterpolableTrackParameters operator+(const InterpolableTrackParameters& other) const {
                // return InterpolableTrackParameters(
                    // trackParameters->fourPosition() + other.trackParameters->fourPosition(),
                    // trackParameters->direction() + other.trackParameters->direction(),
                    // trackParameters->qOverP() + other.trackParameters->qOverP(),
                    // trackParameters->covariance(), trackParameters->particleHypothesis());
            // }

            // InterpolableTrackParameters operator*(Acts::ActsScalar scalar) const {
                // return InterpolableTrackParameters(
                    // trackParameters->fourPosition() * scalar,
                    // trackParameters->direction() * scalar,
                    // trackParameters->qOverP() * scalar,
                    // trackParameters->covariance(), trackParameters->particleHypothesis());
            // }

            // friend InterpolableTrackParameters operator*(Acts::ActsScalar scalar, const InterpolableTrackParameters& tp) {
                // return tp * scalar;
            // }

            // friend std::ostream& operator<<(std::ostream& os, const InterpolableTrackParameters& tp) {
                // os << "FourPosition: " << tp.trackParameters->fourPosition().transpose() << std::endl;
                // os << "Direction: " << tp.trackParameters->direction().transpose() << std::endl;
                // os << "QOverP: " << tp.trackParameters->qOverP() << std::endl;
                // return os;
            // }
        // };

        TrackLookupGridAccumulator(const Config& config)
            : m_cfg(std::move(config)),
            m_ipGrid(std::make_tuple(
                AxisType(m_cfg.xBounds.first, m_cfg.xBounds.second, m_cfg.bins.first),
                AxisType(m_cfg.yBounds.first, m_cfg.yBounds.second, m_cfg.bins.second))),
            m_refGrid(std::make_tuple(
                AxisType(m_cfg.xBounds.first, m_cfg.xBounds.second, m_cfg.bins.first),
                AxisType(m_cfg.yBounds.first, m_cfg.yBounds.second, m_cfg.bins.second))) {
                    std::cout << "\n\n\n BINS " << m_cfg.bins.first << " " << m_cfg.bins.second << std::endl;
                    std::cout << "BOUNDS " << m_cfg.xBounds.first << " " << m_cfg.xBounds.second << std::endl;
                    std::cout << "BOUNDS " << m_cfg.yBounds.first << " " << m_cfg.yBounds.second << std::endl;
                    for (auto it = m_ipGrid.begin(); it != m_ipGrid.end(); ++it) {
                        auto bin = it.localBinsIndices();
                        std::cout << "BIN " << bin[0] << " " << bin[1] << std::endl;
                        std::cout << "CENTER " << m_ipGrid.binCenter(bin).at(0) << " " << m_ipGrid.binCenter(bin).at(1) << std::endl;
                    }
                    std::cout << "\n\n\n";
                }

        void addTrack(
            const Acts::CurvilinearTrackParameters& ipTrackParameters,
            const Acts::CurvilinearTrackParameters& refTrackParameters,
            const Acts::Vector2& position) {
                std::lock_guard<std::mutex> lock(m_writeMutex);
                auto bin = m_ipGrid.localBinsFromPosition(position);
                
                std::cout << "ADDING TRACK AT " << position.transpose() << std::endl;
                std::cout << "refTrackParameters: " << refTrackParameters.fourPosition().transpose() << std::endl;
                std::cout << "ipTrackParameters: " << ipTrackParameters.fourPosition().transpose() << std::endl;
                std::cout << "BIN " << m_ipGrid.localBinsFromPosition(position)[0] << " " << m_ipGrid.localBinsFromPosition(position)[1] << std::endl;
                m_ipGrid.atLocalBins(bin).push_back(ipTrackParameters);
                m_refGrid.atLocalBins(bin).push_back(refTrackParameters);
                std::cout << "SIZE " << m_ipGrid.atPosition(position).size() << std::endl;
        }

        std::pair<TrackLookupGrid, TrackLookupGrid> finalizeLookup() {
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

            TrackLookupGrid ipGrid(std::make_tuple(
                AxisType(m_cfg.xBounds.first, m_cfg.xBounds.second, m_cfg.bins.first),
                AxisType(m_cfg.yBounds.first, m_cfg.yBounds.second, m_cfg.bins.second)));

            TrackLookupGrid refGrid(std::make_tuple(
                AxisType(m_cfg.xBounds.first, m_cfg.xBounds.second, m_cfg.bins.first),
                AxisType(m_cfg.yBounds.first, m_cfg.yBounds.second, m_cfg.bins.second)));

            for (auto it = m_ipGrid.begin(); it != m_ipGrid.end(); ++it) {
                auto bin = it.localBinsIndices();
                std::cout << "BIN " << bin[0] << " " << bin[1] << std::endl;
                if (m_ipGrid.atLocalBins(bin).empty() || 
                    m_refGrid.atLocalBins(bin).empty()) {
                        std::cout << "EMPTY" << std::endl;
                        continue;
                }
                
                ipGrid.atLocalBins(bin) = 
                    std::make_unique<
                        Acts::CurvilinearTrackParameters>(
                            meanTrack(m_ipGrid.atLocalBins(bin)));
                refGrid.atLocalBins(bin) = 
                    std::make_unique<
                        Acts::CurvilinearTrackParameters>(
                            meanTrack(m_refGrid.atLocalBins(bin)));
            }

            return {std::move(ipGrid), std::move(refGrid)};
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;

        std::mutex m_writeMutex;

        Acts::Grid<
            std::vector<Acts::CurvilinearTrackParameters>,
            AxisType, AxisType> m_ipGrid;
        Acts::Grid<
            std::vector<Acts::CurvilinearTrackParameters>,
            AxisType, AxisType> m_refGrid;
};

class TrackLookupEstimationAlgorithm : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Reference tracking layers
            std::vector<const Acts::Surface*> refLayers;
            /// Input SourceLink collection
            std::string inputHits = "InputHits";
            /// Input Measurement collection
            std::string inputParticles = "InputParticles";
            /// Output SourceLink collection
            std::string outputIPTrackParameters = "OutputIPTrackParameters";
            /// Output TrackParameters collection
            std::string outputRefLayerTrackParameters = "OutputRefLayerTrackParameters";
            /// Track lookup grid accumulator
            std::shared_ptr<TrackLookupGridAccumulator> trackLookupGridAccumulator;
            /// Track lookup grid writer
            std::vector<std::shared_ptr<ITrackLookupGridWriter>>
                trackLookupGridWriters{};
        };
        
        /// @brief Constructor
        TrackLookupEstimationAlgorithm(Config config, Acts::Logging::Level level)
            : IAlgorithm("TrackLookupEstimationAlgorithm", level),
                m_cfg(std::move(config)) {
            m_inputParticles.initialize(m_cfg.inputParticles);
            m_inputSimHits.initialize(m_cfg.inputHits);
            m_outputIPTrackParameters.initialize(m_cfg.outputIPTrackParameters);
            m_outputRefLayerTrackParameters.initialize(
                m_cfg.outputRefLayerTrackParameters);
        }
        ~TrackLookupEstimationAlgorithm() {
            auto [ipGrid, refGrid] =
                m_cfg.trackLookupGridAccumulator->finalizeLookup();
            for (auto& writer : m_cfg.trackLookupGridWriters) {
                writer->writeLookup(ipGrid, refGrid);
            }
        };

        /// @brief The execute method
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            const auto& particles = m_inputParticles(ctx);
            const auto& hits = m_inputSimHits(ctx);
        
            std::cout << "REFERENCE LAYERS: " << m_cfg.refLayers.size() << std::endl;
            for (const auto& layer : m_cfg.refLayers) {
                std::cout << "LAYER: " << layer->geometryId() << std::endl;
                std::cout << "LAYER: " << layer->center(ctx.geoContext).transpose() << std::endl;
            }

            std::cout << "PARTICLES: " << particles.size() << std::endl;
            std::cout << "HITS: " << hits.size() << std::endl;

            for (const auto& particle : particles) {
                std::cout << "PARTICLE: " << particle.numberOfHits() << std::endl;
            }

            TrackParametersContainer ipTrackParameters;
            TrackParametersContainer refLayerTrackParameters;
            for (auto& hit : hits) {
                const auto surfIt =
                    std::find_if(m_cfg.refLayers.begin(), m_cfg.refLayers.end(),
                                [&hit](const Acts::Surface* surface) {
                                    return surface->geometryId() == hit.geometryId();
                                });
                if (surfIt == m_cfg.refLayers.end()) {
                    continue;
                }
                std::cout << "FOUND SURFACE" << std::endl;
        
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
            
                m_cfg.trackLookupGridAccumulator->addTrack(
                    ipPars, 
                    refLayerPars, 
                    (*surfIt)->globalToLocal(
                        ctx.geoContext, 
                        hit.position(),
                        Acts::Vector3{0, 1, 0}).value());
    
                ipTrackParameters.push_back(std::move(ipPars));
                refLayerTrackParameters.push_back(std::move(refLayerPars));
            }

            m_outputIPTrackParameters(ctx, std::move(ipTrackParameters));
            m_outputRefLayerTrackParameters(ctx, std::move(refLayerTrackParameters));

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;
        
        std::mutex m_writeMutex;

        ReadDataHandle<SimParticleContainer> m_inputParticles{
            this,
            "InputSimParticles"};
        
        ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
        
        WriteDataHandle<TrackParametersContainer> m_outputIPTrackParameters{
            this, "OutputIPTrackParameters"};
        
        WriteDataHandle<TrackParametersContainer> m_outputRefLayerTrackParameters{
            this, "OutputRefLayerTrackParameters"};
};

}  // namespace Experimental

}  // namespace ActsExamples
