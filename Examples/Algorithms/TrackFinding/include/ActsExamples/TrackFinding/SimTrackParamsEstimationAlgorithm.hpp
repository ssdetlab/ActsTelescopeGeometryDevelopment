#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"    
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "Acts/EventData/TrackParameters.hpp"   

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Seeding/PathSeeder.hpp"

#include "Acts/Detector/Detector.hpp"

#include "Acts/MagneticField/ConstantBField.hpp"

namespace ActsExamples {

namespace Experimental {

using AxisType = Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Open>;

using TrackParameters = Acts::CurvilinearTrackParameters; 
using TrackLookupGrid = Acts::Grid<TrackParameters, AxisType, AxisType>;

using TrackParametersContainer = std::vector<TrackParameters>;
using TrackLookupContainerGrid = Acts::Grid<TrackParametersContainer, AxisType, AxisType>;

class ITrackLookupGridWriter {
    public:
        /// Virtual Destructor
        virtual ~ITrackLookupGridWriter() = default;

        /// Writer method
        ///
        /// @param detMaterial the detector material maps
        virtual void writeLookup(const TrackLookupGrid& trackGrid) = 0;
};

class TrackLookupGridJsonWriter : public ITrackLookupGridWriter {
    public:
        struct Config {
            /// Output file name
            std::string fileName = "track_lookup_grid";
        };

        TrackLookupGridJsonWriter(Config config) : m_cfg(std::move(config)) {}

        void writeLookup(const TrackLookupGrid& trackGrid) override {
            std::ofstream file(m_cfg.fileName);
            if (!file.is_open()) {
                throw std::invalid_argument("Could not open file");
            }

            // Acts::JsonGeometryConverter::writeGrid(file, trackGrid);
        }

    private:
        Config m_cfg;
};

class TrackLookupGridAccumulator {
    public: 
        struct Config {
            /// Binning of the grid
            std::pair<size_t, size_t> bins;
            /// Grid bounds
            std::pair<Acts::ActsScalar, Acts::ActsScalar> xBounds;
            std::pair<Acts::ActsScalar, Acts::ActsScalar> yBounds;
        };

        TrackLookupGridAccumulator(Config config) 
            : m_cfg(std::move(config)),
            m_trackContainerGrid(
                std::make_tuple(
                    AxisType(m_cfg.xBounds.first, m_cfg.xBounds.second, m_cfg.bins.first),
                    AxisType(m_cfg.yBounds.first, m_cfg.yBounds.second, m_cfg.bins.second))) {}

        void addTrack(const TrackParameters& trackParameters, const Acts::Vector2& position) {
            m_trackContainerGrid.atPosition(position).push_back(trackParameters);
        }

        TrackLookupGrid finalizeLookup() {
            auto meanTrack = [](const TrackParametersContainer& container) {
                Acts::Vector4 meanVertex;
                Acts::Vector3 meanDirection;
                Acts::ActsScalar meanQOverP;

                for (const auto& track : container) {
                    meanVertex += track.fourPosition();
                    meanDirection += track.direction();
                    meanQOverP += track.qOverP();
                }

                meanVertex /= container.size();
                meanDirection /= container.size();
                meanQOverP /= container.size();

                return TrackParameters(
                    meanVertex, 
                    meanDirection, 
                    meanQOverP, 
                    std::nullopt, 
                    container.front().particleHypothesis());
            };

            TrackLookupGrid trackGrid(
                std::make_tuple(
                    AxisType(m_cfg.xBounds.first, m_cfg.xBounds.second, m_cfg.bins.first),
                    AxisType(m_cfg.yBounds.first, m_cfg.yBounds.second, m_cfg.bins.second)));

            for (auto it = m_trackContainerGrid.begin(); it != m_trackContainerGrid.end(); it++) {
                auto bin = m_trackContainerGrid.atLocalBins(it.localBinsIndices());
                if (bin.empty()) {
                    continue;
                }

                trackGrid.atLocalBins(it.localBinsIndices()) = meanTrack(bin);
            }

            return trackGrid;
        }

    private:
        Config m_cfg;

        TrackLookupContainerGrid m_trackContainerGrid;
};

class SimTrackParamsEstimationAlgorithm : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Reference tracking layer IDs
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
            std::vector<std::shared_ptr<ITrackLookupGridWriter>> trackLookupGridWriters{};
        };

        /// @brief Constructor
        SimTrackParamsEstimationAlgorithm(Config config, Acts::Logging::Level level)
            : IAlgorithm("SimTrackParamsEstimationAlgorithm", level),
            m_cfg(std::move(config)) {
                m_inputParticles.initialize(m_cfg.inputParticles);
                m_inputSimHits.initialize(m_cfg.inputHits);
                m_outputIPTrackParameters.initialize(m_cfg.outputIPTrackParameters);
                m_outputRefLayerTrackParameters.initialize(m_cfg.outputRefLayerTrackParameters);
        }
        ~SimTrackParamsEstimationAlgorithm() {
            TrackLookupGrid trackGrid = m_cfg.trackLookupGridAccumulator->finalizeLookup();
            for (auto& writer : m_cfg.trackLookupGridWriters) {
                writer->writeLookup(trackGrid);
            }
        };

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            const auto& particles = m_inputParticles(ctx);
            const auto& hits = m_inputSimHits(ctx);

            for (auto& hit : hits) {
                const auto surfIt = std::find_if(
                    m_cfg.refLayers.begin(), 
                    m_cfg.refLayers.end(), 
                    [&hit](const Acts::Surface* surface) {
                        return surface->geometryId() == hit.geometryId();
                    }); 
                if (surfIt == m_cfg.refLayers.end()) {
                    continue;   
                }

                auto id = hit.particleId();
                auto particle = particles.find(id);

                if (particle == particles.end()) {
                    throw std::invalid_argument("Particle not found");
                }

                auto firstLayerPars = Acts::CurvilinearTrackParameters(
                    hit.fourPosition(),
                    hit.direction(),
                    particle->qOverP(),
                    std::nullopt,
                    particle->hypothesis()
                );

                auto ipPars = Acts::CurvilinearTrackParameters(
                    particle->fourPosition(),
                    particle->direction(),
                    particle->qOverP(),
                    std::nullopt,
                    particle->hypothesis());
                
                m_cfg.trackLookupGridAccumulator->addTrack(
                    ipPars, 
                    (*surfIt)->globalToLocal(
                        ctx.geoContext, 
                        hit.position(), 
                        Acts::Vector3{0, 1, 0}).value());
            }

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;

        ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputSimParticles"};
        
        ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};

        WriteDataHandle<TrackParametersContainer> m_outputIPTrackParameters
            {this, "OutputIPTrackParameters"};

        WriteDataHandle<TrackParametersContainer> m_outputRefLayerTrackParameters
            {this, "OutputRefLayerTrackParameters"};
};

}  // namespace Experimental

}  // namespace ActsExamples
