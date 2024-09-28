#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"    
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Seeding/PathSeeder.hpp"

#include "Acts/Detector/Detector.hpp"

#include "Acts/MagneticField/ConstantBField.hpp"

namespace ActsExamples {

namespace Experimental {

class SensitiveGridConstructor {
    public:
        using AxisType = Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Open>;
        using GridType = Acts::Grid<std::vector<Acts::SourceLink>, AxisType, AxisType>;

        struct Config {
            /// Detector to construct the grid for
            std::vector<const Acts::Surface*> surfaces;
            /// Binning of the grid
            std::pair<size_t, size_t> bins = {50, 50};
        };

        SensitiveGridConstructor(Config config) : m_cfg(std::move(config)) {}

        std::unordered_map<
            Acts::GeometryIdentifier, GridType> constructGrid(
                const Acts::GeometryContext& /*geoCtx*/,
                const MeasurementContainer& measurements) const {
                    std::unordered_map<Acts::GeometryIdentifier, GridType> gridLookup;

                    // Construct a binned grid for each sensitive surface
                    for (auto surface : m_cfg.surfaces) {
                        if (!surface->geometryId().sensitive()) {
                            continue;
                        }
                        if (surface->bounds().type() != Acts::SurfaceBounds::eRectangle) {
                            throw std::invalid_argument("Only rectangular bounds are supported");
                        }

                        std::vector<Acts::ActsScalar> bounds = surface->bounds().values();

                        Acts::ActsScalar xMin = bounds.at(0);
                        Acts::ActsScalar xMax = bounds.at(2);

                        Acts::ActsScalar yMin = bounds.at(1);
                        Acts::ActsScalar yMax = bounds.at(3);
    
                        AxisType xAxis(xMin, xMax, m_cfg.bins.first);
                        AxisType yAxis(yMin, yMax, m_cfg.bins.second);

                        GridType grid(std::make_tuple(xAxis, yAxis));

                        gridLookup.insert({surface->geometryId(), grid});
                    }
                    // Fill the grids with source links
                    for (auto meas : measurements) {
                        Acts::ActsScalar locX = meas.parameters()[Acts::eBoundLoc0];
                        Acts::ActsScalar locY = meas.parameters()[Acts::eBoundLoc1];
                        Acts::GeometryIdentifier geoId = 
                            meas.sourceLink().get<IndexSourceLink>().geometryId();

                        auto bin = gridLookup.at(geoId).localBinsFromPosition(
                            Acts::Vector2(locX, locY));
                        gridLookup.at(geoId).atLocalBins(bin).push_back(meas.sourceLink());
                    }

                    return gridLookup;
        }

    private:
        Config m_cfg;
};

class StraightLineIntersectionFinder {
    public:
        Acts::ActsScalar m_tol = 1e-4;

        std::vector<const Acts::Surface*> m_surfaces;

        // Find the intersections along the path
        // and return them in the order of the path
        // length
        std::vector<std::pair<Acts::GeometryIdentifier, Acts::Vector3>> operator()(
            const Acts::GeometryContext& geoCtx, 
            const Acts::FreeTrackParameters& trackParameters) const {
                Acts::Vector3 position = trackParameters.position();
                Acts::Vector3 direction = trackParameters.direction();

                std::vector<std::pair<Acts::GeometryIdentifier, Acts::Vector3>> sIntersections;
                // Intersect the surfaces
                for (auto& surface : m_surfaces) {
                    // Get the intersection
                    auto sMultiIntersection = surface->intersect(
                        geoCtx, position, direction,
                        Acts::BoundaryTolerance::AbsoluteCartesian(m_tol, m_tol));
                
                    // Take the closest
                    auto closestForward = sMultiIntersection.closestForward();
                
                    // Store if the intersection is reachable
                    if (closestForward.status() == Acts::IntersectionStatus::reachable &&
                        closestForward.pathLength() > 0.0) {
                            sIntersections.push_back(
                                {closestForward.object()->geometryId(), closestForward.position()});
                            continue;
                    }
                }
                return sIntersections;
        }
};

class PathSeedingAlgorithm : public IAlgorithm {
    public:
        /// @brief The nested configuration struct
        struct Config {
            /// Path seeder
            std::shared_ptr<Acts::Experimental::PathSeeder> seeder;
            /// SourceLink grid
            std::shared_ptr<SensitiveGridConstructor> sourceLinkGridConstructor;
            /// Input Measurement collection
            std::string inputMeasurements = "InputMeasurements";
            /// Output SourceLink collection
            std::string outputSourceLinks = "OutputSourceLinks";
            /// Output TrackParameters collection
            std::string outputTrackParameters = "OutputTrackParameters";
        };

        /// @brief Constructor
        PathSeedingAlgorithm(Config config, Acts::Logging::Level level)
            : IAlgorithm("PathSeedingAlgorithm", level),
            m_cfg(std::move(config)) {
                if (!m_cfg.seeder) {
                    throw std::invalid_argument("Missing path seeder");
                }
                if (!m_cfg.sourceLinkGridConstructor) {
                    throw std::invalid_argument("Missing source link grid constructor");
                }

                m_inputMeasurements.initialize(m_cfg.inputMeasurements);
                m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
                m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
        }
        ~PathSeedingAlgorithm() = default;

        /// @brief The execute method        
        ProcessCode execute(const AlgorithmContext& ctx) const override {
            // Get the input measurements
            // and source links
            const auto& measurements = m_inputMeasurements(ctx);
            
            auto gridLookup = 
                m_cfg.sourceLinkGridConstructor->constructGrid(
                    ctx.geoContext, measurements);

            std::vector<Acts::Experimental::PathSeeder::PathSeed> pathSeeds;
            m_cfg.seeder->findSeeds(ctx.geoContext, gridLookup, pathSeeds);

            TrackParametersContainer outTrackParameters;
            IndexSourceLinkContainer outSourceLinks;
            for (auto [ipPars, sls] : pathSeeds) {
                outTrackParameters.push_back(ipPars);

                for (auto& sl : sls) {
                    auto idxSl = sl.get<IndexSourceLink>();
                    outSourceLinks.insert(idxSl);
                }
            }

            m_outputSourceLinks(ctx, std::move(outSourceLinks));
            m_outputTrackParameters(ctx, std::move(outTrackParameters));
            
            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;

        ReadDataHandle<MeasurementContainer> m_inputMeasurements
            {this, "InputMeasurements"};

        WriteDataHandle<IndexSourceLinkContainer> m_outputSourceLinks
            {this, "OutputSeeds"};

        WriteDataHandle<TrackParametersContainer> m_outputTrackParameters
            {this, "OutputTrackParameters"};
};

}  // namespace Experimental

}  // namespace ActsExamples
