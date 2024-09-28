// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/GridIterator.hpp"
#include "Acts/Seeding/detail/UtilityFunctions.hpp"
#include "Acts/EventData/TrackParameters.hpp"

namespace Acts::Experimental {

/// @brief Seeding algorigthm that extracts
/// the IP parameters and sorts the source links
/// into possible track candidates
///
/// The algorithm to convert the given source links
/// into seeds -- pairs of IP parameters and the corresponding
/// source links -- as follows: First the source links
/// are sorted into a user-defined grid. Then, iteration over the source links
/// is performed. If a source link is attached to a surface that is
/// in the reference tracking layer, as defined by the user, the IP parameters
/// are estimated and the tracking layers are intersected to construct the
/// core of the "Path". The source links in the subsequent layers are then
/// added to the seed if they lie within the path width of the core.
/// Both the list of source links and the IP parameters are stored in the seed
/// struct.
///
/// @tparam axis_t Type of the axis to bin
/// the source links
///
/// @note The algorithm is designed to be used in the
/// context of a telescope-style geometry. The surfaces
/// are assumed to be planar.
///
/// @note Handling of the rotated surfaces has to happen
/// in the user-defined delegate functions.
class PathSeeder {
    public:
        using PathSeed = 
            std::pair<CurvilinearTrackParameters, std::vector<SourceLink>>;

        /// @brief Delegate to estimate the IP parameters
        /// and the momentum direction at the reference tracking layer
        ///
        /// @arg Geometry context to use
        /// @arg Pivot source link
        ///
        /// @return Pair of the track parameters at the IP and
        /// the reference tracking layer
        using TrackEstimator =
            Delegate<std::pair<CurvilinearTrackParameters, CurvilinearTrackParameters>(
                const GeometryContext&, const SourceLink&)>;

        /// @brief Delegate to find the intersections for the given pivot
        /// source link
        ///
        /// @arg The geometry context to use
        /// @arg Track parameters at the reference tracking layer
        using IntersectionLookup =
            Delegate<std::vector<std::pair<GeometryIdentifier, Vector2>>(
                const GeometryContext&, const CurvilinearTrackParameters&)>;

        /// @brief Delegate to provide the path width around
        /// the intersection point to pull the source links
        /// from the grid
        ///
        /// @arg The geometry context to use
        /// @arg The geometry identifier to use if the
        /// path width is varied across different tracking layers
        ///
        /// @return The path width in the bin0 and bin1 direction
        /// defined with respect to the surface normal
        using PathWidthLookup = Delegate<std::pair<ActsScalar, ActsScalar>(
            const GeometryContext&, const GeometryIdentifier&)>;
        
        /// @brief The nested configuration struct
        struct Config {
            /// Parameters estimator
            TrackEstimator trackEstimator;
            /// Intersection finder
            IntersectionLookup intersectionFinder;
            /// Path width provider
            PathWidthLookup pathWidthProvider;
            /// Reference layer IDs
            std::vector<GeometryIdentifier> refLayerIds;
            /// Direction of the telescope extent
            BinningValue orientation = BinningValue::binX;
        };

        /// @brief Constructor
        PathSeeder(const Config& config) : m_cfg(std::move(config)) {};
        
        /// @brief Destructor
        ~PathSeeder() = default;
        
        /// @brief Extract the IP parameters and
        /// sort the source links into the seeds
        ///
        /// @param gctx The geometry context
        /// @param sourceLinks The source links to seed
        ///
        /// @return The vector of seeds
        template <typename grid_t, typename container_t>
        void findSeeds(
            const GeometryContext& gctx,
            const std::unordered_map<
                GeometryIdentifier,grid_t>& sourceLinkGridLookup,
            container_t& seedCollection) const {
                // Get plane of the telescope
                // sensitive surfaces
                int bin0 = static_cast<int>(BinningValue::binX);
                int bin1 = static_cast<int>(BinningValue::binY);
                if (m_cfg.orientation == BinningValue::binX) {
                    bin0 = static_cast<int>(BinningValue::binY);
                    bin1 = static_cast<int>(BinningValue::binZ);
                } else if (m_cfg.orientation == BinningValue::binY) {
                    bin0 = static_cast<int>(BinningValue::binX);
                    bin1 = static_cast<int>(BinningValue::binZ);
                }
    
                // Create the seeds
                for (auto& firstGeoId : m_cfg.refLayerIds) {
                    auto firstGrid = sourceLinkGridLookup.at(firstGeoId);
    
                    for (auto it = firstGrid.begin(); it != firstGrid.end(); it++) {
                        std::vector<SourceLink> pivotSourceLinks = *it;
    
                        for (const auto& pivot : pivotSourceLinks) {
                            // Get the IP parameters
                            auto [ipParameters, firstLayerParameters] =
                                m_cfg.trackEstimator(gctx, pivot);

                            // Intersect with the surfaces
                            std::vector<std::pair<GeometryIdentifier, Vector2>> intersections =
                                m_cfg.intersectionFinder(gctx, firstLayerParameters);
    
                            // Continue if no intersections
                            if (intersections.empty()) {
                                continue;
                            }

                            // Create the seed
                            std::vector<SourceLink> seedSourceLinks;

                            // Add the pivot source link
                            seedSourceLinks.push_back(pivot);
    
                            // Iterate over the intersections
                            // and get the source links
                            // in the subsequent layers
                            for (auto& [geoId, refPoint] : intersections) {
                                // Get the path width
                                auto [pathWidth0, pathWidth1] = 
                                    m_cfg.pathWidthProvider(gctx, geoId);
                        
                                // Get the bounds of the path
                                ActsScalar top0 = refPoint[bin0] + pathWidth0;
                                ActsScalar bot0 = refPoint[bin0] - pathWidth0;
                                ActsScalar top1 = refPoint[bin1] + pathWidth1;
                                ActsScalar bot1 = refPoint[bin1] - pathWidth1;
                        
                                // Get the lookup table for the source links
                                auto grid = sourceLinkGridLookup.at(geoId);
                        
                                // Get the range of bins to search for source links
                                auto botLeftBin = grid.localBinsFromPosition(Vector2(bot0, bot1));
                                auto topRightBin = grid.localBinsFromPosition(Vector2(top0, top1));
                        
                                // Get the source links from the lookup table
                                // by iterating over the bin ranges
                                auto currentBin = botLeftBin;
                                while (currentBin.at(1) <= topRightBin.at(1)) {
                                    while (currentBin.at(0) <= topRightBin.at(0)) {
                                        auto sourceLinksToAdd = grid.atLocalBins(currentBin);
    
                                        seedSourceLinks.insert(
                                            seedSourceLinks.end(), sourceLinksToAdd.begin(),
                                            sourceLinksToAdd.end());

                                        currentBin.at(0)++;
                                    }
                                    currentBin.at(1)++;
                                    currentBin.at(0) = botLeftBin.at(0);
                                }
                            }
                            PathSeed seed = {ipParameters, seedSourceLinks};

                            // Add the seed to the collection
                            Acts::detail::pushBackOrInsertAtEnd(
                                seedCollection, seed);
                        }
                    }
                }
        }

    private:
        Config m_cfg;
};

}  // namespace Acts::Experimental
