// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/detail/CuboidalDetectorHelper.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Seeding/PathSeeder.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"

BOOST_AUTO_TEST_SUITE(PathSeeder)

using namespace Acts;
using namespace Acts::UnitLiterals;

using Axis = Acts::Axis<AxisType::Equidistant, AxisBoundaryType::Open>;
using Grid = Acts::Grid<std::vector<SourceLink>, Axis, Axis>;

// Parameters for the geometry
const ActsScalar halfY = 10.;
const ActsScalar halfZ = 10.;
const ActsScalar deltaX = 10.;
const ActsScalar deltaYZ = 1.;

using TrackContainer = Acts::TrackContainer<
    Acts::VectorTrackContainer,
    Acts::VectorMultiTrajectory,
    Acts::detail::ValueHolder>;

using TrackStateContainerBackend =
    typename TrackContainer::TrackStateContainerBackend;

/// The map(-like) container accessor
template <typename container_t>
struct TestContainerAccessor {
    using Container = container_t;
    using Key = typename container_t::key_type;
    using Value = typename container_t::mapped_type;
    
    /// This iterator adapter is needed to have the deref operator return a single
    /// source link instead of the map pair <GeometryIdentifier,SourceLink>
    struct Iterator {
        using BaseIterator = typename container_t::const_iterator;
    
        using iterator_category = typename BaseIterator::iterator_category;
        using value_type = typename BaseIterator::value_type;
        using difference_type = typename BaseIterator::difference_type;
        using pointer = typename BaseIterator::pointer;
        using reference = typename BaseIterator::reference;
    
        Iterator& operator++() {
        ++m_iterator;
            return *this;
        }
    
        bool operator==(const Iterator& other) const {
            return m_iterator == other.m_iterator;
        }
    
        bool operator!=(const Iterator& other) const { return !(*this == other); }
    
        Acts::SourceLink operator*() const {
            const auto& sl = m_iterator->second;
            return Acts::SourceLink{sl};
        }
    
        BaseIterator m_iterator;
    };
    
    // pointer to the container
    const Container* container = nullptr;
    
    // get the range of elements with requested key
    std::pair<Iterator, Iterator> range(const Acts::Surface& surface) const {
        assert(container != nullptr);
        auto [begin, end] = container->equal_range(surface.geometryId());
        return {Iterator{begin}, Iterator{end}};
    }
};

class BranchStopper {
    public:
        using BranchStopperResult =
            Acts::CombinatorialKalmanFilterBranchStopperResult;

        struct Config {
            int minMeasurements = 3;
        };

        BranchStopper(Config cfg) : m_cfg(std::move(cfg)) {};

        BranchStopperResult operator()(
            const TrackContainer::TrackProxy& track,
            const TrackContainer::TrackStateProxy& trackState) const {
        
            // bool enoughMeasurements =
                // track.nMeasurements() >= m_cfg.minMeasurements;
            // if (!enoughMeasurements) {
                // return BranchStopperResult::Continue;
            // }
            // else {
                // return BranchStopperResult::StopAndKeep;
            // }
            return BranchStopperResult::Continue;
        }
        
    private:
        const Config m_cfg;
};

BranchStopper branchStopper({2});

using StraightPropagator =
    Acts::Propagator<
        Acts::StraightLineStepper, 
        Acts::Experimental::DetectorNavigator>;
using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<
    ConstantFieldStepper, 
    Acts::Experimental::DetectorNavigator>;

using KalmanUpdater = Acts::GainMatrixUpdater;
using KalmanSmoother = Acts::GainMatrixSmoother;
using CombinatorialKalmanFilter =
    Acts::CombinatorialKalmanFilter<ConstantFieldPropagator, TrackContainer>;
using TestSourceLinkContainer =
    std::unordered_multimap<
        Acts::GeometryIdentifier, 
        Acts::detail::Test::TestSourceLink>;
using TestSourceLinkAccessor = TestContainerAccessor<TestSourceLinkContainer>;
using CombinatorialKalmanFilterOptions =
    Acts::CombinatorialKalmanFilterOptions<
        TestSourceLinkAccessor::Iterator,
        TrackContainer>;

KalmanUpdater kfUpdater;
KalmanSmoother kfSmoother;

Acts::GeometryContext gctx;
Acts::MagneticFieldContext mctx;
Acts::CalibrationContext cctx;

// configuration for the measurement selector
Acts::MeasurementSelector::Config measurementSelectorCfg = {
    // global default: no chi2 cut, only one measurement per surface
    {Acts::GeometryIdentifier(),
    {{}, {10000000000}, {1000u}}},
};

Acts::MeasurementSelector measSel{measurementSelectorCfg};

ConstantFieldPropagator makeConstantFieldPropagator(
    std::shared_ptr<const Acts::Experimental::Detector> geo, double bz) {
        Acts::Experimental::DetectorNavigator::Config cfg;
        cfg.detector = geo.get();
        cfg.resolvePassive = false;
        cfg.resolveMaterial = true;
        cfg.resolveSensitive = true;
        Acts::Experimental::DetectorNavigator navigator(
            cfg,
            Acts::getDefaultLogger("DetectorNavigator", Acts::Logging::VERBOSE));
        auto field =
            std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, bz));
        ConstantFieldStepper stepper(std::move(field));
        return ConstantFieldPropagator(std::move(stepper), std::move(navigator));
}

Acts::CombinatorialKalmanFilterExtensions<TrackContainer> 
    getExtensions() {
        Acts::CombinatorialKalmanFilterExtensions<TrackContainer> extensions;
            extensions.calibrator.template connect<
                &Acts::detail::Test::testSourceLinkCalibrator<TrackStateContainerBackend>>();
        extensions.updater.template connect<
            &KalmanUpdater::operator()<TrackStateContainerBackend>>(&kfUpdater);
        extensions.measurementSelector.template connect<
            &Acts::MeasurementSelector::select<TrackStateContainerBackend>>(
            &measSel);
        extensions.branchStopper.template connect<
            &BranchStopper::operator()>(&branchStopper);
        return extensions;
}

CombinatorialKalmanFilterOptions makeCkfOptions() {
    // leave the accessor empty, this will have to be set before running the CKF
    return CombinatorialKalmanFilterOptions(
        gctx, mctx, cctx,
        Acts::SourceLinkAccessorDelegate<TestSourceLinkAccessor::Iterator>{},
        getExtensions(), Acts::PropagatorPlainOptions(gctx, mctx));
}

// Intersection finding to get the
// region of interest for seeding
class NoFieldIntersectionFinder {
 public:
  ActsScalar m_tol = 1e-4;

  std::vector<const Surface*> m_surfaces;

  // Find the intersections along the path
  // and return them in the order of the path
  // length
  std::vector<std::pair<GeometryIdentifier, Vector3>> operator()(
      const GeometryContext& geoCtx, const Vector3& position,
      const Vector3& direction, [[maybe_unused]] const ActsScalar& Pmag = 0,
      [[maybe_unused]] const ActsScalar& Charge = 0) const {
    std::vector<std::pair<GeometryIdentifier, Vector3>> sIntersections;
    // Intersect the surfaces
    for (auto& surface : m_surfaces) {
      // Get the intersection
      auto sMultiIntersection = surface->intersect(
          geoCtx, position, direction,
          BoundaryTolerance::AbsoluteCartesian(m_tol, m_tol));

      // Take the closest
      auto closestForward = sMultiIntersection.closestForward();

      // Store if the intersection is reachable
      if (closestForward.status() == IntersectionStatus::reachable &&
          closestForward.pathLength() > 0.0) {
        sIntersections.push_back(
            {closestForward.object()->geometryId(), closestForward.position()});
        continue;
      }
    }
    return sIntersections;
  }
};

// A simple path width provider to set
// the grid lookup boundaries around the
// intersection point
class PathWidthProvider {
 public:
    std::pair<ActsScalar, ActsScalar> width;

  std::pair<ActsScalar, ActsScalar> operator()(
      const GeometryContext& /*gctx*/, const GeometryIdentifier& /*geoId*/) const {
    return width;
  }
};

// Calibrator to transform the source links
// to global coordinates
class SourceLinkCalibrator {
 public:
  SourceLinkSurfaceAccessor m_surfaceAccessor;

  Vector3 operator()(const GeometryContext& geoCtx,
                     const SourceLink& sourceLink) const {
    auto ssl = sourceLink.get<detail::Test::TestSourceLink>();
    auto res = m_surfaceAccessor(sourceLink)
                   ->localToGlobal(geoCtx, ssl.parameters, Vector3{0, 1, 0});
    return res;
  }
};

// Estimator of the particle's energy,
// vertex, momentum direction at the IP
// and the direction at the first hit
class TrackEstimator {
 public:
  Vector3 ip;

  std::tuple<ActsScalar, ActsScalar, Vector3, Vector3, Vector3> operator()(
      const GeometryContext& /*geoCtx*/, const Vector3& pivot) const {
    Vector3 direction = (pivot - ip).normalized();
    return {1_e, 1._GeV, ip, direction, direction};
  }
};

// Construct grid with the source links
struct ConstructSourceLinkGrid {
    SourceLinkSurfaceAccessor m_surfaceAccessor;
   
    std::unordered_map<GeometryIdentifier, Grid> construct(
        GeometryContext geoCtx,
        std::vector<SourceLink> sourceLinks) {
            // Lookup table for each layer
            std::unordered_map<GeometryIdentifier, Grid> lookupTable;
        
            // Construct a binned grid for each layer
            for (int i : {14, 15, 16, 17}) {
                Axis xAxis(-halfY, halfY, 50);
                Axis yAxis(-halfZ, halfZ, 50);
            
                Grid grid(std::make_tuple(xAxis, yAxis));

                GeometryIdentifier geoId;

                geoId.setSensitive(i);

                lookupTable.insert({geoId, grid});
            }
            // Fill the grid with source links
            for (auto& sl : sourceLinks) {
                auto ssl = sl.get<detail::Test::TestSourceLink>();
                auto id = ssl.m_geometryId;
            
                // Grid works with global positions
                Vector3 globalPos = m_surfaceAccessor(sl)->localToGlobal(
                    geoCtx, ssl.parameters, Vector3{0, 1, 0});
            
                auto bin =
                    lookupTable.at(id)
                        .localBinsFromPosition(Vector2(globalPos.y(), globalPos.z()));
                lookupTable.at(id).atLocalBins(bin).push_back(sl);
            }
        
            return lookupTable;
    }

};

// Construct a simple telescope detector
std::shared_ptr<Experimental::Detector> constructTelescopeDetector() {
  RotationMatrix3 rotation;
  double angle = 90_degree;
  Vector3 xPos(cos(angle), 0., sin(angle));
  Vector3 yPos(0., 1., 0.);
  Vector3 zPos(-sin(angle), 0., cos(angle));
  rotation.col(0) = xPos;
  rotation.col(1) = yPos;
  rotation.col(2) = zPos;

  // Create a bunch of surface bounds
  auto surf0bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  auto surf1bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  auto surf2bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  auto surf3bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  // Create a bunch of surfaces
  auto transform0 = Transform3::Identity();
  auto surf0 = Surface::makeShared<PlaneSurface>(
      transform0 * Transform3(rotation), std::move(surf0bounds));
  auto geoId0 = GeometryIdentifier();
  geoId0.setSensitive(1);

  auto transform1 =
      Transform3::Identity() * Translation3(Vector3(2 * deltaX, 0, 0));
  auto surf1 = Surface::makeShared<PlaneSurface>(
      transform1 * Transform3(rotation), std::move(surf1bounds));
  auto geoId1 = GeometryIdentifier();
  geoId1.setSensitive(2);

  auto transform2 =
      Transform3::Identity() * Translation3(Vector3(4 * deltaX, 0, 0));
  auto surf2 = Surface::makeShared<PlaneSurface>(
      transform2 * Transform3(rotation), std::move(surf2bounds));
  auto geoId2 = GeometryIdentifier();
  geoId2.setSensitive(3);

  auto transform3 =
      Transform3::Identity() * Translation3(Vector3(6 * deltaX, 0, 0));
  auto surf3 = Surface::makeShared<PlaneSurface>(
      transform3 * Transform3(rotation), std::move(surf3bounds));
  auto geoId3 = GeometryIdentifier();
  geoId3.setSensitive(4);

  // Create a bunch of volume bounds
  auto vol0bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  auto vol1bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  auto vol2bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  auto vol3bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  // Create a bunch of volumes
  auto vol0 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol0",
      transform0, std::move(vol0bounds), {surf0}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  auto vol1 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol1",
      transform1, std::move(vol1bounds), {surf1}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  auto vol2 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol2",
      transform2, std::move(vol2bounds), {surf2}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  auto vol3 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol3",
      transform3, std::move(vol3bounds), {surf3}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes = {
      vol0, vol1, vol2, vol3};

  // Connect the volumes
  auto portalContainer = Experimental::detail::CuboidalDetectorHelper::connect(
      gctx, volumes, BinningValue::binX, {}, Logging::VERBOSE);

  // Make sure that the geometry ids are
  // independent of the potential Id generation
  // changes
  int id = 1;

  // Volume ids
  for (auto& volume : volumes) {
    volume->assignGeometryId(id);
    id++;
  }
  // Intervolume portal ids
  for (auto& volume : volumes) {
    for (auto& port : volume->portalPtrs()) {
      if (port->surface().geometryId() == 0) {
        port->surface().assignGeometryId(id);
        id++;
      }
    }
  }
  // Surface ids
  for (auto& surf : {surf0, surf1, surf2, surf3}) {
    auto geoId = GeometryIdentifier();
    geoId.setSensitive(id);
    surf->assignGeometryId(geoId);
    id++;
  }

  auto detector = Experimental::Detector::makeShared(
      "TelescopeDetector", volumes, Experimental::tryRootVolumes());

  return detector;
}

std::vector<SourceLink> createSourceLinks(
    const GeometryContext& geoCtx, const Experimental::Detector& detector) {
  NoFieldIntersectionFinder intersectionFinder;

  for (auto volume : detector.volumes()) {
    for (auto surface : volume->surfaces()) {
      intersectionFinder.m_surfaces.push_back(surface);
    }
  }

  Vector3 vertex(-5., 0., 0.);
//   std::vector<ActsScalar> phis = {-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15};
  std::vector<ActsScalar> phis = {-0.05, 0.05};

  std::vector<SourceLink> sourceLinks;
  for (ActsScalar phi : phis) {
    Vector3 direction(cos(phi), sin(phi), 0.);

    auto intersections = intersectionFinder(geoCtx, vertex, direction);

    SquareMatrix2 cov = SquareMatrix2::Identity();

    for (auto& [id, refPoint] : intersections) {
      auto surf = *detector.sensitiveHierarchyMap().find(id);
      Vector2 val = surf->globalToLocal(geoCtx, refPoint, direction).value();

      detail::Test::TestSourceLink sourceLink(eBoundLoc0, eBoundLoc1, val, cov,
                                              id, id.value());

      SourceLink sl{sourceLink};
      sourceLinks.push_back(sl);
    }
  }

  return sourceLinks;
}

BOOST_AUTO_TEST_CASE(PathSeederZeroField) {
  using SurfaceAccessor =
      detail::Test::Experimental::TestSourceLinkSurfaceAccessor;

  // Create detector
  auto detector = constructTelescopeDetector();

  // Create source links
  auto sourceLinks = createSourceLinks(gctx, *detector);

  // Prepare the PathSeeder
  auto pathSeederCfg = Acts::Experimental::PathSeeder::Config();

  // Grid to bin the source links
  SurfaceAccessor surfaceAccessor{*detector};
  auto sourceLinkGridConstructor = ConstructSourceLinkGrid();
    sourceLinkGridConstructor.m_surfaceAccessor.connect<&SurfaceAccessor::operator()>(
        &surfaceAccessor);

    // Create the grid
    std::unordered_map<GeometryIdentifier, Grid> sourceLinkGrid =
        sourceLinkGridConstructor.construct(gctx, sourceLinks);

  // Estimator of the IP and first hit
  // parameters of the track
  TrackEstimator trackEstimator;
  trackEstimator.ip = Vector3(-5., 0., 0.);
  pathSeederCfg.trackEstimator.connect<&TrackEstimator::operator()>(
      &trackEstimator);

  // Transforms the source links to global coordinates
  SourceLinkCalibrator sourceLinkCalibrator;
  sourceLinkCalibrator.m_surfaceAccessor.connect<&SurfaceAccessor::operator()>(
      &surfaceAccessor);
  pathSeederCfg.sourceLinkCalibrator.connect<&SourceLinkCalibrator::operator()>(
      &sourceLinkCalibrator);

  // Intersection finder
  NoFieldIntersectionFinder intersectionFinder;
  for (auto volume : detector->volumes()) {
    for (auto surface : volume->surfaces()) {
      intersectionFinder.m_surfaces.push_back(surface);
    }
  }
  pathSeederCfg.intersectionFinder
      .connect<&NoFieldIntersectionFinder::operator()>(&intersectionFinder);

  // Path width provider
    PathWidthProvider pathWidthProvider;

  pathSeederCfg.pathWidthProvider.connect<&PathWidthProvider::operator()>(&pathWidthProvider);

    GeometryIdentifier geoId;
    geoId.setSensitive(14);
    pathSeederCfg.firstLayerIds.push_back(geoId);

  // Create the PathSeeder
  Acts::Experimental::PathSeeder pathSeeder(pathSeederCfg);

  // Get the seeds
    pathWidthProvider.width = {0.01, 0.01};

    std::vector<Acts::Experimental::PathSeeder::Seed> seeds;
    // SeedTreeContainer seeds;
    pathSeeder.getSeeds(gctx, sourceLinkGrid, seeds);

    // Check the seeds
    BOOST_CHECK_EQUAL(seeds.size(), 2);

    auto seed = seeds.at(0);
    TestSourceLinkContainer ckfSourceLinks;
    for (auto& sl : seed.sourceLinks) {
        auto ssl = sl.get<detail::Test::TestSourceLink>();
        ckfSourceLinks.insert({ssl.m_geometryId, ssl});
    }

    auto options = makeCkfOptions();
    options.propagatorPlainOptions.direction = Acts::Direction::Forward;
    
    TestSourceLinkAccessor slAccessor;
    slAccessor.container = &ckfSourceLinks;
    options.sourcelinkAccessor.connect<&TestSourceLinkAccessor::range>(
        &slAccessor);
    
    TrackContainer tc{
        Acts::VectorTrackContainer{},
        Acts::VectorMultiTrajectory{}};
    
    // Create IP covariance matrix from
    // reasonable standard deviations
    Acts::BoundVector ipStdDev;
    ipStdDev[Acts::eBoundLoc0] = 100_um;
    ipStdDev[Acts::eBoundLoc1] = 100_um;
    ipStdDev[Acts::eBoundTime] = 25_ns;
    ipStdDev[Acts::eBoundPhi] = 2_degree;
    ipStdDev[Acts::eBoundTheta] = 2_degree;
    ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
    Acts::BoundSquareMatrix ipCov =
        ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

    // CKF implementation to be tested
    CombinatorialKalmanFilter ckf(
        makeConstantFieldPropagator(detector, 0), 
        Acts::getDefaultLogger("CKF", Acts::Logging::VERBOSE));

    // run the CKF for all initial track states
    for (std::size_t trackId = 0u;  trackId < 1; ++trackId) {
        Acts::Vector4 mPos4 = {seed.ipVertex.x(), seed.ipVertex.y(), seed.ipVertex.z(), 0};

        Acts::ActsScalar p = seed.ipP; 
        Acts::ActsScalar theta = std::acos(seed.ipDir.z()/p);
        Acts::ActsScalar phi = std::atan2(seed.ipDir.y(), seed.ipDir.x());

        Acts::CurvilinearTrackParameters ipParameters(
            mPos4, phi, theta,
            1_e / p, ipCov, 
            Acts::ParticleHypothesis::electron());

        options.propagatorPlainOptions.maxSteps = 10000;

        auto res = ckf.findTracks(ipParameters, options, tc);
        if (!res.ok()) {
            BOOST_TEST_INFO(res.error() << " " << res.error().message());
        }
        BOOST_REQUIRE(res.ok());

        std::cout << "Tracks = " << tc.size() << std::endl;
        std::cout << "\n\n\n" << std::endl;
        for (std::size_t tid = 0u; tid < tc.size(); ++tid) {
            const auto track = tc.getTrack(tid);
        
            // check purity of first found track
            // find the number of hits not originating from the right track
            std::size_t numHits = 0u;
            std::size_t nummismatchedHits = 0u;
            for (const auto trackState : track.trackStatesReversed()) {
                const auto& sl = trackState.getUncalibratedSourceLink();
                const auto& ssl = sl.get<detail::Test::TestSourceLink>();
                std::cout << "Track state: " 
                << trackState.parameters().transpose() << " ----> "
                << ssl.parameters.transpose()
                << std::endl;
            }
            std::cout << "----------------" << std::endl;
        }

    }
}

BOOST_AUTO_TEST_SUITE_END()
