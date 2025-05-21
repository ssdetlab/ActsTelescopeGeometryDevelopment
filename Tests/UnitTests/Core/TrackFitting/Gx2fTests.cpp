// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/detail/CuboidalDetectorHelper.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Visualization/EventDataView3D.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

#include <vector>

#include "FitterTestsCommon.hpp"

using namespace Acts::UnitLiterals;
using namespace Acts::detail::Test;

Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;
const auto gx2fLogger = Acts::getDefaultLogger("Gx2f", logLevel);

namespace Acts::Test {

/// @brief Helper function to visualise measurements in a 3D environment.
///
/// This function iterates through the provided measurements and visualises each
/// one using the specified 3D visualisation helper. The visualisation takes
/// into account the surface transformations and localisation errors.
///
/// @param helper The 3D visualisation helper used to draw the measurements.
/// @param measurements A collection of measurements to be visualised, containing source links with parameters and covariance information.
/// @param geometry A shared pointer to the constant tracking geometry used to find surfaces associated with measurements.
/// @param geoCtx The geometry context used for transformations and accessing geometry-related information.
/// @param locErrorScale Scaling factor for localisation errors. Default value is 1.0.
/// @param viewConfig Configuration settings for the visualisation view. Default value is s_viewMeasurement.
static void drawMeasurements(
    IVisualization3D& helper, const Measurements& measurements,
    const std::shared_ptr<const TrackingGeometry>& geometry,
    const Acts::GeometryContext& geoCtx, double locErrorScale = 1.,
    const ViewConfig& viewConfig = s_viewMeasurement) {
  std::cout << "\n*** Draw measurements ***\n" << std::endl;

  for (auto& singleMeasurement : measurements.sourceLinks) {
    auto cov = singleMeasurement.covariance;
    auto lposition = singleMeasurement.parameters;

    auto surf = geometry->findSurface(singleMeasurement.m_geometryId);
    auto transf = surf->transform(geoCtx);

    EventDataView3D::drawMeasurement(helper, lposition, cov, transf,
                                     locErrorScale, viewConfig);
  }
}

//// Construct initial track parameters.
Acts::CurvilinearTrackParameters makeParameters(
    const ActsScalar x = 0.0_m, const ActsScalar y = 0.0_m,
    const ActsScalar z = 0.0_m, const ActsScalar w = 42_ns,
    const ActsScalar phi = 0_degree, const ActsScalar theta = 90_degree,
    const ActsScalar p = 2_GeV, const ActsScalar q = 1_e) {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  const Acts::BoundSquareMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // define a track in the transverse plane along x
  const Acts::Vector4 mPos4(x, y, z, w);
  return Acts::CurvilinearTrackParameters(mPos4, phi, theta, q / p, cov,
                                          Acts::ParticleHypothesis::pion());
}

static std::vector<Acts::SourceLink> prepareSourceLinks(
    const std::vector<TestSourceLink>& sourceLinks) {
  std::vector<Acts::SourceLink> result;
  std::transform(sourceLinks.begin(), sourceLinks.end(),
                 std::back_inserter(result),
                 [](const auto& sl) { return Acts::SourceLink{sl}; });
  return result;
}

/// @brief Create a simple telescope detector in the Y direction.
///
/// We cannot reuse the previous detector, since the cuboid volume builder only
/// allows merging of YZ-faces.
///
/// @param geoCtx
/// @param nSurfaces Number of surfaces
std::shared_ptr<const Acts::Experimental::Detector> makeToyDetectorYdirection(
    const Acts::GeometryContext& geoContext) {
  // Construct a cubic detector with 3 volumes
  Acts::RotationMatrix3 rotation;
  double angle = 90_degree;
  Acts::Vector3 xPos(cos(angle), 0., sin(angle));
  Acts::Vector3 yPos(0., 1., 0.);
  Acts::Vector3 zPos(-sin(angle), 0., cos(angle));
  rotation.col(0) = xPos;
  rotation.col(1) = yPos;
  rotation.col(2) = zPos;

  auto bounds1 = std::make_unique<Acts::CuboidVolumeBounds>(3, 3, 3);
  auto transform1 = Acts::Transform3::Identity();
  auto surface1 = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform1 * Acts::Transform3(rotation),
      std::make_shared<Acts::RectangleBounds>(2, 2));
  auto volume1 = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), geoContext,
      "volume1", transform1, std::move(bounds1), {surface1}, {},
      Acts::Experimental::tryNoVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  auto bounds2 = std::make_unique<Acts::CuboidVolumeBounds>(3, 3, 3);
  auto transform2 =
      Acts::Transform3::Identity() * Acts::Translation3(Acts::Vector3(6, 0, 0));
  auto surface2 = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform2 * Acts::Transform3(rotation),
      std::make_shared<Acts::RectangleBounds>(2, 2));
  auto volume2 = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), geoContext,
      "volume2", transform2, std::move(bounds2), {surface2}, {},
      Acts::Experimental::tryNoVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  auto bounds3 = std::make_unique<Acts::CuboidVolumeBounds>(3, 3, 3);
  auto transform3 = Acts::Transform3::Identity() *
                    Acts::Translation3(Acts::Vector3(12, 0, 0));
  auto surface3 = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform3 * Acts::Transform3(rotation),
      std::make_shared<Acts::RectangleBounds>(2, 2));
  auto volume3 = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), geoContext,
      "volume3", transform3, std::move(bounds3), {surface3}, {},
      Acts::Experimental::tryNoVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      detectorVolumes = {volume1, volume2, volume3};

  auto portalContainer =
      Acts::Experimental::detail::CuboidalDetectorHelper::connect(
          geoContext, detectorVolumes, Acts::BinningValue::binX, {},
          Acts::Logging::VERBOSE);

  // Make sure that the geometry ids are
  // independent of the potential Id generation
  // changes
  int id = 1;

  // Volume ids: 1-3
  for (auto& volume : detectorVolumes) {
    volume->assignGeometryId(id);
    id++;
  }
  // Intervolume portal ids: 6,7,10,11
  for (auto& volume : detectorVolumes) {
    for (auto& port : volume->portalPtrs()) {
      if (port->surface().geometryId() == 0) {
        port->surface().assignGeometryId(id);
        id++;
      }
    }
  }
  // Surface ids: 12-14
  for (auto& surf : {surface1, surface2, surface3}) {
    Acts::GeometryIdentifier geoId;
    geoId.setSensitive(id);
    surf->assignGeometryId(geoId);
    id++;
  }

  auto detector = Acts::Experimental::Detector::makeShared(
      "cubicDetector", detectorVolumes, Acts::Experimental::tryRootVolumes());
  return detector;
}

// Construct a straight-line propagator.
auto makeStraightPropagator(
    std::shared_ptr<const Acts::Experimental::Detector> geo) {
  Acts::Experimental::DetectorNavigator::Config cfg;
  cfg.detector = geo.get();
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Experimental::DetectorNavigator navigator(
      cfg, Acts::getDefaultLogger("Navigator", Acts::Logging::VERBOSE));
  Acts::StraightLineStepper stepper;
  return Acts::Propagator<Acts::StraightLineStepper,
                          Acts::Experimental::DetectorNavigator>(
      stepper, std::move(navigator));
}

template <typename stepper_t>
auto makeConstantFieldPropagator(
    std::shared_ptr<const Acts::Experimental::Detector> geo, double bz) {
  Acts::Experimental::DetectorNavigator::Config cfg;
  cfg.detector = geo.get();
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Experimental::DetectorNavigator navigator(
      cfg, Acts::getDefaultLogger("Navigator", Acts::Logging::VERBOSE));
  auto field =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, bz));
  stepper_t stepper(std::move(field));
  return Acts::Propagator<decltype(stepper),
                          Acts::Experimental::DetectorNavigator>(
      std::move(stepper), std::move(navigator));
}

struct Detector {
  // geometry
  std::shared_ptr<const Acts::Experimental::Detector> geometry;
};

BOOST_AUTO_TEST_SUITE(Gx2fTest)
ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("Gx2fTests", logLevel))

// Context objects
const Acts::GeometryContext geoCtx;
const Acts::MagneticFieldContext magCtx;
const Acts::CalibrationContext calCtx;

// Measurement resolutions
const MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                        {25_um, 50_um}};
const MeasurementResolutionMap resMapAllPixel = {
    {Acts::GeometryIdentifier().setVolume(0), resPixel}};

// This test checks if the call to the fitter works and no errors occur in the
// framework, without fitting and updating any parameters
//
BOOST_AUTO_TEST_CASE(Fit5Iterations) {
  ACTS_INFO("*** Test: Fit5Iterations -- Start");

  std::default_random_engine rng(42);

  ACTS_DEBUG("Create the detector");
  const std::size_t nSurfaces = 3;
  Detector detector;
  detector.geometry = makeToyDetectorYdirection(geoCtx);

  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
  const auto parametersMeasurements = makeParameters(-3_mm);
  const auto startParametersFit = makeParameters(-3_mm, 2.5_mm, 2.5_mm);

  ACTS_DEBUG("Create the measurements");
  using SimPropagator = Acts::Propagator<Acts::StraightLineStepper,
                                         Acts::Experimental::DetectorNavigator>;
  const SimPropagator simPropagator = makeStraightPropagator(detector.geometry);
  const auto measurements =
      createMeasurements(simPropagator, geoCtx, magCtx, parametersMeasurements,
                         resMapAllPixel, rng);
  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());

  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);

  ACTS_DEBUG("Set up the fitter");
  const Surface* rSurface = &parametersMeasurements.referenceSurface();

  using RecoStepper = EigenStepper<>;
  const auto recoPropagator =
      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);

  using RecoPropagator = decltype(recoPropagator);
  using Gx2Fitter =
      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());

  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  Acts::detail::Test::Experimental::TestSourceLinkSurfaceAccessor
      surfaceAccessor{*detector.geometry};
  extensions.surfaceAccessor
      .connect<&Acts::detail::Test::Experimental::
                   TestSourceLinkSurfaceAccessor::operator()>(&surfaceAccessor);

  const Experimental::Gx2FitterOptions gx2fOptions(
      geoCtx, magCtx, calCtx, extensions,
      PropagatorPlainOptions(geoCtx, magCtx), rSurface, false, false,
      FreeToBoundCorrection(false), 5, 0);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  ACTS_DEBUG("Fit the track");
  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
                              startParametersFit, gx2fOptions, tracks);

  BOOST_REQUIRE(res.ok());

  const auto& track = *res;

  BOOST_CHECK_EQUAL(track.tipIndex(), nSurfaces - 1);
  BOOST_CHECK(track.hasReferenceSurface());

  // Track quantities
  CHECK_CLOSE_ABS(track.chi2(), 8., 2.);
  BOOST_CHECK_EQUAL(track.nDoF(), nSurfaces * 2);
  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
  BOOST_CHECK_EQUAL(track.nMeasurements(), nSurfaces);
  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);

  // Parameters
  // We need quite coarse checks here, since on different builds
  // the created measurements differ in the randomness
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -11., 7e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -15., 6e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1e-5, 1e3);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], M_PI / 2, 1e-3);
  BOOST_CHECK_EQUAL(track.parameters()[eBoundQOverP], 1);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime],
                    startParametersFit.parameters()[eBoundTime], 1e-6);
  BOOST_CHECK_CLOSE(track.covariance().determinant(), 1e-27, 4e0);

  // Convergence
  BOOST_CHECK_EQUAL(
      (track.template component<
          std::uint32_t,
          hashString(Experimental::Gx2fConstants::gx2fnUpdateColumn)>()),
      5);

  ACTS_INFO("*** Test: Fit5Iterations -- Finish");
}

// BOOST_AUTO_TEST_CASE(UpdatePushedToNewVolume) {
//   ACTS_INFO("*** Test: UpdatePushedToNewVolume -- Start");
//
//   std::default_random_engine rng(42);
//
//   ACTS_DEBUG("Create the detector");
//   const std::size_t nSurfaces = 5;
//   Detector detector;
//   detector.geometry = makeToyDetectorYdirection(geoCtx, nSurfaces);
//
//   ACTS_DEBUG("Set the start parameters for measurement creation and fit");
//   const auto parametersMeasurements =
//       makeParameters(0_mm, 0_mm, 0_mm, 42_ns, 90_degree, 90_degree, 1_GeV,
//       1_e);
//   const auto startParametersFit = makeParameters(
//       1500_mm, 0_mm, 0_mm, 42_ns, -45_degree, -90_degree, 1_GeV, 1_e);
//
//   ACTS_DEBUG("Create the measurements");
//   using SimPropagator = Acts::Propagator<Acts::StraightLineStepper,
//                                          Acts::Experimental::DetectorNavigator>;
//   const SimPropagator simPropagator =
//   makeStraightPropagator(detector.geometry); const auto measurements =
//       createMeasurements(simPropagator, geoCtx, magCtx,
//       parametersMeasurements,
//                          resMapAllPixel, rng);
//   const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
//   ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());
//
//   BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);
//
//   ACTS_DEBUG("Set up the fitter");
//   const Surface* rSurface = &parametersMeasurements.referenceSurface();
//
//   using RecoStepper = EigenStepper<>;
//   const auto recoPropagator =
//       makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);
//
//   using RecoPropagator = decltype(recoPropagator);
//   using Gx2Fitter =
//       Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
//   const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());
//
//   Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
//   extensions.calibrator
//       .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
//   TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
//   extensions.surfaceAccessor
//       .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);
//
//   const Experimental::Gx2FitterOptions gx2fOptions(
//       geoCtx, magCtx, calCtx, extensions,
//       PropagatorPlainOptions(geoCtx, magCtx), rSurface, false, false,
//       FreeToBoundCorrection(false), 5, 0);
//
//   Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
//                               Acts::VectorMultiTrajectory{}};
//
//   ACTS_DEBUG("Fit the track");
//  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
//  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
//  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
//                              startParametersFit, gx2fOptions, tracks);
//
//  // Helper to visualise the detector
//  {
//    std::cout << "\n*** Create .obj of Detector ***\n" << std::endl;
//    // Only need for obj
//    ObjVisualization3D obj;
//
//    bool triangulate = true;
//    ViewConfig viewSensitive = ViewConfig({0, 180, 240});
//    viewSensitive.triangulate = triangulate;
//    ViewConfig viewPassive = ViewConfig({240, 280, 0});
//    viewPassive.triangulate = triangulate;
//    ViewConfig viewVolume = ViewConfig({220, 220, 0});
//    viewVolume.triangulate = triangulate;
//    ViewConfig viewContainer = ViewConfig({220, 220, 0});
//    viewContainer.triangulate = triangulate;
//    ViewConfig viewGrid = ViewConfig({220, 0, 0});
//    viewGrid.nSegments = 8;
//    viewGrid.offset = 3.;
//    viewGrid.triangulate = triangulate;
//
//    std::string tag = "gx2f_toydet";
//
//    const Acts::TrackingVolume& tgVolume =
//        *(detector.geometry->highestTrackingVolume());
//
//    GeometryView3D::drawTrackingVolume(obj, tgVolume, geoCtx, viewContainer,
//                                       viewVolume, viewPassive, viewSensitive,
//                                       viewGrid, true, tag);
//  }
//  // Helper to visualise the measurements
//  {
//    std::cout << "\n*** Create .obj of measurements ***\n" << std::endl;
//    ObjVisualization3D obj;
//
//    double localErrorScale = 10000000.;
//    ViewConfig mcolor({255, 145, 48});
//    mcolor.offset = 2;
//    //  mcolor.visible = true;
//
//    drawMeasurements(obj, measurements, detector.geometry, geoCtx,
//                     localErrorScale, mcolor);
//
//    obj.write("meas");
//  }
//
//  BOOST_REQUIRE(!res.ok());
//  BOOST_CHECK_EQUAL(
//      res.error(),
//      Acts::Experimental::GlobalChiSquareFitterError::UpdatePushedToNewVolume);
//
//  ACTS_INFO("*** Test: UpdatePushedToNewVolume -- Finish");
//}
//
// BOOST_AUTO_TEST_CASE(MixedDetector) {
//  ACTS_INFO("*** Test: MixedDetector -- Start");
//
//  std::default_random_engine rng(42);
//
//  ACTS_DEBUG("Create the detector");
//  const std::size_t nSurfaces = 7;
//  Detector detector;
//  detector.geometry = makeToyDetector(geoCtx, nSurfaces);
//
//  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
//  const auto parametersMeasurements = makeParameters();
//  const auto startParametersFit = makeParameters(
//      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);
//
//  ACTS_DEBUG("Create the measurements");
//  const MeasurementResolutionMap resMap = {
//      {Acts::GeometryIdentifier().setVolume(2).setLayer(2), resPixel},
//      {Acts::GeometryIdentifier().setVolume(2).setLayer(4), resStrip0},
//      {Acts::GeometryIdentifier().setVolume(2).setLayer(6), resStrip1},
//      {Acts::GeometryIdentifier().setVolume(2).setLayer(8), resPixel},
//      {Acts::GeometryIdentifier().setVolume(2).setLayer(10), resStrip0},
//      {Acts::GeometryIdentifier().setVolume(2).setLayer(12), resStrip1},
//      {Acts::GeometryIdentifier().setVolume(2).setLayer(14), resPixel},
//  };
//
//  using SimPropagator = Acts::Propagator<Acts::StraightLineStepper,
//                                         Acts::Experimental::DetectorNavigator>;
//  const SimPropagator simPropagator =
//  makeStraightPropagator(detector.geometry); const auto measurements =
//  createMeasurements(
//      simPropagator, geoCtx, magCtx, parametersMeasurements, resMap, rng);
//  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
//  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());
//
//  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);
//
//  ACTS_DEBUG("Set up the fitter");
//  const Surface* rSurface = &parametersMeasurements.referenceSurface();
//
//  using RecoStepper = EigenStepper<>;
//  const auto recoPropagator =
//      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);
//
//  using RecoPropagator = decltype(recoPropagator);
//  using Gx2Fitter =
//      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
//  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());
//
//  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
//  extensions.calibrator
//      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
//  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
//  extensions.surfaceAccessor
//      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);
//
//  const Experimental::Gx2FitterOptions gx2fOptions(
//      geoCtx, magCtx, calCtx, extensions,
//      PropagatorPlainOptions(geoCtx, magCtx), rSurface, false, false,
//      FreeToBoundCorrection(false), 5, 0);
//
//  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
//                              Acts::VectorMultiTrajectory{}};
//
//  ACTS_DEBUG("Fit the track");
//  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
//  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
//  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
//                              startParametersFit, gx2fOptions, tracks);
//
//  BOOST_REQUIRE(res.ok());
//
//  const auto& track = *res;
//  BOOST_CHECK_EQUAL(track.tipIndex(), nSurfaces - 1);
//  BOOST_CHECK(track.hasReferenceSurface());
//
//  // Track quantities
//  CHECK_CLOSE_ABS(track.chi2(), 8.5, 4.);
//  BOOST_CHECK_EQUAL(track.nDoF(), 10u);
//  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
//  BOOST_CHECK_EQUAL(track.nMeasurements(), nSurfaces);
//  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
//  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);
//
//  // Parameters
//  // We need quite coarse checks here, since on different builds
//  // the created measurements differ in the randomness
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -11., 7e0);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -15., 6e0);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1e-5, 1e3);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], M_PI / 2, 1e-3);
//  BOOST_CHECK_EQUAL(track.parameters()[eBoundQOverP], 1);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime],
//                    startParametersFit.parameters()[eBoundTime], 1e-6);
//  BOOST_CHECK_CLOSE(track.covariance().determinant(), 2e-28, 1e0);
//
//  // Convergence
//  BOOST_CHECK_EQUAL(
//      (track.template component<
//          std::uint32_t,
//          hashString(Experimental::Gx2fConstants::gx2fnUpdateColumn)>()),
//      5);
//
//  ACTS_INFO("*** Test: MixedDetector -- Finish");
//}
//
//// This test checks if we can fit QOverP, when a magnetic field is introduced
// BOOST_AUTO_TEST_CASE(FitWithBfield) {
//   ACTS_INFO("*** Test: FitWithBfield -- Start");
//
//   std::default_random_engine rng(42);
//
//   ACTS_DEBUG("Create the detector");
//   const std::size_t nSurfaces = 5;
//   Detector detector;
//   detector.geometry = makeToyDetector(geoCtx, nSurfaces);
//
//   ACTS_DEBUG("Set the start parameters for measurement creation and fit");
//   const auto parametersMeasurements = makeParameters();
//   const auto startParametersFit = makeParameters(
//       7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);
//
//   ACTS_DEBUG("Create the measurements");
//   using SimStepper = EigenStepper<>;
//   const auto simPropagator =
//       makeConstantFieldPropagator<SimStepper>(detector.geometry, 0.3_T);
//
//   const auto measurements =
//       createMeasurements(simPropagator, geoCtx, magCtx,
//       parametersMeasurements,
//                          resMapAllPixel, rng);
//
//   const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
//   ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());
//
//   BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);
//
//   ACTS_DEBUG("Set up the fitter");
//   const Surface* rSurface = &parametersMeasurements.referenceSurface();
//
//   // Reuse the SimPropagator, since it already uses the EigenStepper<>
//   using SimPropagator = decltype(simPropagator);
//   using Gx2Fitter =
//       Experimental::Gx2Fitter<SimPropagator, VectorMultiTrajectory>;
//   const Gx2Fitter fitter(simPropagator, gx2fLogger->clone());
//
//   Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
//   extensions.calibrator
//       .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
//   TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
//   extensions.surfaceAccessor
//       .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);
//
//   const Experimental::Gx2FitterOptions gx2fOptions(
//       geoCtx, magCtx, calCtx, extensions,
//       PropagatorPlainOptions(geoCtx, magCtx), rSurface, false, false,
//       FreeToBoundCorrection(false), 5, 0);
//
//   Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
//                               Acts::VectorMultiTrajectory{}};
//
//   ACTS_DEBUG("Fit the track");
//  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
//  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
//  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
//                              startParametersFit, gx2fOptions, tracks);
//
//  BOOST_REQUIRE(res.ok());
//
//  const auto& track = *res;
//  BOOST_CHECK_EQUAL(track.tipIndex(), nSurfaces - 1);
//  BOOST_CHECK(track.hasReferenceSurface());
//
//  // Track quantities
//  CHECK_CLOSE_ABS(track.chi2(), 7.5, 1.5);
//  BOOST_CHECK_EQUAL(track.nDoF(), nSurfaces * 2);
//  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
//  BOOST_CHECK_EQUAL(track.nMeasurements(), nSurfaces);
//  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
//  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);
//
//  // Parameters
//  // We need quite coarse checks here, since on different builds
//  // the created measurements differ in the randomness
//  // TODO investigate further the reference values for eBoundPhi and det(cov)
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -11., 8e0);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -15., 6e0);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1e-4, 1e3);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], M_PI / 2, 1e-3);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundQOverP], 0.5, 2e-1);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime],
//                    startParametersFit.parameters()[eBoundTime], 1e-6);
//  BOOST_CHECK_CLOSE(track.covariance().determinant(), 8e-35, 4e0);
//
//  // Convergence
//  BOOST_CHECK_EQUAL(
//      (track.template component<
//          std::uint32_t,
//          hashString(Experimental::Gx2fConstants::gx2fnUpdateColumn)>()),
//      5);
//
//  ACTS_INFO("*** Test: FitWithBfield -- Finish");
//}
//
// BOOST_AUTO_TEST_CASE(relChi2changeCutOff) {
//  ACTS_INFO("*** Test: relChi2changeCutOff -- Start");
//
//  std::default_random_engine rng(42);
//
//  ACTS_DEBUG("Create the detector");
//  const std::size_t nSurfaces = 5;
//  Detector detector;
//  detector.geometry = makeToyDetector(geoCtx, nSurfaces);
//
//  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
//  const auto parametersMeasurements = makeParameters();
//  const auto startParametersFit = makeParameters(
//      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);
//
//  ACTS_DEBUG("Create the measurements");
//  // simulation propagator
//  using SimPropagator = Acts::Propagator<Acts::StraightLineStepper,
//                                         Acts::Experimental::DetectorNavigator>;
//  const SimPropagator simPropagator =
//  makeStraightPropagator(detector.geometry); const auto measurements =
//      createMeasurements(simPropagator, geoCtx, magCtx,
//      parametersMeasurements,
//                         resMapAllPixel, rng);
//  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
//  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());
//
//  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);
//
//  ACTS_DEBUG("Set up the fitter");
//  const Surface* rSurface = &parametersMeasurements.referenceSurface();
//
//  using RecoStepper = EigenStepper<>;
//  const auto recoPropagator =
//      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);
//
//  using RecoPropagator = decltype(recoPropagator);
//  using Gx2Fitter =
//      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
//  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());
//
//  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
//  extensions.calibrator
//      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
//  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
//  extensions.surfaceAccessor
//      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);
//
//  const Experimental::Gx2FitterOptions gx2fOptions(
//      geoCtx, magCtx, calCtx, extensions,
//      PropagatorPlainOptions(geoCtx, magCtx), rSurface, false, false,
//      FreeToBoundCorrection(false), 500, 1e-5);
//
//  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
//                              Acts::VectorMultiTrajectory{}};
//
//  ACTS_DEBUG("Fit the track");
//  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
//  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
//  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
//                              startParametersFit, gx2fOptions, tracks);
//
//  BOOST_REQUIRE(res.ok());
//
//  const auto& track = *res;
//  BOOST_CHECK_EQUAL(track.tipIndex(), nSurfaces - 1);
//  BOOST_CHECK(track.hasReferenceSurface());
//
//  // Track quantities
//  CHECK_CLOSE_ABS(track.chi2(), 8., 2.);
//  BOOST_CHECK_EQUAL(track.nDoF(), nSurfaces * 2);
//  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
//  BOOST_CHECK_EQUAL(track.nMeasurements(), nSurfaces);
//  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
//  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);
//
//  // Parameters
//  // We need quite coarse checks here, since on different builds
//  // the created measurements differ in the randomness
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -11., 7e0);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -15., 6e0);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1e-5, 1e3);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], M_PI / 2, 1e-3);
//  BOOST_CHECK_EQUAL(track.parameters()[eBoundQOverP], 1);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime],
//                    startParametersFit.parameters()[eBoundTime], 1e-6);
//  BOOST_CHECK_CLOSE(track.covariance().determinant(), 1e-27, 4e0);
//
//  // Convergence
//  BOOST_CHECK_LT(
//      (track.template component<
//          std::uint32_t,
//          hashString(Experimental::Gx2fConstants::gx2fnUpdateColumn)>()),
//      10);
//
//  ACTS_INFO("*** Test: relChi2changeCutOff -- Finish");
//}
//
// BOOST_AUTO_TEST_CASE(DidNotConverge) {
//  ACTS_INFO("*** Test: DidNotConverge -- Start");
//
//  std::default_random_engine rng(42);
//
//  ACTS_DEBUG("Create the detector");
//  const std::size_t nSurfaces = 5;
//  Detector detector;
//  detector.geometry = makeToyDetector(geoCtx, nSurfaces);
//
//  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
//  const auto parametersMeasurements = makeParameters();
//  const auto startParametersFit = makeParameters(
//      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);
//
//  ACTS_DEBUG("Create the measurements");
//  // simulation propagator
//  using SimPropagator = Acts::Propagator<Acts::StraightLineStepper,
//                                         Acts::Experimental::DetectorNavigator>;
//  const SimPropagator simPropagator =
//  makeStraightPropagator(detector.geometry); const auto measurements =
//      createMeasurements(simPropagator, geoCtx, magCtx,
//      parametersMeasurements,
//                         resMapAllPixel, rng);
//  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
//  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());
//
//  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);
//
//  ACTS_DEBUG("Set up the fitter");
//  const Surface* rSurface = &parametersMeasurements.referenceSurface();
//
//  using RecoStepper = EigenStepper<>;
//  const auto recoPropagator =
//      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);
//
//  using RecoPropagator = decltype(recoPropagator);
//  using Gx2Fitter =
//      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
//  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());
//
//  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
//  extensions.calibrator
//      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
//  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
//  extensions.surfaceAccessor
//      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);
//
//  // The relChi2changeCutOff = 0 prevents to stop the fitter after
//  convergence,
//  // therefore all updates will be done (even if the result does not change).
//  // Since we didn't break due to convergence, we reach nUpdatesMax and
//  // therefore fail the fit.
//  const Experimental::Gx2FitterOptions gx2fOptions(
//      geoCtx, magCtx, calCtx, extensions,
//      PropagatorPlainOptions(geoCtx, magCtx), rSurface, false, false,
//      FreeToBoundCorrection(false), 6, 0);
//
//  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
//                              Acts::VectorMultiTrajectory{}};
//
//  ACTS_DEBUG("Fit the track");
//  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
//  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
//  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
//                              startParametersFit, gx2fOptions, tracks);
//
//  BOOST_REQUIRE(!res.ok());
//  BOOST_CHECK_EQUAL(
//      res.error(),
//      Acts::Experimental::GlobalChiSquareFitterError::DidNotConverge);
//
//  ACTS_INFO("*** Test: DidNotConverge -- Finish");
//}
//
// BOOST_AUTO_TEST_CASE(NotEnoughMeasurements) {
//  ACTS_INFO("*** Test: NotEnoughMeasurements -- Start");
//
//  std::default_random_engine rng(42);
//
//  ACTS_DEBUG("Create the detector");
//  const std::size_t nSurfaces = 2;
//  Detector detector;
//  detector.geometry = makeToyDetector(geoCtx, nSurfaces);
//
//  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
//  const auto parametersMeasurements = makeParameters();
//  const auto startParametersFit = makeParameters(
//      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);
//
//  ACTS_DEBUG("Create the measurements");
//  // simulation propagator
//  using SimPropagator = Acts::Propagator<Acts::StraightLineStepper,
//                                         Acts::Experimental::DetectorNavigator>;
//  const SimPropagator simPropagator =
//  makeStraightPropagator(detector.geometry); const auto measurements =
//      createMeasurements(simPropagator, geoCtx, magCtx,
//      parametersMeasurements,
//                         resMapAllPixel, rng);
//  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
//  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());
//
//  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);
//
//  ACTS_DEBUG("Set up the fitter");
//  const Surface* rSurface = &parametersMeasurements.referenceSurface();
//
//  using RecoStepper = EigenStepper<>;
//  const auto recoPropagator =
//      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);
//
//  using RecoPropagator = decltype(recoPropagator);
//  using Gx2Fitter =
//      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
//  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());
//
//  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
//  extensions.calibrator
//      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
//  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
//  extensions.surfaceAccessor
//      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);
//
//  const Experimental::Gx2FitterOptions gx2fOptions(
//      geoCtx, magCtx, calCtx, extensions,
//      PropagatorPlainOptions(geoCtx, magCtx), rSurface, false, false,
//      FreeToBoundCorrection(false), 6, 0);
//
//  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
//                              Acts::VectorMultiTrajectory{}};
//
//  ACTS_DEBUG("Fit the track");
//  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
//  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
//  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
//                              startParametersFit, gx2fOptions, tracks);
//
//  BOOST_REQUIRE(!res.ok());
//  BOOST_CHECK_EQUAL(
//      res.error(),
//      Acts::Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements);
//
//  ACTS_INFO("*** Test: NotEnoughMeasurements -- Finish");
//}
//
// BOOST_AUTO_TEST_CASE(FindHoles) {
//  ACTS_INFO("*** Test: FindHoles -- Start");
//
//  std::default_random_engine rng(42);
//
//  ACTS_DEBUG("Create the detector");
//  //  const std::size_t nSurfaces = 7;
//  const std::size_t nSurfaces = 8;
//  Detector detector;
//  detector.geometry = makeToyDetector(geoCtx, nSurfaces);
//
//  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
//  const auto parametersMeasurements = makeParameters();
//  const auto startParametersFit = makeParameters(
//      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);
//
//  ACTS_DEBUG("Create the measurements");
//  using SimPropagator = Acts::Propagator<Acts::StraightLineStepper,
//                                         Acts::Experimental::DetectorNavigator>;
//  const SimPropagator simPropagator =
//  makeStraightPropagator(detector.geometry); const auto measurements =
//      createMeasurements(simPropagator, geoCtx, magCtx,
//      parametersMeasurements,
//                         resMapAllPixel, rng);
//  auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
//  ACTS_VERBOSE("sourceLinks.size() [before] = " << sourceLinks.size());
//
//  // We remove the first measurement in the list. This does not create a hole.
//  sourceLinks.erase(std::next(sourceLinks.begin(), 0));
//  ACTS_VERBOSE(
//      "sourceLinks.size() [after first erase] = " << sourceLinks.size());
//
//  // We remove the last measurement in the list. This does not create a hole.
//  sourceLinks.pop_back();
//  ACTS_VERBOSE("sourceLinks.size() [after pop] = " << sourceLinks.size());
//
//  // We remove the second to last measurement in the list. This effectively
//  // creates a hole on that surface.
//  const std::size_t indexHole = sourceLinks.size() - 2;
//  ACTS_VERBOSE("Remove measurement " << indexHole);
//  sourceLinks.erase(std::next(sourceLinks.begin(), indexHole));
//  ACTS_VERBOSE("sourceLinks.size() [after second-to-last erase]= "
//               << sourceLinks.size());
//
//  // We removed 3 measurements
//  //  const std::size_t nMeasurements = nSurfaces - 2;
//  const std::size_t nMeasurements = nSurfaces - 3;
//  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);
//
//  ACTS_DEBUG("Set up the fitter");
//  const Surface* rSurface = &parametersMeasurements.referenceSurface();
//
//  using RecoStepper = EigenStepper<>;
//  const auto recoPropagator =
//      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);
//
//  using RecoPropagator = decltype(recoPropagator);
//  using Gx2Fitter =
//      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
//  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());
//
//  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
//  extensions.calibrator
//      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
//  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
//  extensions.surfaceAccessor
//      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);
//
//  const Experimental::Gx2FitterOptions gx2fOptions(
//      geoCtx, magCtx, calCtx, extensions,
//      PropagatorPlainOptions(geoCtx, magCtx), rSurface, false, false,
//      FreeToBoundCorrection(false), 20, 1e-5);
//
//  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
//                              Acts::VectorMultiTrajectory{}};
//
//  ACTS_DEBUG("Fit the track");
//  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
//  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
//  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
//                              startParametersFit, gx2fOptions, tracks);
//
//  BOOST_REQUIRE(res.ok());
//
//  const auto& track = *res;
//
//  // -1, because the index starts at 0
//  // -2, because the first and the last surface are not part of the track
//  BOOST_CHECK_EQUAL(track.tipIndex(), nSurfaces - 1 - 2);
//  BOOST_CHECK(track.hasReferenceSurface());
//
//  // Track quantities
//  CHECK_CLOSE_ABS(track.chi2(), 6.5, 2.);
//  BOOST_CHECK_EQUAL(track.nDoF(), 10u);
//  BOOST_CHECK_EQUAL(track.nHoles(), 1u);
//  BOOST_CHECK_EQUAL(track.nMeasurements(), nMeasurements);
//  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
//  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);
//
//  // Parameters
//  // We need quite coarse checks here, since on different builds
//  // the created measurements differ in the randomness
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -11., 7e0);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -15., 6e0);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1e-5, 1e3);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], M_PI / 2, 1e-3);
//  BOOST_CHECK_EQUAL(track.parameters()[eBoundQOverP], 1);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime],
//                    startParametersFit.parameters()[eBoundTime], 1e-6);
//  BOOST_CHECK_CLOSE(track.covariance().determinant(), 4.7e-28, 2e0);
//
//  ACTS_INFO("*** Test: FindHoles -- Finish");
//}
//
// BOOST_AUTO_TEST_CASE(Material) {
//  ACTS_INFO("*** Test: Material -- Start");
//
//  std::default_random_engine rng(42);
//
//  ACTS_DEBUG("Create the detector");
//  const std::size_t nSurfaces = 7;
//  const std::set<std::size_t> surfaceIndexWithMaterial = {4};
//  Detector detector;
//  detector.geometry =
//      makeToyDetector(geoCtx, nSurfaces, surfaceIndexWithMaterial);
//
//  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
//  const auto parametersMeasurements = makeParameters();
//  const auto startParametersFit = makeParameters(
//      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);
//
//  ACTS_DEBUG("Create the measurements");
//  using SimPropagator = Acts::Propagator<Acts::StraightLineStepper,
//                                         Acts::Experimental::DetectorNavigator>;
//  const SimPropagator simPropagator =
//  makeStraightPropagator(detector.geometry); auto measurements =
//      createMeasurements(simPropagator, geoCtx, magCtx,
//      parametersMeasurements,
//                         resMapAllPixel, rng);
//
//  const Acts::ActsVector<2> scatterOffset = {100_mm, 100_mm};
//  const std::size_t indexMaterialSurface = 3;
//  for (std::size_t iMeas = indexMaterialSurface; iMeas < nSurfaces; iMeas++) {
//    // This only works, because our detector is evenly spaced
//    const std::size_t offsetFactor = iMeas - indexMaterialSurface;
//
//    auto& sl = measurements.sourceLinks[iMeas];
//    sl.parameters[0] += scatterOffset[0] * offsetFactor;
//    sl.parameters[1] += scatterOffset[1] * offsetFactor;
//  }
//
//  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
//  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());
//
//  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);
//
//  ACTS_DEBUG("Set up the fitter");
//  const Surface* rSurface = &parametersMeasurements.referenceSurface();
//
//  using RecoStepper = EigenStepper<>;
//  const auto recoPropagator =
//      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);
//
//  using RecoPropagator = decltype(recoPropagator);
//  using Gx2Fitter =
//      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
//  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());
//
//  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
//  extensions.calibrator
//      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
//  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
//  extensions.surfaceAccessor
//      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);
//
//  const Experimental::Gx2FitterOptions gx2fOptions(
//      geoCtx, magCtx, calCtx, extensions,
//      PropagatorPlainOptions(geoCtx, magCtx), rSurface, false, false,
//      FreeToBoundCorrection(false), 5, 0);
//
//  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
//                              Acts::VectorMultiTrajectory{}};
//
//  ACTS_DEBUG("Fit the track");
//  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
//  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
//  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
//                              startParametersFit, gx2fOptions, tracks);
//
//  BOOST_REQUIRE(res.ok());
//
//  const auto& track = *res;
//
//  BOOST_CHECK_EQUAL(track.tipIndex(), nSurfaces - 1);
//  BOOST_CHECK(track.hasReferenceSurface());
//
//  // TODO Add material handling to the gx2f, to pass the 6 commented tests
//  // Track quantities
//  //  CHECK_CLOSE_ABS(track.chi2(), 8., 2.);
//  BOOST_CHECK_EQUAL(track.nDoF(), nSurfaces * 2);
//  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
//  BOOST_CHECK_EQUAL(track.nMeasurements(), nSurfaces);
//  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
//  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);
//
//  // Parameters
//  // We need quite coarse checks here, since on different builds
//  // the created measurements differ in the randomness
//  //  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -11., 7e0);
//  //  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -15., 6e0);
//  //  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1e-5, 1e3);
//  //  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], M_PI / 2, 1e-3);
//  BOOST_CHECK_EQUAL(track.parameters()[eBoundQOverP], 1);
//  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime],
//                    startParametersFit.parameters()[eBoundTime], 1e-6);
//  //  BOOST_CHECK_CLOSE(track.covariance().determinant(), 1e-27, 4e0);
//
//  // Convergence
//  BOOST_CHECK_EQUAL(
//      (track.template component<
//          std::uint32_t,
//          hashString(Experimental::Gx2fConstants::gx2fnUpdateColumn)>()),
//      5);
//
//  ACTS_INFO("*** Test: Material -- Finish");
//}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test
