// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/tree/test_unit.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/CompositePortalLink.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/GridIterator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/TrackingGeometryJsonConverter.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"
#include "ActsTests/CommonHelpers/TemporaryDirectory.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <functional>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <nlohmann/json_fwd.hpp>

namespace ActsTests {

using namespace Acts;
using namespace Experimental;
using namespace UnitLiterals;
using enum CylinderVolumeBounds::Face;
using enum AxisDirection;

using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

template <typename navigator_t = Navigator,
          typename stepper_t = ConstantFieldStepper>
struct StateCollector {
  struct this_result {
    std::vector<typename navigator_t::State> navigation;
    std::vector<typename stepper_t::State> stepping;
  };

  using result_type = this_result;

  template <typename propagator_state_t>
  Result<void> act(propagator_state_t& state, const stepper_t& /*stepper*/,
                   const navigator_t& /*navigator*/, result_type& result,
                   const Logger& /*logger*/) const {
    result.navigation.push_back(state.navigation);
    result.stepping.push_back(state.stepping);
    return {};
  }
};

auto logger = getDefaultLogger("UnitTests", Logging::VERBOSE);

BOOST_AUTO_TEST_SUITE(JsonSuite)

// BOOST_AUTO_TEST_CASE(TrackingGeometryJsonConverterRoundTrip) {
//   using namespace Acts;
//
//   GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
//
//   auto root = std::make_shared<TrackingVolume>(
//       Transform3::Identity(),
//       std::make_shared<CuboidVolumeBounds>(5., 5., 5.), "root");
//   root->assignGeometryId(GeometryIdentifier{}.withVolume(1u));
//   TrackingVolume* rootPtr = root.get();
//
//   Transform3 childTransform = Transform3::Identity();
//   childTransform.pretranslate(Vector3{1., 0., 0.});
//   auto child = std::make_unique<TrackingVolume>(
//       childTransform, std::make_shared<CylinderVolumeBounds>(0.5, 1.0, 2.0),
//       "child");
//   child->assignGeometryId(GeometryIdentifier{}.withVolume(2u));
//   TrackingVolume* childPtr = child.get();
//   root->addVolume(std::move(child));
//
//   auto trivialBounds = std::make_shared<const RectangleBounds>(1., 1.);
//   auto trivialSurface =
//       Surface::makeShared<PlaneSurface>(Transform3::Identity(),
//       trivialBounds);
//   trivialSurface->assignGeometryId(
//       GeometryIdentifier{}.withVolume(1u).withExtra(11u));
//   auto trivialLink =
//       std::make_unique<TrivialPortalLink>(trivialSurface, *childPtr);
//   root->addPortal(std::make_shared<Portal>(Direction::AlongNormal(),
//                                            std::move(trivialLink)));
//
//   auto compositeBounds = std::make_shared<const RectangleBounds>(1., 1.);
//   Transform3 compositeTransformA = Transform3::Identity();
//   compositeTransformA.pretranslate(Vector3{-1., 0., 0.});
//   Transform3 compositeTransformB = Transform3::Identity();
//   compositeTransformB.pretranslate(Vector3{1., 0., 0.});
//   auto compositeSurfaceA =
//       Surface::makeShared<PlaneSurface>(compositeTransformA,
//       compositeBounds);
//   auto compositeSurfaceB =
//       Surface::makeShared<PlaneSurface>(compositeTransformB,
//       compositeBounds);
//   compositeSurfaceA->assignGeometryId(
//       GeometryIdentifier{}.withVolume(1u).withExtra(12u));
//   compositeSurfaceB->assignGeometryId(
//       GeometryIdentifier{}.withVolume(1u).withExtra(13u));
//   auto compositeLinkA =
//       std::make_unique<TrivialPortalLink>(compositeSurfaceA, *childPtr);
//   auto compositeLinkB =
//       std::make_unique<TrivialPortalLink>(compositeSurfaceB, *rootPtr);
//   auto compositeLink = std::make_unique<CompositePortalLink>(
//       std::move(compositeLinkA), std::move(compositeLinkB),
//       AxisDirection::AxisX);
//   root->addPortal(std::make_shared<Portal>(Direction::AlongNormal(),
//                                            std::move(compositeLink)));
//
//   auto gridBounds = std::make_shared<const RectangleBounds>(2., 1.);
//   auto gridSurface =
//       Surface::makeShared<PlaneSurface>(Transform3::Identity(), gridBounds);
//   gridSurface->assignGeometryId(
//       GeometryIdentifier{}.withVolume(1u).withExtra(14u));
//   auto gridLink = GridPortalLink::make(gridSurface, AxisDirection::AxisX,
//                                        Axis{AxisBound, -2., 2., 2});
//   AnyGridView<const TrackingVolume*> gridView(gridLink->grid());
//   gridView.atLocalBins({0u}) = rootPtr;
//   gridView.atLocalBins({1u}) = childPtr;
//   gridView.atLocalBins({2u}) = rootPtr;
//   gridView.atLocalBins({3u}) = childPtr;
//
//   std::vector<TrivialPortalLink> artifacts;
//   auto artifactSurface =
//       Surface::makeShared<PlaneSurface>(Transform3::Identity(), gridBounds);
//   artifactSurface->assignGeometryId(
//       GeometryIdentifier{}.withVolume(1u).withExtra(15u));
//   artifacts.emplace_back(artifactSurface, *childPtr);
//   gridLink->setArtifactPortalLinks(std::move(artifacts));
//
//   root->addPortal(
//       std::make_shared<Portal>(Direction::AlongNormal(),
//       std::move(gridLink)));
//
//   auto sharedPortalBounds = std::make_shared<const RectangleBounds>(0.75,
//   0.75); auto sharedPortalSurface = Surface::makeShared<PlaneSurface>(
//       Transform3::Identity(), sharedPortalBounds);
//   sharedPortalSurface->assignGeometryId(
//       GeometryIdentifier{}.withVolume(1u).withExtra(16u));
//
//   auto sharedPortal = std::make_shared<Portal>(
//       gctx, std::make_unique<TrivialPortalLink>(sharedPortalSurface,
//       *childPtr), std::make_unique<TrivialPortalLink>(sharedPortalSurface,
//       *rootPtr));
//   root->addPortal(sharedPortal);
//   childPtr->addPortal(sharedPortal);
//
//   TrackingGeometryJsonConverter converter;
//   nlohmann::json encoded = converter.toJson(gctx, *root);
//   TemporaryDirectory tmpDir{};
//   auto jsonPath = tmpDir.path() / "tracking_geometry_roundtrip.json";
//
//   {
//     std::ofstream out(jsonPath);
//     BOOST_REQUIRE(out.good());
//     out << encoded.dump(2);
//   }
//
//   nlohmann::json encodedFromFile;
//   {
//     std::ifstream in(jsonPath);
//     BOOST_REQUIRE(in.good());
//     in >> encodedFromFile;
//   }
//
//   auto decodedRoot = converter.trackingVolumeFromJson(gctx, encodedFromFile);
//
//   BOOST_REQUIRE(decodedRoot != nullptr);
//   BOOST_CHECK_EQUAL(decodedRoot->volumeName(), "root");
//   BOOST_CHECK_EQUAL(decodedRoot->volumeBounds().type(),
//   VolumeBounds::eCuboid);
//
//   std::vector<TrackingVolume*> decodedChildren;
//   for (auto& decodedChild : decodedRoot->volumes()) {
//     decodedChildren.push_back(&decodedChild);
//   }
//   BOOST_REQUIRE_EQUAL(decodedChildren.size(), 1u);
//   BOOST_CHECK_EQUAL(decodedChildren.front()->volumeName(), "child");
//   BOOST_CHECK_EQUAL(decodedChildren.front()->volumeBounds().type(),
//                     VolumeBounds::eCylinder);
//
//   std::vector<Portal*> decodedPortals;
//   for (auto& portal : decodedRoot->portals()) {
//     decodedPortals.push_back(&portal);
//   }
//   BOOST_REQUIRE_EQUAL(decodedPortals.size(), 4u);
//
//   const auto* decodedTrivial = dynamic_cast<const TrivialPortalLink*>(
//       decodedPortals.at(0)->getLink(Direction::AlongNormal()));
//   BOOST_REQUIRE(decodedTrivial != nullptr);
//   BOOST_CHECK_EQUAL(decodedTrivial->volume().volumeName(), "child");
//
//   const auto* decodedComposite = dynamic_cast<const CompositePortalLink*>(
//       decodedPortals.at(1)->getLink(Direction::AlongNormal()));
//   BOOST_REQUIRE(decodedComposite != nullptr);
//   BOOST_CHECK_EQUAL(decodedComposite->size(), 2u);
//   BOOST_CHECK_EQUAL(decodedComposite->direction(), AxisDirection::AxisX);
//
//   const auto* decodedGrid = dynamic_cast<const GridPortalLink*>(
//       decodedPortals.at(2)->getLink(Direction::AlongNormal()));
//   BOOST_REQUIRE(decodedGrid != nullptr);
//   BOOST_CHECK_EQUAL(decodedGrid->dim(), 1u);
//   BOOST_CHECK_EQUAL(decodedGrid->artifactPortalLinks().size(), 1u);
//
//   AnyGridConstView<const TrackingVolume*>
//   decodedGridView(decodedGrid->grid());
//   BOOST_REQUIRE(decodedGridView.atLocalBins({0u}) != nullptr);
//   BOOST_REQUIRE(decodedGridView.atLocalBins({1u}) != nullptr);
//   BOOST_REQUIRE(decodedGridView.atLocalBins({2u}) != nullptr);
//   BOOST_REQUIRE(decodedGridView.atLocalBins({3u}) != nullptr);
//
//   BOOST_CHECK_EQUAL(decodedGridView.atLocalBins({0u})->volumeName(), "root");
//   BOOST_CHECK_EQUAL(decodedGridView.atLocalBins({1u})->volumeName(),
//   "child");
//   BOOST_CHECK_EQUAL(decodedGridView.atLocalBins({2u})->volumeName(), "root");
//   BOOST_CHECK_EQUAL(decodedGridView.atLocalBins({3u})->volumeName(),
//   "child");
//
//   std::vector<Portal*> decodedChildPortals;
//   for (auto& portal : decodedChildren.front()->portals()) {
//     decodedChildPortals.push_back(&portal);
//   }
//   BOOST_REQUIRE_EQUAL(decodedChildPortals.size(), 1u);
//
//   bool sharedPortalPreserved = false;
//   for (Portal* rootPortal : decodedPortals) {
//     if (rootPortal == decodedChildPortals.front()) {
//       sharedPortalPreserved = true;
//       break;
//     }
//   }
//   BOOST_CHECK(sharedPortalPreserved);
//
//   auto decodedGeometry =
//       converter.trackingGeometryFromJson(gctx, encodedFromFile);
//   BOOST_REQUIRE(decodedGeometry != nullptr);
//   BOOST_REQUIRE(decodedGeometry->highestTrackingVolume() != nullptr);
//   BOOST_CHECK_EQUAL(decodedGeometry->highestTrackingVolume()->volumeName(),
//                     "root");
// }

// BOOST_AUTO_TEST_CASE(TrackingGeometryJsonConverterRoundTripGen3Cylindrical) {
//   using namespace Acts;
//
//   GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
//
//   CylindricalTrackingGeometry cylindricalGeometryBuilder(gctx, true);
//   auto sourceGeometry = cylindricalGeometryBuilder();
//
//   BOOST_REQUIRE(sourceGeometry != nullptr);
//   BOOST_REQUIRE(sourceGeometry->highestTrackingVolume() != nullptr);
//   BOOST_CHECK(sourceGeometry->geometryVersion() ==
//               TrackingGeometry::GeometryVersion::Gen3);
//
//   auto countVolumesAndPortals = [](const TrackingVolume& world) {
//     std::size_t volumeCount = 0u;
//     std::size_t portalCount = 0u;
//
//     std::function<void(const TrackingVolume&)> traverse =
//         [&](const TrackingVolume& volume) {
//           ++volumeCount;
//           for (const auto& portal : volume.portals()) {
//             static_cast<void>(portal);
//             ++portalCount;
//           }
//           for (const auto& child : volume.volumes()) {
//             traverse(child);
//           }
//         };
//
//     traverse(world);
//     return std::pair{volumeCount, portalCount};
//   };
//
//   const auto [sourceVolumeCount, sourcePortalCount] =
//       countVolumesAndPortals(*sourceGeometry->highestTrackingVolume());
//   BOOST_CHECK_GT(sourceVolumeCount, 0u);
//   BOOST_CHECK_GT(sourcePortalCount, 0u);
//
//   TrackingGeometryJsonConverter converter;
//   nlohmann::json encoded = converter.toJson(gctx, *sourceGeometry);
//
//   TemporaryDirectory tmpDir{};
//   auto jsonPath = "tracking_geometry_gen3_roundtrip.json";
//
//   {
//     std::ofstream out(jsonPath);
//     BOOST_REQUIRE(out.good());
//     out << encoded.dump(2);
//   }
//
//   nlohmann::json encodedFromFile;
//   {
//     std::ifstream in(jsonPath);
//     BOOST_REQUIRE(in.good());
//     in >> encodedFromFile;
//   }
//
//   auto decodedGeometry =
//       converter.trackingGeometryFromJson(gctx, encodedFromFile);
//   BOOST_REQUIRE(decodedGeometry != nullptr);
//   BOOST_REQUIRE(decodedGeometry->highestTrackingVolume() != nullptr);
//
//   const auto [decodedVolumeCount, decodedPortalCount] =
//       countVolumesAndPortals(*decodedGeometry->highestTrackingVolume());
//
//   BOOST_CHECK_EQUAL(decodedVolumeCount, sourceVolumeCount);
//   BOOST_CHECK_EQUAL(decodedPortalCount, sourcePortalCount);
//   BOOST_CHECK_EQUAL(decodedGeometry->highestTrackingVolume()->volumeName(),
//                     sourceGeometry->highestTrackingVolume()->volumeName());
// }

void checkHierarchy(const GeometryContext& gctx,
                    const TrackingVolume::VolumeRange& volsA,
                    const TrackingVolume::VolumeRange& volsB) {
  auto checkSurfaces = [&](const auto& surfA, const auto& surfB) {
    std::cout << surfA.type() << " -- " << surfB.type() << "\n";
    BOOST_CHECK_EQUAL(surfA.type(), surfB.type());
    std::cout << surfA.geometryId() << " -- " << surfB.geometryId() << "\n";
    BOOST_CHECK_EQUAL(surfA.geometryId(), surfB.geometryId());
    std::cout << surfA.bounds() << " -- " << surfB.bounds() << "\n";
    BOOST_CHECK_EQUAL(surfA.bounds(), surfB.bounds());
    std::cout << (surfA.surfacePlacement() != nullptr) << " -- "
              << (surfB.surfacePlacement() != nullptr) << "\n";

    std::cout << "\n"
              << surfA.localToGlobalTransform(gctx).matrix() << " --\n"
              << surfB.localToGlobalTransform(gctx).matrix() << "\n";
    BOOST_CHECK_LT((surfA.localToGlobalTransform(gctx).matrix() -
                    surfB.localToGlobalTransform(gctx).matrix())
                       .norm(),
                   1e-10);
    const auto* matA = surfA.surfaceMaterial();
    const auto* matB = surfB.surfaceMaterial();
    if (matA != nullptr) {
      BOOST_CHECK_NE(matB, nullptr);
      BOOST_CHECK_EQUAL(surfA.surfaceMaterial()->materialSlab(Vector2{0, 0}),
                        surfB.surfaceMaterial()->materialSlab(Vector2{0, 0}));
    }
  };
  std::function<void(const Acts::PortalLinkBase* linkA,
                     const Acts::PortalLinkBase* linkB)>
      checkLinks = [&](const auto* linkA, const auto* linkB) {
        const auto* trivialA =
            dynamic_cast<const Acts::TrivialPortalLink*>(linkA);
        const auto* trivialB =
            dynamic_cast<const Acts::TrivialPortalLink*>(linkB);
        if (trivialA != nullptr) {
          BOOST_CHECK_NE(trivialB, nullptr);
          BOOST_CHECK(trivialA->volume() == trivialB->volume());
          checkSurfaces(trivialA->surface(), trivialB->surface());
        }

        const auto* gridA = dynamic_cast<const Acts::GridPortalLink*>(linkA);
        const auto* gridB = dynamic_cast<const Acts::GridPortalLink*>(linkB);
        if (gridA != nullptr) {
          BOOST_CHECK_NE(gridB, nullptr);
          BOOST_CHECK_EQUAL(gridA->dim(), gridB->dim());
          BOOST_CHECK_EQUAL(gridA->direction(), gridB->direction());
          BOOST_CHECK(gridA->grid() == gridB->grid());
          checkSurfaces(gridA->surface(), gridB->surface());

          const auto& childrenA = gridA->artifactPortalLinks();
          const auto& childrenB = gridA->artifactPortalLinks();
          for (std::size_t i = 0; i < childrenA.size(); i++) {
            const auto& childA = *(childrenA.begin() + i);
            const auto& childB = *(childrenB.begin() + i);
            checkLinks(&childA, &childB);
          }
        }

        const auto* compositeA =
            dynamic_cast<const Acts::CompositePortalLink*>(linkA);
        const auto* compositeB =
            dynamic_cast<const Acts::CompositePortalLink*>(linkB);
        if (compositeA != nullptr) {
          BOOST_CHECK_NE(compositeB, nullptr);
          const auto& childrenA = compositeA->links();
          const auto& childrenB = compositeB->links();
          for (std::size_t i = 0; i < childrenA.size(); i++) {
            const auto& childA = childrenA.at(i);
            const auto& childB = childrenB.at(i);
            checkLinks(&childA, &childB);
          }
        }
      };

  BOOST_CHECK_EQUAL(volsA.size(), volsB.size());
  for (std::size_t i = 0; i < volsA.size(); i++) {
    std::cout << "------------------------------\n";
    const auto& volA = volsA.at(i);
    const auto& volB = volsB.at(i);
    BOOST_CHECK(volA == volB);
    std::cout << volA.volumeName() << " -- " << volB.volumeName() << "\n";

    const auto& surfsA = volA.surfaces();
    const auto& surfsB = volB.surfaces();
    BOOST_CHECK_EQUAL(surfsA.size(), surfsB.size());
    for (std::size_t j = 0; j < surfsA.size(); j++) {
      std::cout << surfsA.at(i).geometryId() << " -- "
                << surfsB.at(i).geometryId() << "\n";
      checkSurfaces(surfsA.at(i), surfsB.at(i));
    }

    std::cout << "PORTALS CHECK\n";
    const auto& portsA = volA.portals();
    const auto& portsB = volB.portals();
    BOOST_CHECK_EQUAL(portsA.size(), portsB.size());
    for (std::size_t j = 0; j < portsA.size(); j++) {
      const auto& portA = portsA.at(j);
      const auto& portB = portsB.at(j);

      const auto& surfA = portA.surface();
      const auto& surfB = portB.surface();
      checkSurfaces(surfA, surfB);

      const auto* alongA = portA.getLink(Acts::Direction::AlongNormal());
      const auto* alongB = portB.getLink(Acts::Direction::AlongNormal());
      checkLinks(alongA, alongB);

      const auto* oppositeA = portA.getLink(Acts::Direction::OppositeNormal());
      const auto* oppositeB = portB.getLink(Acts::Direction::OppositeNormal());
      checkLinks(oppositeA, oppositeB);
    }
    checkHierarchy(gctx, volA.volumes(), volB.volumes());
  }
}

class DummySurfacePlacement : public SurfacePlacementBase {
 public:
  explicit DummySurfacePlacement(const std::shared_ptr<Surface>& surface,
                                 std::size_t id, const GeometryContext& gctx)
      : m_transform(surface->localToGlobalTransform(gctx)),
        m_surface(surface),
        m_id(id) {}
  ~DummySurfacePlacement() override = default;

  const Transform3& localToGlobalTransform(
      const GeometryContext& /*gctx*/) const override {
    std::cout << "DUMMY SUFACE PLACEMENT TRANSFORM CALL\n";
    std::cout << m_transform.matrix() << "\n";
    return m_transform;
  };

  const Surface& surface() const override {
    std::cout << "DUMMY SURFACE PLACEMENT SURFACE CALL CONST "
              << (m_surface == nullptr) << "\n";
    return *m_surface;
  }

  Surface& surface() override {
    std::cout << "DUMMY SURFACE PLACEMENT SURFACE CALL "
              << (m_surface == nullptr) << "\n";

    return *m_surface;
  }

  bool isSensitive() const override { return true; }

  std::size_t id() const { return m_id; }

 private:
  Transform3 m_transform;
  std::shared_ptr<Surface> m_surface;
  std::size_t m_id;
};

nlohmann::json encodeDummySurfacePlacement(
    const DummySurfacePlacement& placement, const GeometryContext& /*gctx*/) {
  nlohmann::json jPlacement;
  jPlacement["kind"] = "DummySurfacePlacement";
  jPlacement["id"] = placement.id();
  return jPlacement;
}

std::shared_ptr<SurfacePlacementBase> decodeDummySurfacePlacement(
    const nlohmann::json& encoded, const GeometryContext& gctx,
    const std::shared_ptr<Surface>& surface) {
  auto id = encoded["id"].get<std::size_t>();
  return std::make_shared<DummySurfacePlacement>(surface, id, gctx);
}

BOOST_AUTO_TEST_CASE(test) {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  MagneticFieldContext mctx;

  auto field = std::make_shared<ConstantBField>(Vector3{0, 0, 0_T});
  EigenStepper<> stepper{field};

  BoundTrackParameters start = BoundTrackParameters::createCurvilinear(
      Vector4::Zero(), Vector3::UnitX(), 1. / 1_GeV, std::nullopt,
      ParticleHypothesis::pion());

  // Action list and abort list
  using EndOfWorld = EndOfWorldReached;
  using ReferenceActorList = ActorList<StateCollector<>, EndOfWorld>;
  using PropagatorOptions =
      typename ConstantFieldPropagator::template Options<ReferenceActorList>;

  // Options definition
  PropagatorOptions options(gctx, mctx);

  // -----------------------------------------------------------------------------
  // Build geomtery
  TryAllNavigationPolicy::Config tryAllConfig;
  tryAllConfig.portals = true;
  tryAllConfig.sensitives = false;

  // Create the navigation policy factory
  std::unique_ptr<NavigationPolicyFactory> navPolicyFactory =
      NavigationPolicyFactory{}
          .add<TryAllNavigationPolicy>(tryAllConfig)
          .asUniquePtr();

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {2_mm, 20_mm};
  auto root = std::make_unique<Blueprint>(cfg);

  auto& cyl = root->addCylinderContainer("Container", AxisDirection::AxisZ);

  double rStep = 10_mm;
  double hlZ = 30_mm;
  std::vector<std::shared_ptr<DummySurfacePlacement>> placements;
  for (std::size_t i = 0; i < 1; i++) {
    auto cylBounds =
        std::make_shared<CylinderVolumeBounds>(rStep * i, rStep * (i + 1), hlZ);

    auto childCylL = std::make_unique<TrackingVolume>(
        Transform3::Identity() * Translation3{Vector3{0, 0, -hlZ}}, cylBounds,
        "childL" + std::to_string(i));
    childCylL->setNavigationPolicy(
        navPolicyFactory->build(gctx, *childCylL, *logger));

    auto childCylR = std::make_unique<TrackingVolume>(
        Transform3::Identity() * Translation3{Vector3{0, 0, hlZ}}, cylBounds,
        "childR" + std::to_string(i));
    childCylR->setNavigationPolicy(
        navPolicyFactory->build(gctx, *childCylR, *logger));

    auto planeBounds = std::make_shared<Acts::RectangleBounds>(5, 5);
    auto planeL = Acts::Surface::makeShared<PlaneSurface>(
        Transform3::Identity(), planeBounds);
    placements.push_back(std::make_shared<DummySurfacePlacement>(
        planeL, placements.size(), gctx));
    planeL->assignSurfacePlacement(*placements.back());

    childCylL->addSurface(planeL);

    auto childNodeL =
        std::make_shared<StaticBlueprintNode>(std::move(childCylL));
    auto childNodeR =
        std::make_shared<StaticBlueprintNode>(std::move(childCylR));

    cyl.addChild(childNodeL);
    cyl.addChild(childNodeR);
  }

  auto sourceGeometry = root->construct({}, gctx, *logger);

  std::cout << "-----------------------------------\n";
  std::cout << "SOURCE GEO PLACEMENTS:\n";
  sourceGeometry->visitSurfaces([&](const auto* surf) {
    std::cout << "SOURCE GEO ID " << surf->geometryId() << "\n";
    std::cout << "PLACEMENT " << (surf->surfacePlacement() == nullptr) << "\n";
    if (surf->surfacePlacement() != nullptr) {
      std::cout << "SOURCE DUMMY PLACEMENT "
                << dynamic_cast<const DummySurfacePlacement*>(
                       surf->surfacePlacement())
                       ->id()
                << "\n";
      std::cout << "SOURCE DUMMY TRANSFORM\n"
                << surf->localToGlobalTransform(gctx).matrix() << "\n";
    }
  });
  std::cout << "-----------------------------------\n";

  // -----------------------------------------------------
  // Encode-decode geometry
  // auto converterCfg = TrackingGeometryJsonConverter::Config::defaultConfig();
  // converterCfg.surfacePlacementEncoder.registerFunction(
  //     encodeDummySurfacePlacement);
  // converterCfg.surfacePlacementDecoder.registerKind(
  //     "DummySurfacePlacement", decodeDummySurfacePlacement);
  // -----------------------------------------------------
  TrackingGeometryJsonConverter converter;
  nlohmann::json encoded = converter.toJson(gctx, *sourceGeometry);

  auto jsonPath = "tracking_geometry_gen3_roundtrip.json";

  {
    std::ofstream out(jsonPath);
    BOOST_REQUIRE(out.good());
    out << encoded.dump(2);
  }

  nlohmann::json encodedFromFile;
  {
    std::ifstream in(jsonPath);
    BOOST_REQUIRE(in.good());
    in >> encodedFromFile;
  }

  auto decodedGeometry =
      converter.trackingGeometryFromJson(gctx, encodedFromFile);

  std::cout << "-----------------------------------\n";
  std::cout << "DECODED GEO PLACEMENTS:\n";
  decodedGeometry->visitSurfaces([&](const auto* surf) {
    std::cout << "DECODED GEO ID " << surf->geometryId() << "\n";
    std::cout << "DECODED PLACEMENT " << (surf->surfacePlacement() == nullptr)
              << "\n";
    if (surf->surfacePlacement() != nullptr) {
      std::cout << "DECODED DUMMY PLACEMENT "
                << dynamic_cast<const DummySurfacePlacement*>(
                       surf->surfacePlacement())
                       ->id()
                << "\n";
      std::cout << "DECODED DUMMY TRANSFORM\n"
                << surf->localToGlobalTransform(gctx).matrix() << "\n";
    }
  });
  std::cout << "-----------------------------------\n";

  const auto* htvSource = sourceGeometry->highestTrackingVolume();
  const auto* htvDecoded = decodedGeometry->highestTrackingVolume();
  BOOST_CHECK(*htvSource == *htvDecoded);

  std::cout << "HIERARCHY CHECK START\n";
  checkHierarchy(gctx, htvSource->volumes(), htvDecoded->volumes());
  std::cout << "HIERARCHY CHECK END\n";

  // // -----------------------------------------------------
  // Navigator::Config sourceNavCfg;
  // sourceNavCfg.trackingGeometry =
  //     std::make_shared<TrackingGeometry>(*sourceGeometry);
  // sourceNavCfg.resolveSensitive = true;
  // sourceNavCfg.resolveMaterial = true;
  // sourceNavCfg.resolvePassive = false;
  // Navigator sourceNavigator{
  //     sourceNavCfg,
  //     Acts::getDefaultLogger("SourceNavigator", Acts::Logging::VERBOSE)};
  // ConstantFieldPropagator sourcePropagator(stepper, sourceNavigator);

  // auto sourceRes = sourcePropagator.propagate(start, options).value();

  // // -----------------------------------------------------
  // Navigator::Config decodedNavCfg;
  // decodedNavCfg.trackingGeometry =
  //     std::make_shared<TrackingGeometry>(*decodedGeometry);
  // decodedNavCfg.resolveSensitive = true;
  // decodedNavCfg.resolveMaterial = true;
  // decodedNavCfg.resolvePassive = false;
  // Navigator decodedNavigator{
  //     sourceNavCfg,
  //     Acts::getDefaultLogger("DecodedNavigator", Acts::Logging::VERBOSE)};
  // ConstantFieldPropagator decodedPropagator(stepper, decodedNavigator);

  // auto decodedRes = decodedPropagator.propagate(start, options).value();

  // // -----------------------------------------------------
  // // Compare propagation
  // auto& sourceStates = sourceRes.template
  // get<StateCollector<>::result_type>(); auto& decodedStates =
  //     decodedRes.template get<StateCollector<>::result_type>();
  // for (std::size_t i = 0; i < sourceStates.navigation.size(); i++) {
  //   auto& navigationA = sourceStates.navigation.at(i);
  //   auto& navigationB = decodedStates.navigation.at(i);
  //   if (!navigationA.navigationBreak) {
  //     BOOST_CHECK(*navigationA.currentVolume == *navigationB.currentVolume);
  //     BOOST_CHECK(*navigationA.startVolume == *navigationB.startVolume);
  //     BOOST_CHECK(*navigationA.currentSurface ==
  //     *navigationB.currentSurface); BOOST_CHECK(*navigationA.startSurface ==
  //     *navigationB.startSurface);
  //   } else {
  //     BOOST_CHECK_EQUAL(navigationA.currentVolume, nullptr);
  //     BOOST_CHECK_EQUAL(navigationB.currentVolume, nullptr);
  //   }

  //   auto& steppingA = sourceStates.stepping.at(i);
  //   auto& steppingB = decodedStates.stepping.at(i);
  //   BOOST_CHECK(steppingA.pars == steppingB.pars);
  //   BOOST_CHECK(steppingA.nSteps == steppingB.nSteps);
  //   BOOST_CHECK(steppingA.pathAccumulated == steppingB.pathAccumulated);
  // }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
