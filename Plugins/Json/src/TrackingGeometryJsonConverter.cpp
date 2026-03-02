// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/TrackingGeometryJsonConverter.hpp"

#include "Acts/Geometry/CompositePortalLink.hpp"
#include "Acts/Geometry/ConeVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GenericCuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/SurfacePlacementBase.hpp"
#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/AlgebraJsonConverter.hpp"
#include "ActsPlugins/Json/GeometryIdentifierJsonConverter.hpp"
#include "ActsPlugins/Json/SurfaceJsonConverter.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <ranges>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <nlohmann/json_fwd.hpp>

namespace {

/// Internal notes on the conversion flow used in this translation unit:
/// - Encoder side:
///   - normalize enum-like values (axis/shape metadata) to stable JSON strings
///   - encode bounds and portal-link polymorphic payloads using kind tags
///   - emit volumes in deterministic depth-first order and connect them by IDs
///   - emit unique portals in a top-level table and refer to them from volumes
/// - Decoder side:
///   - read and validate schema metadata
///   - construct all volumes first, then attach children and portals
///   - rebuild polymorphic portal links from kind-tagged payloads

constexpr const char* kHeaderKey = "acts-geometry-volumes";
constexpr const char* kVersionKey = "format-version";
constexpr const char* kScopeKey = "scope";
constexpr const char* kScopeValue = "volumes-bounds-portals";
constexpr int kFormatVersion = 1;

constexpr const char* kPortalsKey = "portals";
constexpr const char* kSurfacesKey = "surfaces";
constexpr const char* kVolumesKey = "volumes";

constexpr const char* kRootVolumeIdKey = "root_volume_id";
constexpr const char* kVolumeIdKey = "volume_id";
constexpr const char* kPortalIdKey = "portal_id";
constexpr const char* kSurfaceIdKey = "surface_id";

constexpr const char* kNameKey = "name";
constexpr const char* kGeometryIdKey = "geometry_id";
constexpr const char* kTransformKey = "transform";
constexpr const char* kBoundsKey = "bounds";
constexpr const char* kChildrenKey = "children";

constexpr const char* kPortalIdsKey = "portal_ids";
constexpr const char* kAlongNormalKey = "along_normal";
constexpr const char* kOppositeNormalKey = "opposite_normal";

constexpr const char* kNavigationPolicyKey = "navigation_policy";

constexpr const char* kKindKey = "kind";
constexpr const char* kValuesKey = "values";
constexpr const char* kTargetVolumeIdKey = "target_volume_id";
constexpr const char* kDirectionKey = "direction";
constexpr const char* kAxesKey = "axes";
constexpr const char* kBinsKey = "bins";
constexpr const char* kLocalKey = "local";
constexpr const char* kArtifactLinksKey = "artifact_links";
constexpr const char* kBoundaryTypeKey = "boundary_type";
constexpr const char* kAxisTypeKey = "axis_type";
constexpr const char* kMinKey = "min";
constexpr const char* kMaxKey = "max";
constexpr const char* kNBinsKey = "n_bins";
constexpr const char* kEdgesKey = "edges";

constexpr const char* kTrivialPortalLinkKind = "TrivialPortalLink";
constexpr const char* kCompositePortalLinkKind = "CompositePortalLink";
constexpr const char* kGridPortalLinkKind = "GridPortalLink";

constexpr const char* kConeVolumeBoundsKind = "ConeVolumeBounds";
constexpr const char* kCuboidVolumeBoundsKind = "CuboidVolumeBounds";
constexpr const char* kCutoutCylinderVolumeBoundsKind =
    "CutoutCylinderVolumeBounds";
constexpr const char* kCylinderVolumeBoundsKind = "CylinderVolumeBounds";
constexpr const char* kGenericCuboidVolumeBoundsKind =
    "GenericCuboidVolumeBounds";
constexpr const char* kTrapezoidVolumeBoundsKind = "TrapezoidVolumeBounds";

constexpr const char* kTryAllNavigationPolicyKind = "TryAllNavigationPolicy";
constexpr const char* kSurfaceArrayNavigationPolicyKind =
    "SurfaceArrayNavigationPolicy";
constexpr const char* kMultiNavigationPolicyKind = "MultiNavigationPolicy";

/// Convert axis boundary enum to stable JSON string representation.
std::string axisBoundaryTypeToString(Acts::AxisBoundaryType boundaryType) {
  using enum Acts::AxisBoundaryType;
  switch (boundaryType) {
    case Open:
      return "Open";
    case Bound:
      return "Bound";
    case Closed:
      return "Closed";
    default:
      throw std::invalid_argument("Unknown AxisBoundaryType");
  }
}

/// Convert axis boundary string from JSON back to enum representation.
Acts::AxisBoundaryType axisBoundaryTypeFromString(
    const std::string& encodedBoundaryType) {
  if (encodedBoundaryType == "Open") {
    return Acts::AxisBoundaryType::Open;
  }
  if (encodedBoundaryType == "Bound") {
    return Acts::AxisBoundaryType::Bound;
  }
  if (encodedBoundaryType == "Closed") {
    return Acts::AxisBoundaryType::Closed;
  }
  throw std::invalid_argument("Unknown AxisBoundaryType string: " +
                              encodedBoundaryType);
}

/// Convert axis type enum to stable JSON string representation.
std::string axisTypeToString(Acts::AxisType axisType) {
  using enum Acts::AxisType;
  switch (axisType) {
    case Equidistant:
      return "Equidistant";
    case Variable:
      return "Variable";
    default:
      throw std::invalid_argument("Unknown AxisType");
  }
}

/// Convert axis type string from JSON back to enum representation.
Acts::AxisType axisTypeFromString(const std::string& encodedAxisType) {
  if (encodedAxisType == "Equidistant") {
    return Acts::AxisType::Equidistant;
  }
  if (encodedAxisType == "Variable") {
    return Acts::AxisType::Variable;
  }
  throw std::invalid_argument("Unknown AxisType string: " + encodedAxisType);
}

/// Serialize one generic axis payload used by grid portal links.
nlohmann::json axisToJson(const Acts::IAxis& axis) {
  nlohmann::json jAxis;
  jAxis[kAxisTypeKey] = axisTypeToString(axis.getType());
  jAxis[kBoundaryTypeKey] = axisBoundaryTypeToString(axis.getBoundaryType());
  jAxis[kMinKey] = axis.getMin();
  jAxis[kMaxKey] = axis.getMax();
  jAxis[kNBinsKey] = axis.getNBins();
  jAxis[kEdgesKey] = axis.getBinEdges();
  return jAxis;
}

/// Deserialize one generic axis payload used by grid portal links.
std::unique_ptr<Acts::IAxis> axisFromJson(const nlohmann::json& jAxis) {
  Acts::AxisType axisType = axisTypeFromString(jAxis.at(kAxisTypeKey));
  Acts::AxisBoundaryType boundaryType =
      axisBoundaryTypeFromString(jAxis.at(kBoundaryTypeKey));

  if (axisType == Acts::AxisType::Equidistant) {
    return Acts::IAxis::createEquidistant(boundaryType, jAxis.at(kMinKey),
                                          jAxis.at(kMaxKey),
                                          jAxis.at(kNBinsKey));
  }

  return Acts::IAxis::createVariable(
      boundaryType, jAxis.at(kEdgesKey).get<std::vector<double>>());
}

/// Generic helper to decode concrete volume bounds from JSON "values".
template <typename bounds_t>
std::unique_ptr<Acts::VolumeBounds> decodeVolumeBoundsT(
    const nlohmann::json& jBounds) {
  constexpr std::size_t kValues = bounds_t::BoundValues::eSize;
  const auto values = jBounds.at(kValuesKey).get<std::vector<double>>();
  if (values.size() != kValues) {
    throw std::invalid_argument("Invalid number of values for volume bounds");
  }
  std::array<double, kValues> boundValues{};
  std::copy_n(values.begin(), kValues, boundValues.begin());
  return std::make_unique<bounds_t>(boundValues);
}

/// Generic helper to encode concrete volume bounds to kind + values payload.
template <typename bounds_t>
nlohmann::json encodeVolumeBoundsT(const bounds_t& bounds,
                                   std::string_view kind) {
  nlohmann::json jBounds;
  jBounds[kKindKey] = kind;
  jBounds[kValuesKey] = bounds.values();
  return jBounds;
}

/// Encode `ConeVolumeBounds` payload.
nlohmann::json encodeConeVolumeBounds(const Acts::ConeVolumeBounds& bounds) {
  return encodeVolumeBoundsT(bounds, kConeVolumeBoundsKind);
}

/// Encode `CuboidVolumeBounds` payload.
nlohmann::json encodeCuboidVolumeBounds(
    const Acts::CuboidVolumeBounds& bounds) {
  return encodeVolumeBoundsT(bounds, kCuboidVolumeBoundsKind);
}

/// Encode `CutoutCylinderVolumeBounds` payload.
nlohmann::json encodeCutoutCylinderVolumeBounds(
    const Acts::CutoutCylinderVolumeBounds& bounds) {
  return encodeVolumeBoundsT(bounds, kCutoutCylinderVolumeBoundsKind);
}

/// Encode `CylinderVolumeBounds` payload.
nlohmann::json encodeCylinderVolumeBounds(
    const Acts::CylinderVolumeBounds& bounds) {
  return encodeVolumeBoundsT(bounds, kCylinderVolumeBoundsKind);
}

/// Encode `GenericCuboidVolumeBounds` payload.
nlohmann::json encodeGenericCuboidVolumeBounds(
    const Acts::GenericCuboidVolumeBounds& bounds) {
  return encodeVolumeBoundsT(bounds, kGenericCuboidVolumeBoundsKind);
}

/// Encode `TrapezoidVolumeBounds` payload.
nlohmann::json encodeTrapezoidVolumeBounds(
    const Acts::TrapezoidVolumeBounds& bounds) {
  return encodeVolumeBoundsT(bounds, kTrapezoidVolumeBoundsKind);
}

std::unique_ptr<Acts::INavigationPolicy> decodeMultiNavigationPolicy(
    const nlohmann::json& encoded, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& converter,
    const Acts::TrackingVolume& volume, const Acts::Logger& logger) {
  std::vector<std::unique_ptr<Acts::INavigationPolicy>> children;
  for (const auto& child : encoded.at(kChildrenKey)) {
    children.push_back(converter.navigationPolicyFromJson(
        gctx, child[kNavigationPolicyKey], volume, logger));
  }
  return std::make_unique<Acts::MultiNavigationPolicy>(std::move(children));
}

std::unique_ptr<Acts::INavigationPolicy> decodeTryAllNavigationPolicy(
    const nlohmann::json& encoded, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& /*converter*/,
    const Acts::TrackingVolume& volume, const Acts::Logger& logger) {
  Acts::TryAllNavigationPolicy::Config cfg;
  cfg.passives = encoded.at("passives").get<bool>();
  cfg.sensitives = encoded.at("sensitives").get<bool>();
  cfg.portals = encoded.at("portals").get<bool>();

  return std::make_unique<Acts::TryAllNavigationPolicy>(gctx, volume, logger,
                                                        cfg);
}

std::unique_ptr<Acts::INavigationPolicy> decodeSurfaceArrayNavigationPolicy(
    const nlohmann::json& encoded, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& /*converter*/,
    const Acts::TrackingVolume& volume, const Acts::Logger& logger) {
  Acts::SurfaceArrayNavigationPolicy::Config cfg;
  cfg.layerType = encoded.at("layerType")
                      .get<Acts::SurfaceArrayNavigationPolicy::LayerType>();
  cfg.bins = {encoded.at("bins0").get<std::size_t>(),
              encoded.at("bins1").get<std::size_t>()};

  return std::make_unique<Acts::SurfaceArrayNavigationPolicy>(gctx, volume,
                                                              logger, cfg);
}

nlohmann::json encodeMultiNavigationPolicy(
    const Acts::MultiNavigationPolicy& policy,
    const Acts::TrackingGeometryJsonConverter& converter) {
  nlohmann::json jPolicy;
  jPolicy[kKindKey] = "MultiNavigationPolicy";
  for (const auto& pol : policy.policies()) {
    nlohmann::json jPol;
    jPol["navigation_policy"] = converter.navigationPolicyToJson(*pol);
    jPolicy[kChildrenKey].push_back(jPol);
  }
  return jPolicy;
}

nlohmann::json encodeTryAllNavigationPolicy(
    const Acts::TryAllNavigationPolicy& policy,
    const Acts::TrackingGeometryJsonConverter& /*converter*/) {
  const auto& cfg = policy.config();

  nlohmann::json jPolicy;
  jPolicy[kKindKey] = "TryAllNavigationPolicy";
  jPolicy["portals"] = cfg.portals;
  jPolicy["sensitives"] = cfg.sensitives;
  jPolicy["passives"] = cfg.passives;
  return jPolicy;
}

nlohmann::json encodeSurfaceArrayNavigationPolicy(
    const Acts::SurfaceArrayNavigationPolicy& policy,
    const Acts::TrackingGeometryJsonConverter& /*converter*/) {
  const auto& cfg = policy.config();

  nlohmann::json jPolicy;
  jPolicy[kKindKey] = "SurfaceArrayNavigationPolicy";
  jPolicy["layerType"] = cfg.layerType;
  jPolicy["bins0"] = cfg.bins.first;
  jPolicy["bins1"] = cfg.bins.second;
  return jPolicy;
}

/// Deserialize a portal surface and enforce `RegularSurface` type.
std::shared_ptr<Acts::RegularSurface> regularSurfaceFromJson(
    const nlohmann::json& jSurface) {
  auto surface = Acts::SurfaceJsonConverter::fromJson(jSurface);
  std::cout << "INSIDE TRACKING GEO CHECK "
            << (surface->surfacePlacement() == nullptr) << "\n";
  auto regular = std::dynamic_pointer_cast<Acts::RegularSurface>(surface);
  if (regular == nullptr) {
    throw std::invalid_argument("Portal link surface is not a RegularSurface");
  }
  return regular;
}

/// Build a concrete `GridPortalLink` from decoded surface and axis
/// payload(s).
std::unique_ptr<Acts::GridPortalLink> makeGridPortalLink(
    const std::shared_ptr<Acts::RegularSurface>& surface,
    Acts::AxisDirection direction, const Acts::IAxis& axis0,
    const Acts::IAxis* axis1) {
  std::unique_ptr<Acts::GridPortalLink> grid;

  if (axis1 == nullptr) {
    axis0.visit([&](const auto& a0) {
      using axis_t = std::remove_cvref_t<decltype(a0)>;
      axis_t axisCopy = a0;
      grid =
          Acts::GridPortalLink::make(surface, direction, std::move(axisCopy));
    });
  } else {
    axis0.visit([&](const auto& a0) {
      using axis0_t = std::remove_cvref_t<decltype(a0)>;
      axis0_t axis0Copy = a0;
      axis1->visit([&](const auto& a1) {
        using axis1_t = std::remove_cvref_t<decltype(a1)>;
        axis1_t axis1Copy = a1;
        grid = Acts::GridPortalLink::make(surface, std::move(axis0Copy),
                                          std::move(axis1Copy));
      });
    });
  }

  if (grid == nullptr) {
    throw std::invalid_argument("Could not construct GridPortalLink from axes");
  }

  if (grid->direction() != direction) {
    throw std::invalid_argument(
        "Decoded grid direction does not match payload");
  }

  return grid;
}

/// Encode a `TrivialPortalLink` payload including target volume ID.
nlohmann::json encodeTrivialPortalLink(
    const Acts::TrivialPortalLink& link, const Acts::GeometryContext& /*gctx*/,
    const Acts::TrackingGeometryJsonConverter& /*converter*/,
    const Acts::TrackingGeometryJsonConverter::SurfaceIdLookup& surfaceIds,
    const Acts::TrackingGeometryJsonConverter::VolumeIdLookup& volumeIds) {
  nlohmann::json jLink;
  jLink[kKindKey] = kTrivialPortalLinkKind;
  jLink[kSurfaceIdKey] = surfaceIds.at(link.surface());
  jLink[kTargetVolumeIdKey] = volumeIds.at(link.volume());
  return jLink;
}

/// Encode a `CompositePortalLink` payload including recursively encoded
/// children.
nlohmann::json encodeCompositePortalLink(
    const Acts::CompositePortalLink& link, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& converter,
    const Acts::TrackingGeometryJsonConverter::SurfaceIdLookup& surfaceIds,
    const Acts::TrackingGeometryJsonConverter::VolumeIdLookup& volumeIds) {
  nlohmann::json jLink;
  jLink[kKindKey] = kCompositePortalLinkKind;
  jLink[kSurfaceIdKey] = surfaceIds.at(link.surface());
  jLink[kDirectionKey] = link.direction();
  jLink[kChildrenKey] = nlohmann::json::array();

  for (const auto& child : link.links()) {
    jLink[kChildrenKey].push_back(
        converter.portalLinkToJson(gctx, child, surfaceIds, volumeIds));
  }
  return jLink;
}

/// Encode a `GridPortalLink` payload with axes, bins and artifact links.
nlohmann::json encodeGridPortalLink(
    const Acts::GridPortalLink& link, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& converter,
    const Acts::TrackingGeometryJsonConverter::SurfaceIdLookup& surfaceIds,
    const Acts::TrackingGeometryJsonConverter::VolumeIdLookup& volumeIds) {
  nlohmann::json jLink;
  jLink[kKindKey] = kGridPortalLinkKind;
  jLink[kDirectionKey] = link.direction();
  jLink[kSurfaceIdKey] = surfaceIds.at(link.surface());
  jLink[kAxesKey] = nlohmann::json::array();

  for (const auto* axis : link.grid().axes()) {
    jLink[kAxesKey].push_back(axisToJson(*axis));
  }

  Acts::AnyGridConstView<const Acts::TrackingVolume*> view(link.grid());
  const auto nBins = view.numLocalBins();
  const auto dim = view.dimensions();

  jLink[kBinsKey] = nlohmann::json::array();
  if (dim == 1u) {
    for (std::size_t i0 = 0u; i0 <= nBins.at(0) + 1u; ++i0) {
      nlohmann::json jBin;
      jBin[kLocalKey] = std::vector<std::size_t>{i0};
      const auto* target = view.atLocalBins({i0});
      if (target == nullptr) {
        jBin[kTargetVolumeIdKey] = nullptr;
      } else {
        jBin[kTargetVolumeIdKey] = volumeIds.at(*target);
      }
      jLink[kBinsKey].push_back(std::move(jBin));
    }
  } else if (dim == 2u) {
    for (std::size_t i0 = 0u; i0 <= nBins.at(0) + 1u; ++i0) {
      for (std::size_t i1 = 0u; i1 <= nBins.at(1) + 1u; ++i1) {
        nlohmann::json jBin;
        jBin[kLocalKey] = std::vector<std::size_t>{i0, i1};
        const auto* target = view.atLocalBins({i0, i1});
        if (target == nullptr) {
          jBin[kTargetVolumeIdKey] = nullptr;
        } else {
          jBin[kTargetVolumeIdKey] = volumeIds.at(*target);
        }
        jLink[kBinsKey].push_back(std::move(jBin));
      }
    }
  } else {
    throw std::invalid_argument("Unsupported GridPortalLink dimensionality");
  }

  jLink[kArtifactLinksKey] = nlohmann::json::array();
  for (const auto& artifact : link.artifactPortalLinks()) {
    jLink[kArtifactLinksKey].push_back(
        converter.portalLinkToJson(gctx, artifact, surfaceIds, volumeIds));
  }

  return jLink;
}

/// Decode a `TrivialPortalLink` payload.
std::unique_ptr<Acts::PortalLinkBase> decodeTrivialPortalLink(
    const nlohmann::json& encoded, const Acts::GeometryContext& /*gctx*/,
    const Acts::TrackingGeometryJsonConverter& /*converter*/,
    const Acts::TrackingGeometryJsonConverter::SurfacePointerLookup& surfaces,
    const Acts::TrackingGeometryJsonConverter::VolumePointerLookup& volumes) {
  const auto linkSurfaceId = encoded.at(kSurfaceIdKey).get<std::size_t>();
  const auto targetVolumeId = encoded.at(kTargetVolumeIdKey).get<std::size_t>();
  return std::make_unique<Acts::TrivialPortalLink>(surfaces.at(linkSurfaceId),
                                                   *volumes.at(targetVolumeId));
}

/// Decode a `CompositePortalLink` payload with recursive child decoding.
std::unique_ptr<Acts::PortalLinkBase> decodeCompositePortalLink(
    const nlohmann::json& encoded, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& converter,
    const Acts::TrackingGeometryJsonConverter::SurfacePointerLookup& surfaces,
    const Acts::TrackingGeometryJsonConverter::VolumePointerLookup& volumes) {
  const auto direction = encoded.at(kDirectionKey).get<Acts::AxisDirection>();
  std::vector<std::unique_ptr<Acts::PortalLinkBase>> children;
  for (const auto& child : encoded.at(kChildrenKey)) {
    children.push_back(
        converter.portalLinkFromJson(gctx, child, surfaces, volumes));
  }
  return std::make_unique<Acts::CompositePortalLink>(std::move(children),
                                                     direction);
}

/// Decode a `GridPortalLink` payload including axis/bin assignment and
/// artifacts.
std::unique_ptr<Acts::PortalLinkBase> decodeGridPortalLink(
    const nlohmann::json& encoded, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& converter,
    const Acts::TrackingGeometryJsonConverter::SurfacePointerLookup& surfaces,
    const Acts::TrackingGeometryJsonConverter::VolumePointerLookup& volumes) {
  auto linkSurfaceId = encoded.at(kSurfaceIdKey).get<std::size_t>();
  const auto direction = encoded.at(kDirectionKey).get<Acts::AxisDirection>();

  std::vector<std::unique_ptr<Acts::IAxis>> axes;
  for (const auto& jAxis : encoded.at(kAxesKey)) {
    axes.push_back(axisFromJson(jAxis));
  }
  if (axes.empty() || axes.size() > 2u) {
    throw std::invalid_argument("GridPortalLink requires 1 or 2 axes");
  }

  auto grid =
      makeGridPortalLink(surfaces.at(linkSurfaceId), direction, *axes.at(0),
                         axes.size() == 2u ? axes.at(1).get() : nullptr);

  Acts::AnyGridView<const Acts::TrackingVolume*> view(grid->grid());
  const auto dim = view.dimensions();
  for (const auto& jBin : encoded.at(kBinsKey)) {
    auto local = jBin.at(kLocalKey).get<std::vector<std::size_t>>();
    if (local.size() != dim) {
      throw std::invalid_argument("Grid bin dimensionality mismatch");
    }
    Acts::IGrid::AnyIndexType localIndices(local.begin(), local.end());

    const Acts::TrackingVolume* target = nullptr;
    if (!jBin.at(kTargetVolumeIdKey).is_null()) {
      const auto targetVolumeId =
          jBin.at(kTargetVolumeIdKey).get<std::size_t>();
      target = volumes.at(targetVolumeId);
    }
    view.atLocalBins(localIndices) = target;
  }

  std::vector<Acts::TrivialPortalLink> artifacts;
  if (encoded.contains(kArtifactLinksKey)) {
    for (const auto& jArtifact : encoded.at(kArtifactLinksKey)) {
      auto decodedArtifact =
          converter.portalLinkFromJson(gctx, jArtifact, surfaces, volumes);
      auto* trivial =
          dynamic_cast<Acts::TrivialPortalLink*>(decodedArtifact.get());
      if (trivial == nullptr) {
        throw std::invalid_argument(
            "GridPortalLink artifact link is not trivial");
      }
      artifacts.push_back(std::move(*trivial));
    }
  }
  grid->setArtifactPortalLinks(std::move(artifacts));

  return grid;
}

/// Temporary decoded representation of one serialized surface entry.
struct SurfaceRecord {
  std::size_t surfaceId = 0u;
  nlohmann::json payload;
};

/// Temporary decoded representation of one serialized portal entry.
struct PortalRecord {
  std::size_t portalId = 0u;
  nlohmann::json payload;
};

/// Temporary decoded representation of one serialized volume entry.
struct VolumeRecord {
  std::size_t volumeId = 0u;
  std::string name;
  Acts::GeometryIdentifier::Value geometryId = 0u;
  Acts::Transform3 transform{Acts::Transform3::Identity()};
  nlohmann::json bounds;
  std::vector<std::size_t> children;
  std::vector<std::size_t> portalIds;
  std::vector<std::size_t> surfaceIds;
  nlohmann::json navigationPolicy;
};

/// Validate top-level schema metadata for tracking-geometry JSON payload.
void verifySchemaHeader(const nlohmann::json& encoded) {
  if (!encoded.contains(kHeaderKey)) {
    throw std::invalid_argument("Missing geometry JSON header");
  }
  const auto& header = encoded.at(kHeaderKey);
  if (header.at(kVersionKey).get<int>() != kFormatVersion) {
    throw std::invalid_argument("Unsupported geometry JSON format version");
  }
  if (header.at(kScopeKey).get<std::string>() != kScopeValue) {
    throw std::invalid_argument("Unexpected geometry JSON scope");
  }
}

void collectGeometry(
    const Acts::TrackingVolume& volume,
    std::vector<const Acts::Surface*>& orderedSurfaces,
    std::vector<const Acts::Portal*>& orderedPortals,
    std::vector<const Acts::TrackingVolume*>& orderedVolumes,
    Acts::TrackingGeometryJsonConverter::SurfaceIdLookup& surfaceIds,
    Acts::TrackingGeometryJsonConverter::PortalIdLookup& portalIds,
    Acts::TrackingGeometryJsonConverter::VolumeIdLookup& volumeIds) {
  auto insertSurface = [&](const auto& surf) {
    if (surfaceIds.emplace(surf, orderedSurfaces.size())) {
      orderedSurfaces.push_back(&surf);
    }
  };
  auto insertPortal = [&](const auto& port) {
    if (portalIds.emplace(port, orderedPortals.size())) {
      orderedPortals.push_back(&port);
    }
  };
  auto insertVolume = [&](const auto& vol) {
    if (volumeIds.emplace(vol, orderedVolumes.size())) {
      orderedVolumes.push_back(&vol);
    }
  };

  insertVolume(volume);

  for (const auto& portal : volume.portals()) {
    insertPortal(portal);

    insertSurface(portal.surface());

    const auto* along = portal.getLink(Acts::Direction::AlongNormal());
    if (const auto* alongComposite =
            dynamic_cast<const Acts::CompositePortalLink*>(along)) {
      for (const auto& link : alongComposite->links()) {
        insertSurface(link.surface());
      }
    }

    const auto* opposite = portal.getLink(Acts::Direction::OppositeNormal());
    if (const auto* oppositeComposite =
            dynamic_cast<const Acts::CompositePortalLink*>(opposite)) {
      for (const auto& link : oppositeComposite->links()) {
        insertSurface(link.surface());
      }
    }
  }
  for (const auto& surf : volume.surfaces()) {
    insertSurface(surf);
  }

  for (const auto& child : volume.volumes()) {
    collectGeometry(child, orderedSurfaces, orderedPortals, orderedVolumes,
                    surfaceIds, portalIds, volumeIds);
  }
}

/// Ensure stable geometry identifiers exist for volumes, boundaries and
/// portal surfaces before emitting `TrackingGeometry`.
void ensureIdentifiers(Acts::TrackingVolume& volume,
                       Acts::GeometryIdentifier::Value& nextVolumeId) {
  Acts::GeometryIdentifier volumeId = volume.geometryId();
  if (volumeId == Acts::GeometryIdentifier{}) {
    volumeId = Acts::GeometryIdentifier{}.withVolume(nextVolumeId++);
    volume.assignGeometryId(volumeId);
  }

  for (const auto [ib, boundary] : Acts::enumerate(volume.boundarySurfaces())) {
    auto& mutableBoundarySurface =
        const_cast<Acts::RegularSurface&>(boundary->surfaceRepresentation());
    if (mutableBoundarySurface.geometryId() == Acts::GeometryIdentifier{}) {
      mutableBoundarySurface.assignGeometryId(
          Acts::GeometryIdentifier(volumeId).withBoundary(ib + 1u));
    }
  }

  std::size_t portalIndex = 0u;
  for (auto& portal : volume.portals()) {
    auto& mutablePortalSurface = portal.surface();
    if (mutablePortalSurface.geometryId() == Acts::GeometryIdentifier{}) {
      mutablePortalSurface.assignGeometryId(
          Acts::GeometryIdentifier(volumeId).withExtra(portalIndex + 1u));
    }
    ++portalIndex;
  }

  for (auto& child : volume.volumes()) {
    ensureIdentifiers(child, nextVolumeId);
  }
}

}  // namespace

/// Build default converter configuration with all supported bounds and portal
/// link encoders/decoders registered.
Acts::TrackingGeometryJsonConverter::Config
Acts::TrackingGeometryJsonConverter::Config::defaultConfig() {
  Config cfg;

  cfg.encodeVolumeBounds.registerFunction(encodeConeVolumeBounds)
      .registerFunction(encodeCuboidVolumeBounds)
      .registerFunction(encodeCutoutCylinderVolumeBounds)
      .registerFunction(encodeCylinderVolumeBounds)
      .registerFunction(encodeGenericCuboidVolumeBounds)
      .registerFunction(encodeTrapezoidVolumeBounds);

  cfg.encodeNavigationPolicy.registerFunction(encodeTryAllNavigationPolicy)
      .registerFunction(encodeSurfaceArrayNavigationPolicy)
      .registerFunction(encodeMultiNavigationPolicy);

  cfg.encodePortalLink.registerFunction(encodeTrivialPortalLink)
      .registerFunction(encodeCompositePortalLink)
      .registerFunction(encodeGridPortalLink);

  cfg.decodePortalLink
      .registerKind(kTrivialPortalLinkKind, decodeTrivialPortalLink)
      .registerKind(kCompositePortalLinkKind, decodeCompositePortalLink)
      .registerKind(kGridPortalLinkKind, decodeGridPortalLink);

  cfg.decodeVolumeBounds
      .registerKind(kConeVolumeBoundsKind,
                    decodeVolumeBoundsT<ConeVolumeBounds>)
      .registerKind(kCuboidVolumeBoundsKind,
                    decodeVolumeBoundsT<CuboidVolumeBounds>)
      .registerKind(kCutoutCylinderVolumeBoundsKind,
                    decodeVolumeBoundsT<CutoutCylinderVolumeBounds>)
      .registerKind(kCylinderVolumeBoundsKind,
                    decodeVolumeBoundsT<CylinderVolumeBounds>)
      .registerKind(kGenericCuboidVolumeBoundsKind,
                    decodeVolumeBoundsT<GenericCuboidVolumeBounds>)
      .registerKind(kTrapezoidVolumeBoundsKind,
                    decodeVolumeBoundsT<TrapezoidVolumeBounds>);

  cfg.decodeNavigationPolicy
      .registerKind(kTryAllNavigationPolicyKind, decodeTryAllNavigationPolicy)
      .registerKind(kSurfaceArrayNavigationPolicyKind,
                    decodeSurfaceArrayNavigationPolicy)
      .registerKind(kMultiNavigationPolicyKind, decodeMultiNavigationPolicy);

  return cfg;
}

/// Construct converter from a configuration object.
Acts::TrackingGeometryJsonConverter::TrackingGeometryJsonConverter(
    Config config)
    : m_cfg(std::move(config)) {}

/// Serialize a full `TrackingGeometry` by delegating to world volume
/// serialization.
nlohmann::json Acts::TrackingGeometryJsonConverter::toJson(
    const GeometryContext& gctx, const TrackingGeometry& geometry,
    const Options& options) const {
  return toJson(gctx, *geometry.highestTrackingVolume(), options);
}

/// Serialize a single portal link using the configured type dispatcher.
nlohmann::json Acts::TrackingGeometryJsonConverter::portalLinkToJson(
    const GeometryContext& gctx, const PortalLinkBase& link,
    const SurfaceIdLookup& surfaceIds, const VolumeIdLookup& volumeIds) const {
  return m_cfg.encodePortalLink(link, gctx, *this, surfaceIds, volumeIds);
}

/// Deserialize a single portal link using the configured kind dispatcher.
std::unique_ptr<Acts::PortalLinkBase>
Acts::TrackingGeometryJsonConverter::portalLinkFromJson(
    const GeometryContext& gctx, const nlohmann::json& encoded,
    const SurfacePointerLookup& surfaces,
    const VolumePointerLookup& volumes) const {
  return m_cfg.decodePortalLink(encoded, gctx, *this, surfaces, volumes);
}

nlohmann::json Acts::TrackingGeometryJsonConverter::navigationPolicyToJson(
    const Acts::INavigationPolicy& policy) const {
  return m_cfg.encodeNavigationPolicy(policy, *this);
}

std::unique_ptr<Acts::INavigationPolicy>
Acts::TrackingGeometryJsonConverter::navigationPolicyFromJson(
    const Acts::GeometryContext& gctx, const nlohmann::json& encoded,
    const Acts::TrackingVolume& volume, const Acts::Logger& logger) const {
  return m_cfg.decodeNavigationPolicy(encoded, gctx, *this, volume, logger);
}

/// Serialize one world volume hierarchy to JSON.
///
/// This writes schema metadata, emits all reachable volumes in deterministic
/// depth-first order, stores per-volume portal references, and writes unique
/// portal payloads/links with explicit kind tags in a top-level portal table.
nlohmann::json Acts::TrackingGeometryJsonConverter::toJson(
    const GeometryContext& gctx, const TrackingVolume& world,
    const Options& /*options*/) const {
  nlohmann::json encoded;
  encoded[kHeaderKey] = nlohmann::json::object();
  encoded[kHeaderKey][kVersionKey] = kFormatVersion;
  encoded[kHeaderKey][kScopeKey] = kScopeValue;

  // Collect object
  std::vector<const Surface*> orderedSurfaces;
  std::vector<const Portal*> orderedPortals;
  std::vector<const TrackingVolume*> orderedVolumes;

  SurfaceIdLookup surfaceIds;
  PortalIdLookup portalIds;
  VolumeIdLookup volumeIds;
  collectGeometry(world, orderedSurfaces, orderedPortals, orderedVolumes,
                  surfaceIds, portalIds, volumeIds);

  encoded[kSurfacesKey] = nlohmann::json::array();
  encoded[kPortalsKey] = nlohmann::json::array();
  encoded[kVolumesKey] = nlohmann::json::array();

  encoded[kRootVolumeIdKey] = volumeIds.at(world);

  // Encode surfaces
  // ---------------------------------------------------
  // SurfaceJsonConverter::Options surfaceConverterOpt;
  // surfaceConverterOpt.placementEncoder = m_cfg.surfacePlacementEncoder;
  // ---------------------------------------------------
  for (const auto* surf : orderedSurfaces) {
    nlohmann::json jSurface = SurfaceJsonConverter::toJson(gctx, *surf);
    jSurface[kSurfaceIdKey] = surfaceIds.at(*surf);
    encoded[kSurfacesKey].push_back(std::move(jSurface));
  }

  // Encode portals
  for (const auto* portal : orderedPortals) {
    nlohmann::json jPortal;
    jPortal[kPortalIdKey] = portalIds.at(*portal);

    if (const auto* along = portal->getLink(Direction::AlongNormal());
        along != nullptr) {
      jPortal[kAlongNormalKey] =
          portalLinkToJson(gctx, *along, surfaceIds, volumeIds);
    } else {
      jPortal[kAlongNormalKey] = nullptr;
    }

    if (const auto* opposite = portal->getLink(Direction::OppositeNormal());
        opposite != nullptr) {
      jPortal[kOppositeNormalKey] =
          portalLinkToJson(gctx, *opposite, surfaceIds, volumeIds);
    } else {
      jPortal[kOppositeNormalKey] = nullptr;
    }

    jPortal[kSurfaceIdKey] = surfaceIds.at(portal->surface());
    encoded[kPortalsKey].push_back(std::move(jPortal));
  }

  // Encode volumes
  for (const auto* volume : orderedVolumes) {
    nlohmann::json jVolume;
    jVolume[kVolumeIdKey] = volumeIds.at(*volume);
    jVolume[kNameKey] = volume->volumeName();
    jVolume[kGeometryIdKey] = nlohmann::json(volume->geometryId());
    jVolume[kTransformKey] =
        Transform3JsonConverter::toJson(volume->localToGlobalTransform(gctx));
    jVolume[kBoundsKey] = m_cfg.encodeVolumeBounds(volume->volumeBounds());
    jVolume[kNavigationPolicyKey] =
        navigationPolicyToJson(*volume->navigationPolicy());

    jVolume[kChildrenKey] = nlohmann::json::array();
    for (const auto& child : volume->volumes()) {
      jVolume[kChildrenKey].push_back(volumeIds.at(child));
    }

    jVolume[kPortalIdsKey] = nlohmann::json::array();
    for (const auto& portal : volume->portals()) {
      jVolume[kPortalIdsKey].push_back(portalIds.at(portal));
    }

    jVolume[kSurfaceIdKey] = nlohmann::json::array();
    for (const auto& surface : volume->surfaces()) {
      jVolume[kSurfaceIdKey].push_back(surfaceIds.at(surface));
    }

    encoded[kVolumesKey].push_back(std::move(jVolume));
  }

  return encoded;
}

/// Deserialize one world `TrackingVolume` hierarchy from JSON.
///
/// The reconstruction is performed in phases: validate header, decode all
/// volume records, instantiate volumes, attach child hierarchy, then decode and
/// attach shared portals by ID.
std::shared_ptr<Acts::TrackingVolume>
Acts::TrackingGeometryJsonConverter::trackingVolumeFromJson(
    const GeometryContext& gctx, const nlohmann::json& encoded,
    const Options& /*options*/) {
  verifySchemaHeader(encoded);

  if (!encoded.contains(kVolumesKey) || !encoded.contains(kRootVolumeIdKey) ||
      !encoded.contains(kPortalsKey)) {
    throw std::invalid_argument(
        "Missing volume payload in tracking geometry JSON");
  }

  // Collect surface data
  std::unordered_map<std::size_t, SurfaceRecord> surfaceRecords;
  for (const auto& jSurface : encoded.at(kSurfacesKey)) {
    SurfaceRecord record;
    record.surfaceId = jSurface.at(kSurfaceIdKey).get<std::size_t>();
    record.payload = jSurface;
    const auto inserted =
        surfaceRecords.emplace(record.surfaceId, std::move(record));
    if (!inserted.second) {
      throw std::invalid_argument("Duplicate serialized surface ID");
    }
  }

  // Collect portal data
  std::unordered_map<std::size_t, PortalRecord> portalRecords;
  for (const auto& jPortal : encoded.at(kPortalsKey)) {
    PortalRecord record;
    record.portalId = jPortal.at(kPortalIdKey).get<std::size_t>();
    record.payload = jPortal;
    const auto inserted =
        portalRecords.emplace(record.portalId, std::move(record));
    if (!inserted.second) {
      throw std::invalid_argument("Duplicate serialized portal ID");
    }
  }

  // Collect volume data
  std::unordered_map<std::size_t, VolumeRecord> volumeRecords;
  for (const auto& jVolume : encoded.at(kVolumesKey)) {
    VolumeRecord record;
    record.volumeId = jVolume.at(kVolumeIdKey).get<std::size_t>();
    record.name = jVolume.at(kNameKey).get<std::string>();
    record.transform =
        Transform3JsonConverter::fromJson(jVolume.at(kTransformKey));
    record.bounds = jVolume.at(kBoundsKey);
    record.children = jVolume.value(kChildrenKey, std::vector<std::size_t>{});
    record.portalIds = jVolume.value(kPortalIdsKey, std::vector<std::size_t>{});
    record.surfaceIds =
        jVolume.value(kSurfaceIdKey, std::vector<std::size_t>{});
    record.navigationPolicy = jVolume.at(kNavigationPolicyKey);

    if (!jVolume["geometry_id"].is_null()) {
      GeometryIdentifier geoID =
          jVolume["geometry_id"].get<GeometryIdentifier>();
      record.geometryId = geoID.value();
    } else {
      record.geometryId = 0;
    }

    const auto inserted =
        volumeRecords.emplace(record.volumeId, std::move(record));
    if (!inserted.second) {
      throw std::invalid_argument("Duplicate serialized volume ID");
    }
  }

  // Get root volume id
  const std::size_t rootVolumeId =
      encoded.at(kRootVolumeIdKey).get<std::size_t>();
  if (!volumeRecords.contains(rootVolumeId)) {
    throw std::invalid_argument("Serialized root volume ID does not exist");
  }

  // Collect surface pointers
  SurfacePointerLookup surfacePointers;

  // ---------------------------------------------------
  // SurfaceJsonConverter::Options surfaceConverterOpt;
  // surfaceConverterOpt.placementDecoder = m_cfg.surfacePlacementDecoder;
  // ---------------------------------------------------
  for (const auto& [surfaceId, record] : surfaceRecords) {
    auto surface = regularSurfaceFromJson(record.payload);
    surfacePointers.emplace(surfaceId, surface);
  }

  // Collect volume pointers
  std::unordered_map<std::size_t, std::unique_ptr<TrackingVolume>>
      volumeStorage;
  VolumePointerLookup volumePointers;

  for (const auto& [volumeId, record] : volumeRecords) {
    auto volumeBounds = m_cfg.decodeVolumeBounds(record.bounds);
    auto volume = std::make_unique<TrackingVolume>(
        record.transform, std::move(volumeBounds), record.name);

    GeometryIdentifier geometryId(record.geometryId);
    if (geometryId == GeometryIdentifier{}) {
      geometryId = GeometryIdentifier{}.withVolume(volumeId + 1u);
    }
    volume->assignGeometryId(geometryId);
    volumePointers.emplace(volumeId, volume.get());
    volumeStorage.emplace(volumeId, std::move(volume));
  }

  std::unordered_set<std::size_t> visiting;
  std::unordered_set<std::size_t> built;

  // Assemble the volume hierarchy
  std::function<void(std::size_t)> attachChildren =
      [&](std::size_t volumeId) -> void {
    if (built.contains(volumeId)) {
      return;
    }
    if (!visiting.insert(volumeId).second) {
      throw std::invalid_argument(
          "Cycle detected in serialized volume hierarchy");
    }

    auto& parent = volumeStorage.at(volumeId);
    if (parent == nullptr) {
      throw std::invalid_argument("Volume was already moved unexpectedly");
    }

    for (std::size_t childId : volumeRecords.at(volumeId).children) {
      attachChildren(childId);
      auto& child = volumeStorage.at(childId);
      if (child == nullptr) {
        throw std::invalid_argument(
            "Serialized child volume has already been attached");
      }
      parent->addVolume(std::move(child));
    }

    visiting.erase(volumeId);
    built.insert(volumeId);
  };

  attachChildren(rootVolumeId);

  auto decodePortal = [&](const nlohmann::json& jPortal) {
    std::unique_ptr<PortalLinkBase> along = nullptr;
    std::unique_ptr<PortalLinkBase> opposite = nullptr;

    if (!jPortal.at(kAlongNormalKey).is_null()) {
      along = portalLinkFromJson(gctx, jPortal.at(kAlongNormalKey),
                                 surfacePointers, volumePointers);
    }
    if (!jPortal.at(kOppositeNormalKey).is_null()) {
      opposite = portalLinkFromJson(gctx, jPortal.at(kOppositeNormalKey),
                                    surfacePointers, volumePointers);
    }
    if (along == nullptr && opposite == nullptr) {
      throw std::invalid_argument("Portal has no links");
    }

    auto portal =
        std::make_shared<Portal>(gctx, std::move(along), std::move(opposite));
    portal->surface().assignGeometryId(
        surfacePointers.at(jPortal.at(kSurfaceIdKey))->geometryId());
    return portal;
  };

  PortalPointerLookup portalPointers;
  for (const auto& [portalId, record] : portalRecords) {
    const auto inserted =
        portalPointers.emplace(portalId, decodePortal(record.payload));
    if (!inserted) {
      throw std::invalid_argument("Portal pointer reconstruction failed");
    }
  }

  auto logger =
      Acts::getDefaultLogger("navigationPolicyLogger", Acts::Logging::VERBOSE);
  for (const auto& [volumeId, record] : volumeRecords) {
    auto* volume = volumePointers.find(volumeId);
    if (volume == nullptr) {
      throw std::invalid_argument("Volume pointer reconstruction failed");
    }

    for (const std::size_t portalId : record.portalIds) {
      volume->addPortal(portalPointers.at(portalId));
    }
    for (const std::size_t surfaceId : record.surfaceIds) {
      std::cout << "ADDING SURFACE "
                << surfacePointers.at(surfaceId)->geometryId()
                << " WITH PLACEMENT "
                << (surfacePointers.at(surfaceId)->surfacePlacement() ==
                    nullptr)
                << "\n";
      volume->addSurface(surfacePointers.at(surfaceId));
    }

    volume->setNavigationPolicy(navigationPolicyFromJson(
        gctx, record.navigationPolicy, *volume, *logger));
  }

  auto root = std::move(volumeStorage.at(rootVolumeId));
  if (root == nullptr) {
    throw std::invalid_argument("Root volume reconstruction failed");
  }

  return std::shared_ptr<TrackingVolume>(std::move(root));
}

/// Deserialize a full `TrackingGeometry` from serialized world-volume payload.
///
/// After world volume reconstruction this also ensures missing geometry
/// identifiers are assigned consistently.
std::shared_ptr<Acts::TrackingGeometry>
Acts::TrackingGeometryJsonConverter::trackingGeometryFromJson(
    const GeometryContext& gctx, const nlohmann::json& encoded,
    const Options& options) {
  auto world = trackingVolumeFromJson(gctx, encoded, options);

  GeometryIdentifier::Value nextVolumeId = 1u;
  ensureIdentifiers(*world, nextVolumeId);

  return std::make_shared<TrackingGeometry>(
      world, nullptr, GeometryIdentifierHook{}, getDummyLogger(), false);
}
