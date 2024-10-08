// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Detray/DetrayMaterialConverter.hpp"

#include "Acts/Detector/Detector.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Plugins/Detray/DetrayConversionUtils.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"

namespace {

struct MaterialSurfaceSelector {
  std::vector<const Acts::Surface*> surfaces = {};

  /// @param surface is the test surface
  void operator()(const Acts::Surface* surface) {
    if (surface->surfaceMaterial() != nullptr) {
      if (std::find(surfaces.begin(), surfaces.end(), surface) ==
          surfaces.end()) {
        surfaces.push_back(surface);
      }
    }
  }
};

/// This creates dummy axes to allow homogeneous material for the moment
/// to be represented as grid surface material
std::vector<detray::io::axis_payload> homogeneousAxesPayloads() {
  Acts::BinningData bDataX(Acts::BinningValue::binX, -1, 1);
  bDataX.option = Acts::BinningOption::closed;
  Acts::BinningData bDataY(Acts::BinningValue::binY, -1, 1);
  bDataY.option = Acts::BinningOption::closed;
  auto axisPayloadX = Acts::DetrayConversionUtils::convertBinningData(bDataX);
  auto axisPayloadY = Acts::DetrayConversionUtils::convertBinningData(bDataY);

  return {axisPayloadX, axisPayloadY};
}

}  // namespace

detray::io::material_slab_payload
Acts::DetrayMaterialConverter::convertMaterialSlab(
    const MaterialSlab& materialSlab) {
  detray::io::material_slab_payload slab;
  // Fill the material parameters and the thickness
  const auto& material = materialSlab.material();
  slab.thickness = materialSlab.thickness();
  slab.mat = detray::io::material_payload{
      {material.X0(), material.L0(), material.Ar(), material.Z(),
       material.massDensity(), material.molarDensity(), 0.}};
  slab.type = detray::io::material_id::slab;
  return slab;
}

detray::io::grid_payload<detray::io::material_slab_payload,
                         detray::io::material_id>
Acts::DetrayMaterialConverter::convertSurfaceMaterial(
    const ISurfaceMaterial& material, const Logger& logger) {
  detray::io::grid_payload<detray::io::material_slab_payload,
                           detray::io::material_id>
      materialGrid;

  // Check the material types
  // (1) homogeneous -> 1 x 1 bin grid with closed axes
  auto homogeneousMaterial =
      dynamic_cast<const HomogeneousSurfaceMaterial*>(&material);
  if (homogeneousMaterial != nullptr) {
    ACTS_VERBOSE(
        "DetrayMaterialConverter: found homogeneous surface material, this "
        "will be modelled as a 1x1 bin grid");
    // A single bin entry: convert it and fill it
    detray::io::material_slab_payload slab = convertMaterialSlab(
        homogeneousMaterial->materialSlab(Vector3{0., 0., 0.}));
    detray::io::grid_bin_payload<detray::io::material_slab_payload> slabBin{
        {0, 0}, {slab}};
    // Filling axes and bins
    materialGrid.axes = homogeneousAxesPayloads();
    materialGrid.bins = {slabBin};
    return materialGrid;
  }
  // (2) - binned material -> convert into grid structure
  auto binnedMaterial = dynamic_cast<const BinnedSurfaceMaterial*>(&material);
  if (binnedMaterial != nullptr) {
    ACTS_VERBOSE("DetrayMaterialConverter: found binned surface material");

    // BinUtility modifications
    bool swapped = false;
    // Get the bin utility (make a copy as we may modify it)
    // Detray expects 2-dimensional grid, currently supported are
    // x-y, r-phi, phi-z
    BinUtility bUtility = binnedMaterial->binUtility();
    // Turn the bin value into a 2D grid
    if (bUtility.dimensions() == 1u) {
      if (bUtility.binningData()[0u].binvalue == BinningValue::binR) {
        // Turn to R-Phi
        bUtility += BinUtility(1u, -M_PI, M_PI, closed, BinningValue::binR);
      } else if (bUtility.binningData()[0u].binvalue == BinningValue::binZ) {
        // Turn to Phi-Z - swap needed
        BinUtility nbUtility(1u, -M_PI, M_PI, closed, BinningValue::binPhi);
        nbUtility += bUtility;
        bUtility = std::move(nbUtility);
        swapped = true;
      } else {
        std::runtime_error("Unsupported binning for Detray");
      }
    } else if (bUtility.dimensions() == 2u &&
               bUtility.binningData()[0u].binvalue == BinningValue::binZ &&
               bUtility.binningData()[1u].binvalue == BinningValue::binPhi) {
      BinUtility nbUtility(bUtility.binningData()[1u]);
      nbUtility += bUtility.binningData()[0u];
      bUtility = std::move(nbUtility);
      swapped = true;
    }

    BinningValue bVal0 = bUtility.binningData()[0u].binvalue;
    BinningValue bVal1 = bUtility.binningData()[1u].binvalue;

    // Translate into grid index type
    detray::io::material_id gridIndexType = detray::io::material_id::unknown;
    if (bVal0 == BinningValue::binR && bVal1 == BinningValue::binPhi) {
      gridIndexType = detray::io::material_id::ring2_map;
    } else if (bVal0 == BinningValue::binPhi && bVal1 == BinningValue::binZ) {
      gridIndexType = detray::io::material_id::concentric_cylinder2_map;
    } else if (bVal0 == BinningValue::binX && bVal1 == BinningValue::binY) {
      gridIndexType = detray::io::material_id::rectangle2_map;
    } else {
      std::runtime_error("Unsupported binning for Detray");
    }

    detray::io::typed_link_payload<detray::io::material_id> linkPayload{
        gridIndexType, 0u};
    materialGrid.grid_link = linkPayload;

    // Now convert the (modified) bin utility
    for (const auto& bData : bUtility.binningData()) {
      auto axisPayload = DetrayConversionUtils::convertBinningData(bData);
      materialGrid.axes.push_back(axisPayload);
    }

    // Convert the material slabs from the material matrix
    auto materialMatrix = binnedMaterial->fullMaterial();
    for (std::size_t ib1 = 0; ib1 < materialMatrix.size(); ++ib1) {
      for (std::size_t ib0 = 0; ib0 < materialMatrix[0u].size(); ++ib0) {
        // Translate into a local bin
        std::size_t lb0 = swapped ? ib1 : ib0;
        std::size_t lb1 = swapped ? ib0 : ib1;
        detray::io::material_slab_payload slab =
            convertMaterialSlab(materialMatrix[ib1][ib0]);
        detray::io::grid_bin_payload<detray::io::material_slab_payload> slabBin{
            {static_cast<unsigned int>(lb0), static_cast<unsigned int>(lb1)},
            {slab}};
        // Fill into the grid
        materialGrid.bins.push_back(slabBin);
      }
    }
    return materialGrid;
  }

  throw std::invalid_argument(
      "DetrayMaterialConverter: unknown surface material type detected.");
}

detray::io::detector_grids_payload<detray::io::material_slab_payload,
                                   detray::io::material_id>
Acts::DetrayMaterialConverter::convertSurfaceMaterialGrids(
    const DetrayConversionUtils::GeometryIdCache& geoIdCache,
    const Experimental::Detector& detector, const Logger& logger) {
  // The material grid payload
  detray::io::detector_grids_payload<detray::io::material_slab_payload,
                                     detray::io::material_id>
      materialGrids;

  using DetrayMaterialGrid =
      detray::io::grid_payload<detray::io::material_slab_payload,
                               detray::io::material_id>;

  // Loop over the volumes in order to assign the right volume links
  for (const auto& volume : detector.volumes()) {
    // Per volume surface selector
    MaterialSurfaceSelector selector;
    volume->visitSurfaces(selector);
    ACTS_DEBUG("DetrayMaterialConverter: found "
               << selector.surfaces.size()
               << " surfaces/portals with material in volume "
               << volume->name());
    // Find the voluem index first
    auto volumeIndex = geoIdCache.volumeLinks.find(volume->geometryId());
    if (volumeIndex != geoIdCache.volumeLinks.end()) {
      std::vector<DetrayMaterialGrid> volumeMaterialGrids = {};
      // Now convert the surfaces
      for (const auto& surface : selector.surfaces) {
        // Find the surface index
        auto surfaceIndices =
            geoIdCache.localSurfaceLinks.equal_range(surface->geometryId());
        DetrayMaterialGrid materialGrid =
            convertSurfaceMaterial(*surface->surfaceMaterial(), logger);
        // Loop over the equal range and fill one grid each, this is needed
        // as the initial portal could be split into multiple surfaces
        for (auto itr = surfaceIndices.first; itr != surfaceIndices.second;
             ++itr) {
          // Fill the surface index
          materialGrid.owner_link =
              detray::io::single_link_payload{itr->second};
          // Fill the grid
          volumeMaterialGrids.push_back(materialGrid);
        }
      }
      // Register the grids of this volume
      materialGrids.grids.insert({volumeIndex->second, volumeMaterialGrids});

    } else {
      ACTS_WARNING(
          "DetrayMaterialConverter: volume not found in cache, should not "
          "happen.");
    }
  }
  // Return the material grids payload
  return materialGrids;
}
