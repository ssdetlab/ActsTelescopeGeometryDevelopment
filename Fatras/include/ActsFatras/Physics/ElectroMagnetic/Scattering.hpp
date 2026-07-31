// This file is part of the Acts project.
//
// Copyright (C) 2018-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/detail/GaussianMixture.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/detail/GeneralMixture.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/detail/Highland.hpp"

#include <array>
#include <random>

namespace ActsFatras {

/// Simulate (multiple) scattering using a configurable scattering model.
///
/// @tparam scattering_model_t Model implementation to draw a scattering angle.
template <typename scattering_model_t>
struct GenericScattering {
  /// The scattering formula
  scattering_model_t angle;

  /// Simulate scattering and update the particle parameters.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      defines the passed material
  /// @param[in,out] particle  is the particle being updated
  /// @return Empty secondaries containers.
  ///
  /// @tparam generator_t is a RandomNumberEngine
  template <typename generator_t>
  std::array<Particle, 0> operator()(generator_t &generator,
                                     const Acts::MaterialSlab &slab,
                                     Particle &particle) const {
    // Get the scattering angle
    const auto theta0 = Acts::computeMultipleScatteringTheta0(
        slab, particle.absolutePdg(), particle.mass(), particle.qOverP(),
        particle.absoluteCharge());

    // Draw to independent plane-projected normal angles
    double dx = std::normal_distribution<double>(0, 1)(generator) * theta0;
    double dy = std::normal_distribution<double>(0, 1)(generator) * theta0;
    double theta = std::hypot(dx, dy);
    double psi = atan2(dy, dx);

    Acts::Vector3 direction = particle.direction();
    // construct the combined rotation to the scattered direction
    Acts::RotationMatrix3 rotation(
        // rotation of the scattering deflector axis relative to the reference
        Acts::AngleAxis3(psi, direction) *
        // rotation by the scattering angle around the deflector axis
        Acts::AngleAxis3(theta, Acts::makeCurvilinearUnitU(direction)));
    direction.applyOnTheLeft(rotation);
    particle.setDirection(direction);

    // scattering is non-destructive and produces no secondaries
    return {};
  }
};

using GaussianMixtureScattering = GenericScattering<detail::GaussianMixture>;
using GeneralMixtureScattering = GenericScattering<detail::GeneralMixture>;
using HighlandScattering = GenericScattering<detail::Highland>;

}  // namespace ActsFatras
