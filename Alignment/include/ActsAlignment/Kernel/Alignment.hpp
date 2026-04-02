// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"

#include <limits>
#include <unordered_map>
#include <vector>

namespace ActsAlignment {

/// @brief Alignment transform updater
using AlignmentTransformUpdater = Acts::Delegate<bool(
    Acts::DetectorElementBase*, const Acts::GeometryContext&,
    const Acts::Vector3&, const Acts::Vector3&)>;

/// @brief Alignment result struct
struct AlignmentResult {
  // The change of alignment parameters
  Acts::ActsDynamicVector deltaAlignmentParameters;
  // The aligned parameters for detector elements
  std::unordered_map<Acts::DetectorElementBase*, Acts::Transform3>
      alignedParameters;
  // The covariance of alignment parameters
  Acts::ActsDynamicMatrix alignmentCovariance;
  // The average chi2/ndf (ndf is the measurement dim)
  double averageChi2ONdf = std::numeric_limits<double>::max();
  // The delta chi2
  double deltaChi2 = std::numeric_limits<double>::max();
  // The chi2
  double chi2 = 0;
  // The measurement dimension from all tracks
  std::size_t measurementDim = 0;
  // The alignment degree of freedom
  std::size_t alignmentDof = 0;
  // The number of tracks used for alignment
  std::size_t numTracks = 0;
  // The indexed alignable surfaces
  std::unordered_map<const Acts::Surface*, std::size_t> idxedAlignSurfaces;
};

/// @brief Alignment parameters solver
using AlignmentParametersSolver = Acts::Delegate<void(
    const Acts::GeometryContext& gctx, AlignmentResult& alignRes,
    const Acts::ActsDynamicVector& sumChi2Derivative,
    const Acts::ActsDynamicMatrix& sumChi2SecondDerivative)>;

/// @brief Options for align() call
///
/// @tparam fit_options_t The fit options type
template <typename fit_options_t>
struct AlignmentOptions {
  /// Deleted default constructor
  AlignmentOptions() = delete;

  /// AlignmentOptions
  ///
  /// @param fOptions KF fit options
  /// @param aTransformUpdater updater for the alignment transform
  /// @param aParametersSolver solver for the alignment parameters delta
  /// @param aDetElements alignable detector elements
  /// @param chi2CufOff chi2 alignment convergence tolerance
  /// @param deltaChi2CutOff convergence threshold of chi2/ndf/nIt
  /// @param maxIters maximum number of iteration of the alignment fit
  AlignmentOptions(const fit_options_t& fOptions,
                   const AlignmentTransformUpdater& aTransformUpdater,
                   const AlignmentParametersSolver& aParametersSolver,
                   const std::vector<Acts::DetectorElementBase*>& aDetElements,
                   double chi2CutOff,
                   const std::pair<std::size_t, double>& deltaChi2CutOff,
                   std::size_t maxIters, AlignmentMask mask)
      : fitOptions(fOptions),
        alignmentTransformUpdater(aTransformUpdater),
        alignmentParametersSolver(aParametersSolver),
        alignedDetElements(aDetElements),
        averageChi2ONdfCutOff(chi2CutOff),
        deltaAverageChi2ONdfCutOff(deltaChi2CutOff),
        maxIterations(maxIters),
        alignmentMask(mask) {}

  /// Fit options
  fit_options_t fitOptions;
  /// Alignment transform updater
  AlignmentTransformUpdater alignmentTransformUpdater;
  /// Alignment parameters solver
  AlignmentParametersSolver alignmentParametersSolver;
  /// The detector elements to be aligned
  std::vector<Acts::DetectorElementBase*> alignedDetElements;
  /// The alignment tolerance to determine if the alignment is covered
  double averageChi2ONdfCutOff;
  /// The therhold of average chi2/ndf delta within a number of iterations to
  /// determine if alignment is converged
  std::pair<std::size_t, double> deltaAverageChi2ONdfCutOff;
  /// The maximum number of iterations to run alignment
  std::size_t maxIterations;
  /// The alignment mask
  AlignmentMask alignmentMask;
};

/// @brief KalmanFitter-based alignment implementation
///
/// @tparam fitter_t Type of the fitter class
template <typename fitter_t>
struct Alignment {
  /// Default constructor is deleted
  Alignment() = delete;

  /// Constructor from arguments
  Alignment(fitter_t fitter,
            std::unique_ptr<const Acts::Logger> _logger =
                Acts::getDefaultLogger("Alignment", Acts::Logging::INFO))
      : m_fitter(std::move(fitter)), m_logger{std::move(_logger)} {}

  /// @brief evaluate alignment state for a single track
  ///
  /// @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam fit_options_t The fit options type
  ///
  /// @param gctx The current geometry context object
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param fitOptions The fit Options steering the fit
  /// @param idxedAlignSurfaces The idxed surfaces to be aligned
  /// @param alignMask The alignment mask (same for all detector element for the
  /// moment)
  ///
  /// @result The alignment state for a single track
  template <typename source_link_t, typename start_parameters_t,
            typename fit_options_t>
  Acts::Result<detail::TrackAlignmentState> evaluateTrackAlignmentState(
      const Acts::GeometryContext& gctx,
      const std::vector<source_link_t>& sourcelinks,
      const start_parameters_t& sParameters, const fit_options_t& fitOptions,
      const std::unordered_map<const Acts::Surface*, std::size_t>&
          idxedAlignSurfaces,
      const AlignmentMask& alignMask) const;

  /// @brief calculate the alignment parameters delta
  ///
  /// @tparam trajectory_container_t The trajectories container type
  /// @tparam start_parameters_t The initial parameters container type
  /// @tparam fit_options_t The fit options type
  ///
  /// @param trajectoryCollection The collection of trajectories as input of
  /// fitting
  /// @param startParametersCollection The collection of starting parameters as
  /// input of fitting
  /// @param fitOptions The fit Options steering the fit
  /// @param alignResult [in, out] The aligned result
  /// @param alignMask The alignment mask (same for all measurements now)
  template <typename trajectory_container_t,
            typename start_parameters_container_t, typename fit_options_t>
  void calculateAlignmentParameters(
      const Acts::GeometryContext& gctx,
      const trajectory_container_t& trajectoryCollection,
      const start_parameters_container_t& startParametersCollection,
      const fit_options_t& fitOptions, AlignmentResult& alignResult,
      const AlignmentMask& alignMask,
      const ActsAlignment::AlignmentParametersSolver& alignmentParametersSolver)
      const;

  /// @brief update the detector element alignment parameters
  ///
  /// @param gctx The geometry context
  /// @param alignedDetElements The detector elements to be aligned
  /// @param alignedTransformUpdater The updater for updating the aligned
  /// @param alignResult [in, out] The aligned result
  Acts::Result<void> updateAlignmentParameters(
      const Acts::GeometryContext& gctx,
      const std::vector<Acts::DetectorElementBase*>& alignedDetElements,
      const AlignmentTransformUpdater& alignedTransformUpdater,
      AlignmentResult& alignResult) const;

  /// @brief Alignment implementation
  ///
  /// @tparam trajectory_container_t The trajectories container type
  /// @tparam start_parameters_t The initial parameters container type
  /// @tparam fit_options_t The fit options type
  ///
  /// @param trajectoryCollection The collection of trajectories as input of
  /// fitting
  /// @param startParametersCollection The collection of starting parameters as
  /// input of fitting
  /// @param alignOptions The alignment options
  ///
  /// @result The alignment result
  template <typename trajectory_container_t,
            typename start_parameters_container_t, typename fit_options_t>
  Acts::Result<AlignmentResult> align(
      const trajectory_container_t& trajectoryCollection,
      const start_parameters_container_t& startParametersCollection,
      const AlignmentOptions<fit_options_t>& alignOptions) const;

 private:
  // The fitter
  fitter_t m_fitter;

  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};
}  // namespace ActsAlignment

#include "ActsAlignment/Kernel/Alignment.ipp"
