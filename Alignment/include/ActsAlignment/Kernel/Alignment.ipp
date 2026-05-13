// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/TrackFitting/detail/KalmanGlobalCovariance.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/AlignmentError.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"

#include <cstddef>
#include <queue>

template <typename fitter_t>
template <typename source_link_container_t, typename start_parameters_t,
          typename fit_options_t>
Acts::Result<ActsAlignment::detail::TrackAlignmentState>
ActsAlignment::Alignment<fitter_t>::evaluateTrackAlignmentState(
    const Acts::GeometryContext& gctx,
    const source_link_container_t& sourceLinks,
    const start_parameters_t& sParameters, const fit_options_t& fitOptions,
    const std::unordered_map<const Acts::Surface*, std::size_t>&
        idxedAlignSurfaces,
    const AlignmentMask& alignMask) const {
  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  // Perform the fit
  auto fitRes = m_fitter.fit(sourceLinks.begin(), sourceLinks.end(),
                             sParameters, fitOptions, tracks);

  if (!fitRes.ok()) {
    ACTS_WARNING("Fit failure");
    return fitRes.error();
  }
  // The fit results
  const auto& track = fitRes.value();
  // Calculate the global track parameters covariance with the fitted track
  const auto& globalTrackParamsCov =
      Acts::detail::globalTrackParametersCovariance(
          tracks.trackStateContainer(), track.tipIndex());
  // Calculate the alignment state
  const auto alignState = detail::trackAlignmentState(
      gctx, tracks.trackStateContainer(), track.tipIndex(),
      globalTrackParamsCov, idxedAlignSurfaces, alignMask);
  if (alignState.alignmentDof == 0) {
    ACTS_VERBOSE("No alignment dof on track!");
    return AlignmentError::NoAlignmentDofOnTrack;
  }
  return alignState;
}

template <typename fitter_t>
template <typename trajectory_container_t,
          typename start_parameters_container_t, typename fit_options_t>
void ActsAlignment::Alignment<fitter_t>::calculateAlignmentParameters(
    const Acts::GeometryContext& gctx,
    const trajectory_container_t& trajectoryCollection,
    const start_parameters_container_t& startParametersCollection,
    const fit_options_t& fitOptions, AlignmentResult& alignResult,
    const AlignmentMask& alignMask,
    const ActsAlignment::AlignmentParametersSolver& alignmentParametersSolver)
    const {
  // The number of trajectories must be equal to the number of starting
  // parameters
  assert(trajectoryCollection.size() == startParametersCollection.size());

  // The total alignment degree of freedom
  alignResult.alignmentDof =
      alignResult.idxedAlignSurfaces.size() * Acts::eAlignmentSize;
  // Initialize derivative of chi2 w.r.t. alignment parameters for all tracks
  Acts::ActsDynamicVector sumChi2Derivative =
      Acts::ActsDynamicVector::Zero(alignResult.alignmentDof);
  Acts::ActsDynamicMatrix sumChi2SecondDerivative =
      Acts::ActsDynamicMatrix::Zero(alignResult.alignmentDof,
                                    alignResult.alignmentDof);
  // Copy the fit options
  fit_options_t fitOptionsWithRefSurface = fitOptions;
  // Calculate contribution to chi2 derivatives from all input trajectories
  alignResult.chi2 = 0;
  alignResult.measurementDim = 0;
  alignResult.numTracks = trajectoryCollection.size();
  double sumChi2ONdf = 0;
  for (unsigned int iTraj = 0; iTraj < trajectoryCollection.size(); iTraj++) {
    const auto& sourcelinks = trajectoryCollection.at(iTraj);
    const auto& sParameters = startParametersCollection.at(iTraj);
    auto evaluateRes = evaluateTrackAlignmentState(
        fitOptions.geoContext, sourcelinks, sParameters,
        fitOptionsWithRefSurface, alignResult.idxedAlignSurfaces, alignMask);
    if (!evaluateRes.ok()) {
      ACTS_DEBUG("Evaluation of alignment state for track " << iTraj
                                                            << " failed");
      continue;
    }
    const auto& alignState = evaluateRes.value();
    for (const auto& [rowSurface, rows] : alignState.alignedSurfaces) {
      const auto& [dstRow, srcRow] = rows;
      // Fill the results into full chi2 derivative matrix
      sumChi2Derivative.segment<Acts::eAlignmentSize>(dstRow *
                                                      Acts::eAlignmentSize) +=
          alignState.alignmentToChi2Derivative.segment(
              srcRow * Acts::eAlignmentSize, Acts::eAlignmentSize);

      for (const auto& [colSurface, cols] : alignState.alignedSurfaces) {
        const auto& [dstCol, srcCol] = cols;
        sumChi2SecondDerivative
            .block<Acts::eAlignmentSize, Acts::eAlignmentSize>(
                dstRow * Acts::eAlignmentSize, dstCol * Acts::eAlignmentSize) +=
            alignState.alignmentToChi2SecondDerivative.block(
                srcRow * Acts::eAlignmentSize, srcCol * Acts::eAlignmentSize,
                Acts::eAlignmentSize, Acts::eAlignmentSize);
      }
    }
    alignResult.chi2 += alignState.chi2;
    alignResult.measurementDim += alignState.measurementDim;
    sumChi2ONdf += alignState.chi2 / alignState.measurementDim;
  }
  alignResult.averageChi2ONdf = sumChi2ONdf / alignResult.numTracks;

  alignmentParametersSolver(gctx, alignResult, sumChi2Derivative,
                            sumChi2SecondDerivative);
}

template <typename fitter_t>
Acts::Result<void>
ActsAlignment::Alignment<fitter_t>::updateAlignmentParameters(
    const Acts::GeometryContext& gctx,
    const std::vector<Acts::DetectorElementBase*>& alignedDetElements,
    const ActsAlignment::AlignmentTransformUpdater& alignedTransformUpdater,
    ActsAlignment::AlignmentResult& alignResult) const {
  // Update the aligned transform
  Acts::AlignmentVector deltaAlignmentParam = Acts::AlignmentVector::Zero();
  for (const auto& [surface, index] : alignResult.idxedAlignSurfaces) {
    // Delta transform
    deltaAlignmentParam = alignResult.deltaAlignmentParameters.segment(
        Acts::eAlignmentSize * index, Acts::eAlignmentSize);
    // The delta translation
    Acts::Vector3 deltaCenter =
        deltaAlignmentParam.segment<3>(Acts::eAlignmentCenter0);
    // The delta Euler angles
    Acts::Vector3 deltaEulerAngles =
        deltaAlignmentParam.segment<3>(Acts::eAlignmentRotation0);

    // Update the aligned transform
    ACTS_VERBOSE("Delta of alignment parameters at element "
                 << index << "= \n"
                 << deltaAlignmentParam);
    bool updated = alignedTransformUpdater(alignedDetElements.at(index), gctx,
                                           deltaCenter, deltaEulerAngles);
    if (!updated) {
      ACTS_ERROR("Update alignment parameters for detector element failed");
      return AlignmentError::AlignmentParametersUpdateFailure;
    }
  }

  return Acts::Result<void>::success();
}

template <typename fitter_t>
template <typename trajectory_container_t,
          typename start_parameters_container_t, typename fit_options_t>
Acts::Result<ActsAlignment::AlignmentResult>
ActsAlignment::Alignment<fitter_t>::align(
    const trajectory_container_t& trajectoryCollection,
    const start_parameters_container_t& startParametersCollection,
    const ActsAlignment::AlignmentOptions<fit_options_t>& alignOptions) const {
  // Construct an AlignmentResult object
  AlignmentResult alignResult;

  // Assign index to the alignable surface
  for (unsigned int iDetElement = 0;
       iDetElement < alignOptions.alignedDetElements.size(); iDetElement++) {
    alignResult.idxedAlignSurfaces.emplace(
        &alignOptions.alignedDetElements.at(iDetElement)->surface(),
        iDetElement);
  }
  ACTS_VERBOSE("There are " << alignResult.idxedAlignSurfaces.size()
                            << " detector elements to be aligned");

  // Start the iteration to minimize the chi2
  bool converged = false;
  std::queue<double> recentChi2ONdf;

  // Perform the fit to the trajectories and update alignment parameters
  ACTS_INFO("Max number of iterations: " << alignOptions.maxIterations);
  for (unsigned int iIter = 0; iIter < alignOptions.maxIterations; iIter++) {
    // Calculate the alignment parameters delta etc.
    calculateAlignmentParameters(
        alignOptions.fitOptions.geoContext, trajectoryCollection,
        startParametersCollection, alignOptions.fitOptions, alignResult,
        alignOptions.alignmentMask, alignOptions.alignmentParametersSolver);
    // Screen out the information
    ACTS_INFO("iIter = " << iIter << ", total chi2 = " << alignResult.chi2
                         << ", total measurementDim = "
                         << alignResult.measurementDim
                         << " and average chi2/ndf = "
                         << alignResult.averageChi2ONdf);
    // Check if it has converged against the provided precision
    // 1. either the delta average chi2/ndf in the last few
    // iterations is within tolerance
    if (recentChi2ONdf.size() >=
        alignOptions.deltaAverageChi2ONdfCutOff.first) {
      if (std::abs(recentChi2ONdf.front() - alignResult.averageChi2ONdf) <=
          alignOptions.deltaAverageChi2ONdfCutOff.second) {
        ACTS_INFO(
            "Alignment has converged with change of chi2/ndf < "
            << alignOptions.deltaAverageChi2ONdfCutOff.second << " in the last "
            << alignOptions.deltaAverageChi2ONdfCutOff.first << " iterations"
            << " after " << iIter << " iteration(s)");
        converged = true;
        break;
      }
      recentChi2ONdf.pop();
    }
    // 2. or the average chi2/ndf
    if (alignResult.averageChi2ONdf <= alignOptions.averageChi2ONdfCutOff) {
      ACTS_INFO("Alignment has converged with average chi2/ndf < "
                << alignOptions.averageChi2ONdfCutOff << " after " << iIter
                << " iteration(s)");
      converged = true;
      break;
    }
    // Remove the first element
    // Store the result in the queue
    recentChi2ONdf.push(alignResult.averageChi2ONdf);

    // Not coveraged yet, update the detector element alignment parameters
    auto updateRes = updateAlignmentParameters(
        alignOptions.fitOptions.geoContext, alignOptions.alignedDetElements,
        alignOptions.alignmentTransformUpdater, alignResult);
    if (!updateRes.ok()) {
      ACTS_ERROR("Update alignment parameters failed: " << updateRes.error());
      return updateRes.error();
    }
  }
  // Write out the result
  for (const auto& det : alignOptions.alignedDetElements) {
    const auto& surface = &det->surface();
    const auto& transform = det->transform(alignOptions.fitOptions.geoContext);
    alignResult.alignedParameters.emplace(det, transform);
  }

  // Alignment failure if not converged
  if (!converged) {
    ACTS_ERROR("Alignment is not converged.");
  }

  return alignResult;
}
