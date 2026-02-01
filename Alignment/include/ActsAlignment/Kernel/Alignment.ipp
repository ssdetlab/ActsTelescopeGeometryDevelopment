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
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"

#include <Eigen/src/Core/Matrix.h>

template <typename fitter_t>
template <typename source_link_t, typename start_parameters_t,
          typename fit_options_t>
Acts::Result<ActsAlignment::detail::TrackAlignmentState>
ActsAlignment::Alignment<fitter_t>::evaluateTrackAlignmentState(
    const Acts::GeometryContext& gctx,
    const std::vector<source_link_t>& sourcelinks,
    const start_parameters_t& sParameters, const fit_options_t& fitOptions,
    const std::unordered_map<const Acts::Surface*, std::size_t>&
        idxedAlignSurfaces,
    const ActsAlignment::AlignmentMask& alignMask) const {
  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  // Convert to Acts::SourceLink during iteration
  Acts::SourceLinkAdapterIterator begin{sourcelinks.begin()};
  Acts::SourceLinkAdapterIterator end{sourcelinks.end()};

  // Perform the fit
  auto fitRes = m_fitter.fit(begin, end, sParameters, fitOptions, tracks);

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
    const fit_options_t& fitOptions,
    ActsAlignment::AlignmentResult& alignResult, const AlignmentMask& alignMask,
    const ActsAlignment::AlignmentMode& alignMode) const {
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
  // @Todo: How to update the source link error iteratively?
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

  // Get the inverse of chi2 second derivative matrix (we need this to
  // calculate the covariance of the alignment parameters)
  // @Todo: use more stable method for solving the inverse
  std::size_t alignDof = alignResult.alignmentDof;
  Acts::ActsDynamicMatrix sumChi2SecondDerivativeInverse =
      Acts::ActsDynamicMatrix::Zero(alignDof, alignDof);
  sumChi2SecondDerivativeInverse = sumChi2SecondDerivative.inverse();
  if (sumChi2SecondDerivativeInverse.hasNaN()) {
    ACTS_DEBUG("Chi2 second derivative inverse has NaN");
    // return AlignmentError::AlignmentParametersUpdateFailure;
  }

  //-----------------------------------------------------------------------------
  alignResult.deltaAlignmentParameters =
      Acts::ActsDynamicVector::Zero(alignDof);
  alignResult.alignmentCovariance =
      Acts::ActsDynamicMatrix::Zero(alignDof, alignDof);
  if (alignMode == ActsAlignment::AlignmentMode::local) {
    const int nConstraints = 3;

    // Constraint handling
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(nConstraints, alignDof);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(nConstraints);

    for (const auto& [surf, idx] : alignResult.idxedAlignSurfaces) {
      C(0, idx * Acts::eAlignmentSize + Acts::eAlignmentCenter1) = 1.0;
      C(1, idx * Acts::eAlignmentSize + Acts::eAlignmentCenter2) = 1.0;
      C(2, idx * Acts::eAlignmentSize + Acts::eAlignmentRotation2) = 1.0;
    }

    Eigen::MatrixXd A =
        Eigen::MatrixXd::Zero(sumChi2SecondDerivative.rows() + C.rows(),
                              sumChi2SecondDerivative.cols() + C.rows());
    A.topLeftCorner(sumChi2SecondDerivative.rows(),
                    sumChi2SecondDerivative.cols()) = sumChi2SecondDerivative;
    A.topRightCorner(sumChi2SecondDerivative.rows(), C.rows()) = C.transpose();
    A.bottomLeftCorner(C.rows(), sumChi2SecondDerivative.cols()) = C;

    Eigen::VectorXd rhs(sumChi2SecondDerivative.rows() + C.rows());
    rhs.head(sumChi2SecondDerivative.rows()) = -sumChi2Derivative;
    rhs.tail(C.rows()) = b;

    // Solve
    Eigen::VectorXd sol = A.fullPivLu().solve(rhs);

    alignResult.deltaAlignmentParameters = sol.head(alignDof);

    // Alignment parameters covariance
    alignResult.alignmentCovariance = 2 * sumChi2SecondDerivativeInverse;
  } else if (alignMode == ActsAlignment::AlignmentMode::global) {
    const int nRigid = 6;

    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(alignDof, nRigid);

    bool doCenter0 = ACTS_CHECK_BIT(alignMask, AlignmentMask::Center0);
    bool doCenter1 = ACTS_CHECK_BIT(alignMask, AlignmentMask::Center1);
    bool doCenter2 = ACTS_CHECK_BIT(alignMask, AlignmentMask::Center2);
    bool doAngle0 = ACTS_CHECK_BIT(alignMask, AlignmentMask::Rotation0);
    bool doAngle1 = ACTS_CHECK_BIT(alignMask, AlignmentMask::Rotation1);
    bool doAngle2 = ACTS_CHECK_BIT(alignMask, AlignmentMask::Rotation2);

    Acts::Vector3 centerOfMass = Acts::Vector3::Zero();
    for (const auto& [surf, idx] : alignResult.idxedAlignSurfaces) {
      centerOfMass += surf->center(gctx);
    }
    centerOfMass /= alignResult.idxedAlignSurfaces.size();

    for (const auto& [surf, idx] : alignResult.idxedAlignSurfaces) {
      if (doAngle0) {
        U.block(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter0, 3, 3, 1) =
            Acts::Vector3::UnitX().cross(surf->center(gctx) - centerOfMass);
        U(idx * Acts::eAlignmentSize + Acts::eAlignmentRotation0, 3) = 1;
      }
      if (doAngle1) {
        U.block(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter0, 4, 3, 1) =
            Acts::Vector3::UnitY().cross(surf->center(gctx) - centerOfMass);
        U(idx * Acts::eAlignmentSize + Acts::eAlignmentRotation1, 4) = 1;
      }
      if (doAngle2) {
        U.block(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter0, 5, 3, 1) =
            Acts::Vector3::UnitZ().cross(surf->center(gctx) - centerOfMass);
        U(idx * Acts::eAlignmentSize + Acts::eAlignmentRotation2, 5) = 1;
      }

      if (doCenter0) {
        U(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter0, 0) = 1;
      }
      if (doCenter1) {
        U(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter1, 1) = 1;
      }
      if (doCenter2) {
        U(idx * Acts::eAlignmentSize + Acts::eAlignmentCenter2, 2) = 1;
      }
    }

    std::cout << "RIGID PROJECTOR\n" << U << "\n";
    Eigen::MatrixXd rigidA = U.transpose() * sumChi2SecondDerivative * U;
    Eigen::VectorXd rigidb = -U.transpose() * sumChi2Derivative;

    std::cout << "RIGID DERIVATIVE\n" << -rigidb << "\n";

    double eps = 1e-12 * std::max(1.0, rigidA.norm());
    if (rigidA.fullPivLu().rank() < nRigid) {
      rigidA.diagonal().array() += eps;
    }

    // Solve for z
    Eigen::VectorXd rigidDelta = rigidA.fullPivLu().solve(rigidb);

    alignResult.deltaAlignmentParameters = U * rigidDelta;

    // Alignment parameters covariance
    alignResult.alignmentCovariance = 2 * rigidA.inverse();
  }
  ACTS_VERBOSE("sumChi2SecondDerivative = \n" << sumChi2SecondDerivative);
  ACTS_VERBOSE("sumChi2Derivative = \n" << sumChi2Derivative);
  std::cout << "sumChi2Derivative = \n" << sumChi2Derivative << "\n";
  ACTS_VERBOSE("alignResult.deltaAlignmentParameters \n");

  // chi2 change
  alignResult.deltaChi2 = 0.5 * sumChi2Derivative.transpose() *
                          alignResult.deltaAlignmentParameters;
}

template <typename fitter_t>
Acts::Result<void>
ActsAlignment::Alignment<fitter_t>::updateAlignmentParameters(
    const Acts::GeometryContext& gctx,
    const std::vector<Acts::DetectorElementBase*>& alignedDetElements,
    const ActsAlignment::AlignedTransformUpdater& alignedTransformUpdater,
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
    //@Todo: use a better way to handle this (need dynamic cast to inherited
    // detector element type)
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
  bool alignmentParametersUpdated = false;
  std::queue<double> recentChi2ONdf;

  // Perform the fit to the trajectories and update alignment parameters
  ACTS_INFO("Max number of iterations: " << alignOptions.maxIterations);
  for (unsigned int iIter = 0; iIter < alignOptions.maxIterations; iIter++) {
    // Calculate the alignment parameters delta etc.
    calculateAlignmentParameters(
        alignOptions.fitOptions.geoContext, trajectoryCollection,
        startParametersCollection, alignOptions.fitOptions, alignResult,
        alignOptions.alignmentMask, alignOptions.alignmentMode);
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
    // 2. or the average chi2/ndf (is this correct?)
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

    ACTS_INFO("The solved delta of alignmentParameters = \n "
              << alignResult.deltaAlignmentParameters);
    // Not coveraged yet, update the detector element alignment parameters
    auto updateRes = updateAlignmentParameters(
        alignOptions.fitOptions.geoContext, alignOptions.alignedDetElements,
        alignOptions.alignedTransformUpdater, alignResult);
    if (!updateRes.ok()) {
      ACTS_ERROR("Update alignment parameters failed: " << updateRes.error());
      return updateRes.error();
    }
    alignmentParametersUpdated = true;
  }  // end of all iterations

  // Alignment failure if not converged
  if (!converged) {
    ACTS_ERROR("Alignment is not converged.");
    alignResult.result = AlignmentError::ConvergeFailure;
  }

  // Screen out the final aligned parameters
  // @todo
  if (alignmentParametersUpdated) {
    for (const auto& det : alignOptions.alignedDetElements) {
      const auto& surface = &det->surface();
      const auto& transform =
          det->transform(alignOptions.fitOptions.geoContext);
      // write it to the result
      alignResult.alignedParameters.emplace(det, transform);
      const auto& translation = transform.translation();
      const auto& rotation = transform.rotation();
      const Acts::Vector3 rotAngles = rotation.eulerAngles(2, 1, 0);
      ACTS_VERBOSE("Detector element with surface "
                   << surface->geometryId()
                   << " has aligned geometry position as below:");
      ACTS_VERBOSE("Center (cenX, cenY, cenZ) = " << translation.transpose());
      ACTS_VERBOSE(
          "Euler angles (rotZ, rotY, rotX) = " << rotAngles.transpose());
      ACTS_VERBOSE("Rotation matrix = \n" << rotation);
    }
  } else {
    ACTS_DEBUG("Alignment parameters is not updated.");
  }

  return alignResult;
}
