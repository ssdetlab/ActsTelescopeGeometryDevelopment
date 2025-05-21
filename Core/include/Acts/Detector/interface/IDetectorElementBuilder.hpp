// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include <memory>
#include <vector>

namespace Acts::Experimental {

/// @brief  This is the interface for detector component builders;
/// such a builder could be a simple detector volume builder, with
/// or without internal structure, or more complicated objects.
///
class IDetectorElementBuilder {
 public:
  virtual ~IDetectorElementBuilder() = default;

  /// The interface method to be implemented by all detector
  /// component builder
  ///
  /// @param gctx The geometry context for this call
  ///
  /// @return an outgoing detector component
  virtual std::vector<std::shared_ptr<DetectorElementBase>> construct(
      std::vector<std::shared_ptr<Acts::Surface>> surfaces,
      const GeometryContext& gctx) const = 0;
};

}  // namespace Acts::Experimental
