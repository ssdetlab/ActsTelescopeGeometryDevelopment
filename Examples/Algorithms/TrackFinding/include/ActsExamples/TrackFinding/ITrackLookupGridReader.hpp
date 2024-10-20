// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"

#include <memory>

namespace ActsExamples {

namespace Experimental {

using AxisType =
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Open>;
using TrackLookupGrid =
    Acts::Grid<std::unique_ptr<Acts::CurvilinearTrackParameters>, AxisType, AxisType>;

class ITrackLookupGridReader {
 public:
  /// Virtual Destructor
  virtual ~ITrackLookupGridReader() = default;

  /// Reader method
  virtual TrackLookupGrid readLookup(const std::string& path) = 0;
};

}  // namespace Experimental

}  // namespace ActsExamples
