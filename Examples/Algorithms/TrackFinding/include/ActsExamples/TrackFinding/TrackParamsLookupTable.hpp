// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"

namespace ActsExamples {

using LookupAxis =
    Acts::Axis<
        Acts::AxisType::Equidistant, 
        Acts::AxisBoundaryType::Open>;

using LookupAxisGen = 
    Acts::GridAxisGenerators::EqOpenEqOpen;

using LookupAccumGrid =
    Acts::Grid<
        std::vector<Acts::CurvilinearTrackParameters>, 
        LookupAxis, 
        LookupAxis>;

using LookupAccumGridContainer = 
    std::unordered_map<
        Acts::GeometryIdentifier, LookupAccumGrid>;

using LookupPair = std::pair<
    std::shared_ptr<Acts::CurvilinearTrackParameters>,
    std::shared_ptr<Acts::CurvilinearTrackParameters>>;

using LookupGrid = 
    Acts::Grid<LookupPair, LookupAxis, LookupAxis>;

using Lookup = 
    std::unordered_map<
        Acts::GeometryIdentifier, LookupGrid>;

}  // namespace ActsExamples
