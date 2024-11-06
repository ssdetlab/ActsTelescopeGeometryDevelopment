// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "ActsExamples/TrackFinding/TrackParamsLookupTable.hpp"

#include <memory>

namespace ActsExamples {

class ITrackParamsLookupWriter {
    public:
        /// Virtual Destructor
        virtual ~ITrackParamsLookupWriter() = default;
        
        /// Writer method
        ///
        /// @param lookup the estimated track grid to write
        virtual void writeLookup(
            const Lookup& lookup) = 0;
};

}  // namespace ActsExamples
