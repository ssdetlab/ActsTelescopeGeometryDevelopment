// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/TrackParametersJsonConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Acts {

void to_json(nlohmann::json& j, const CurvilinearTrackParameters& t) {
    j["pos"] = t.fourPosition(); 
    j["dir"] = t.direction();
}

void from_json(nlohmann::json& j, const CurvilinearTrackParameters& t) {
}

}  // namespace Acts
