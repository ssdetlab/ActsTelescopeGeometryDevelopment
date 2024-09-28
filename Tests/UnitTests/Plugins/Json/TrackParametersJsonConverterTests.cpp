// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/Json/TrackParametersJsonConverter.hpp"

#include <fstream>
#include <memory>
#include <string>

#include <nlohmann/json.hpp>

using namespace Acts;

namespace {
std::ofstream out;

Acts::GeometryContext gctx;
}  // namespace

BOOST_AUTO_TEST_SUITE(TrackParametersJsonConversion)

BOOST_AUTO_TEST_CASE(CurvilinearTrackParametersTests) {
    Acts::Vector4 pos4(1., 2., 3., 4.);
    Acts::Vector3 dir(0.1, 0.2, 0.3);
    Acts::ActsScalar qOverP = 0.1;
    Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Identity();
    Acts::ParticleHypothesis particle = Acts::ParticleHypothesis::electron();

    Acts::CurvilinearTrackParameters pars(pos4, dir, qOverP, cov, particle);

    nlohmann::json jParsOut;
    Acts::to_json(jParsOut, pars);

    out.open("CurvilinearTrackParameters.json");
    out << jParsOut.dump(4);
    out.close();

    // std::ifstream in("CurvilinearTrackParameters.json",
                    //  std::ifstream::in | std::ifstream::binary);
    // BOOST_CHECK(in.good());

    // nlohmann::json jParsIn;
    // in >> jParsIn;
    // in.close();

    // Acts::CurvilinearTrackParameters* parsIn = nullptr;

    // Acts::from_json(jParsIn, *parsIn);

    // BOOST_CHECK(parsIn->position() == pars.position());
    // BOOST_CHECK(parsIn->direction() == pars.direction());
    // BOOST_CHECK(parsIn->qOverP() == pars.qOverP());
    // BOOST_CHECK(parsIn->covariance() == pars.covariance());
    // BOOST_CHECK(parsIn->particleHypothesis() == pars.particleHypothesis());
}

BOOST_AUTO_TEST_SUITE_END()
