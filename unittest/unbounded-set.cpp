//
// Copyright (c) 2024 INRIA
//

#include <iostream>
#include "pinocchio/algorithm/constraints/unbounded-set.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_proj)
{
  const int num_tests = int(1e6);
  const int dim = 10;

  const UnboundedSet set(dim);

  BOOST_CHECK(set.dim() == dim);

  BOOST_CHECK(set.isInside(UnboundedSet::Vector::Zero(dim)));
  BOOST_CHECK(set.project(UnboundedSet::Vector::Zero(dim)) == UnboundedSet::Vector::Zero(dim));

  for (int k = 0; k < num_tests; ++k)
  {
    const UnboundedSet::Vector x = UnboundedSet::Vector::Random(dim);

    const auto proj_x = set.project(x);
    const auto proj_proj_x = set.project(proj_x);

    BOOST_CHECK(set.isInside(proj_x, 1e-12));
    BOOST_CHECK(set.isInside(proj_proj_x, 1e-12));
    BOOST_CHECK(proj_x == proj_proj_x);
    if (set.isInside(x))
      BOOST_CHECK(x == proj_x);

    BOOST_CHECK(fabs((x - proj_x).dot(proj_x)) <= 1e-12); // orthogonal projection
  }
}

BOOST_AUTO_TEST_SUITE_END()
