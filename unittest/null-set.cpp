//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/algorithm/constraints/null-set.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_proj)
{
  const int num_tests = int(1000);
  const int dim = 10;

  const NullSet set(dim);

  BOOST_CHECK(set.dim() == dim);

  BOOST_CHECK(set.isInside(NullSet::Vector::Zero(dim)));
  BOOST_CHECK(!set.isInside(NullSet::Vector::Ones(dim)));
  BOOST_CHECK(set.project(NullSet::Vector::Zero(dim)) == NullSet::Vector::Zero(dim));
  BOOST_CHECK(set.project(NullSet::Vector::Ones(dim)) == NullSet::Vector::Zero(dim));

  for (int k = 0; k < num_tests; ++k)
  {
    const NullSet::Vector x = NullSet::Vector::Random(dim);
    if (!x.isZero())
      BOOST_CHECK(!set.isInside(x));

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
