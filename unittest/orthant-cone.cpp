//
// Copyright (c) 2024 INRIA
//

#include <iostream>
#include "pinocchio/algorithm/constraints/orthant-cone.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

template<typename Orthant>
void main_test(const int num_tests, const int dim)
{

  const Orthant orthant(dim);
  typedef typename Orthant::Vector Vector;

  BOOST_CHECK(orthant.dim() == dim);

  BOOST_CHECK(orthant.isInside(Vector::Zero(dim)));
  BOOST_CHECK(orthant.project(Vector::Zero(dim)) == Vector::Zero(dim));
  BOOST_CHECK(&orthant.dual() == &orthant);

  for (int k = 0; k < num_tests; ++k)
  {
    const Vector x = Vector::Random(dim);

    // Cone
    const auto proj_x = orthant.project(x);
    const auto proj_proj_x = orthant.project(proj_x);

    BOOST_CHECK(orthant.isInside(proj_x, 1e-12));
    BOOST_CHECK(orthant.isInside(proj_proj_x, 1e-12));
    BOOST_CHECK(proj_x == proj_proj_x);
    if (orthant.isInside(x))
      BOOST_CHECK(x == proj_x);

    BOOST_CHECK(fabs((x - proj_x).dot(proj_x)) <= 1e-12); // orthogonal projection
  }
}

BOOST_AUTO_TEST_CASE(test_positive_orthant)
{
  const int num_tests = int(1e6);
  const int dim = 10;

  main_test<PositiveOrthantCone>(num_tests, dim);
}

BOOST_AUTO_TEST_CASE(test_negative_orthant)
{
  const int num_tests = int(1e6);
  const int dim = 10;

  main_test<NegativeOrthantCone>(num_tests, dim);
}

BOOST_AUTO_TEST_SUITE_END()
