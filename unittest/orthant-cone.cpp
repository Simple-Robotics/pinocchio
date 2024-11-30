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
void common_test(const int num_tests, const int dim)
{

  const Orthant orthant(dim);
  typedef typename Orthant::Vector Vector;

  BOOST_CHECK(orthant.dim() == dim);

  BOOST_CHECK(orthant.isInside(Vector::Zero(dim)));
  BOOST_CHECK(orthant.project(Vector::Zero(dim)) == Vector::Zero(dim));
  BOOST_CHECK(orthant.dual() == orthant);

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

  common_test<PositiveOrthantCone>(num_tests, dim);

  const PositiveOrthantCone positive_orthant(dim);
  typedef PositiveOrthantCone::Vector Vector;

  for (int k = 0; k < num_tests; ++k)
  {
    const Vector x = Vector::Random(dim);
    const Vector x_proj = positive_orthant.project(x);

    BOOST_CHECK((x_proj.array() >= 0).all());
  }
}

BOOST_AUTO_TEST_CASE(test_negative_orthant)
{
  const int num_tests = int(1e6);
  const int dim = 10;

  common_test<NegativeOrthantCone>(num_tests, dim);

  const NegativeOrthantCone negative_orthant(dim);
  typedef NegativeOrthantCone::Vector Vector;

  for (int k = 0; k < num_tests; ++k)
  {
    const Vector x = Vector::Random(dim);
    const Vector x_proj = negative_orthant.project(x);

    BOOST_CHECK((x_proj.array() <= 0).all());
  }
}

BOOST_AUTO_TEST_CASE(test_decomposition)
{
  const int num_tests = int(1e6);
  const int dim = 10;

  const NegativeOrthantCone negative_orthant(dim);
  const PositiveOrthantCone positive_orthant(dim);
  typedef NegativeOrthantCone::Vector Vector;

  for (int k = 0; k < num_tests; ++k)
  {
    const Vector x = Vector::Random(dim);
    const Vector x_proj_neg = negative_orthant.project(x);
    const Vector x_proj_pos = positive_orthant.project(x);

    BOOST_CHECK(x == (x_proj_neg + x_proj_pos));
  }
}

BOOST_AUTO_TEST_SUITE_END()
