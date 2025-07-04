//
// Copyright (c) 2024 INRIA
//

#include <iostream>
#include "pinocchio/algorithm/constraints/box-set.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_proj)
{
  const int num_tests = int(1e6);
  const int dim = 10;

  const BoxSet box_constraint(-BoxSet::Vector::Ones(dim), BoxSet::Vector::Ones(dim));

  BOOST_CHECK(box_constraint.dim() == dim);

  BOOST_CHECK(box_constraint.isInside(BoxSet::Vector::Zero(dim)));
  BOOST_CHECK(box_constraint.project(BoxSet::Vector::Zero(dim)) == BoxSet::Vector::Zero(dim));

  for (int k = 0; k < num_tests; ++k)
  {
    const BoxSet::Vector x = BoxSet::Vector::Random(dim);

    // Cone
    const auto proj_x = box_constraint.project(x);
    const auto proj_proj_x = box_constraint.project(proj_x);

    BOOST_CHECK(box_constraint.isInside(proj_x, 1e-12));
    BOOST_CHECK(box_constraint.isInside(proj_proj_x, 1e-12));
    BOOST_CHECK(proj_x == proj_proj_x);
    if (box_constraint.isInside(x))
      BOOST_CHECK(x == proj_x);

    BOOST_CHECK(fabs((x - proj_x).dot(proj_x)) <= 1e-12); // orthogonal projection
  }
}

BOOST_AUTO_TEST_CASE(test_scaled_proj)
{
  const int num_tests = int(1e6);
  const int dim = 10;

  BoxSet box_constraint(-BoxSet::Vector::Ones(dim), BoxSet::Vector::Ones(dim));
  BoxSet::Vector scale = BoxSet::Vector::Random(dim);
  scale.array() =
    scale.array().max(BoxSet::Vector::Ones(dim).array() * 0.1); // ensure scale is positive
  const BoxSet scaled_box_constraint(
    -BoxSet::Vector::Ones(dim).cwiseQuotient(scale),
    BoxSet::Vector::Ones(dim).cwiseQuotient(scale));

  BOOST_CHECK(box_constraint.dim() == dim);

  BOOST_CHECK(box_constraint.isInside(BoxSet::Vector::Zero(dim)));
  BoxSet::Vector res(dim);
  box_constraint.scaledProject(BoxSet::Vector::Zero(dim), scale, res);
  BOOST_CHECK(res == BoxSet::Vector::Zero(dim));

  for (int k = 0; k < num_tests; ++k)
  {
    const BoxSet::Vector x = BoxSet::Vector::Random(dim);

    // Cone
    box_constraint.scaledProject(x, scale, res);
    const auto proj_x = res;
    box_constraint.scaledProject(proj_x, scale, res);
    const auto proj_proj_x = res;

    BOOST_CHECK(scaled_box_constraint.isInside(proj_x, 1e-12));
    BOOST_CHECK(scaled_box_constraint.isInside(proj_proj_x, 1e-12));
    BOOST_CHECK(proj_x == proj_proj_x);
    if (scaled_box_constraint.isInside(x))
      BOOST_CHECK(x == proj_x);

    BOOST_CHECK(fabs((x - proj_x).dot(proj_x)) <= 1e-12); // orthogonal projection
  }
}

BOOST_AUTO_TEST_SUITE_END()
