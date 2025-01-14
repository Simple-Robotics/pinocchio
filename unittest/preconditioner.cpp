//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/algorithm/preconditioner-diagonal.hpp"

#include <boost/test/unit_test.hpp>

using namespace pinocchio;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(diagonal_preconditioner)
{
  std::size_t n = 10;
  Eigen::VectorXd precond_vec = Eigen::VectorXd::Ones(n);
  precond_vec.head(n / 2).setConstant(0.1);
  precond_vec.tail(n / 2).setConstant(12.4);
  PreconditionerDiagonal<Eigen::VectorXd> precond(precond_vec);
  Eigen::VectorXd x = Eigen::VectorXd::Ones(n);
  Eigen::VectorXd x_scaled = x;
  precond.scale(x, x_scaled);
  BOOST_CHECK((x_scaled - precond_vec).isZero());
  Eigen::VectorXd x_unscaled = x;
  precond.unscale(x, x_unscaled);
  x.array() /= precond_vec.array();
  BOOST_CHECK((x_unscaled - x).isZero());
}

BOOST_AUTO_TEST_SUITE_END()
