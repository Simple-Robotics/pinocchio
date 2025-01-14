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
  Eigen::VectorXd x_scaled_true = Eigen::VectorXd::Ones(n);
  x_scaled_true.array() /= precond_vec.array();
  precond.scale(x, x_scaled);
  BOOST_CHECK((x_scaled - x_scaled_true).isZero());
  Eigen::VectorXd x_unscaled = x;
  precond.unscale(x, x_unscaled);
  BOOST_CHECK((x_unscaled - precond_vec).isZero());
}

BOOST_AUTO_TEST_SUITE_END()
