//
// Copyright (c) 2024 INRIA
//

#define PINOCCHIO_EIGEN_CHECK_MALLOC
#include <iostream>

#include <pinocchio/fwd.hpp>

#include <boost/variant.hpp> // to avoid C99 warnings

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

#include <pinocchio/algorithm/delassus-operator-dense.hpp>
#include <pinocchio/math/matrix.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

using namespace pinocchio;

BOOST_AUTO_TEST_CASE(test_memory_allocation)
{
  const Eigen::DenseIndex mat_size = 100;
  const Eigen::MatrixXd mat_ = Eigen::MatrixXd::Random(mat_size,mat_size);
  const Eigen::MatrixXd symmetric_mat = mat_.transpose() * mat_;
  
  BOOST_CHECK(isSymmetric(symmetric_mat));
  
  DelassusOperatorDense delassus(symmetric_mat);
  
  Eigen::VectorXd res(mat_size);
  const Eigen::VectorXd rhs = Eigen::VectorXd::Random(mat_size);
  res = delassus * rhs;
  BOOST_CHECK(res.isApprox((symmetric_mat * rhs).eval()));
  
  PowerIterationAlgoTpl<Eigen::VectorXd> power_iteration(mat_size);
  
  // Check memory allocations
  Eigen::internal::set_is_malloc_allowed(false);
  res = delassus * rhs;
  (delassus * rhs).evalTo(res);
  res.noalias() = symmetric_mat * rhs;
  res.noalias() = delassus * rhs;
  evalTo(symmetric_mat * rhs, res);
  power_iteration.run(delassus);
  power_iteration.run(symmetric_mat);
  Eigen::internal::set_is_malloc_allowed(true);
}

BOOST_AUTO_TEST_SUITE_END()

