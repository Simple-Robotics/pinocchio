#define PINOCCHIO_EIGEN_CHECK_MALLOC
#include <iostream>

#include <pinocchio/fwd.hpp>

#include <boost/variant.hpp> // to avoid C99 warnings

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

#include <pinocchio/algorithm/delassus-operator-dense.hpp>
#include <pinocchio/algorithm/delassus-operator-preconditioned.hpp>
#include "pinocchio/algorithm/preconditioner-diagonal.hpp"
#include <pinocchio/math/matrix.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

using namespace pinocchio;

BOOST_AUTO_TEST_CASE(delassus_dense_preconditioned)
{
  const Eigen::DenseIndex mat_size = 50;
  const Eigen::MatrixXd mat = Eigen::MatrixXd::Random(mat_size, mat_size);
  const Eigen::MatrixXd symmetric_mat = mat.transpose() * mat;
  const Eigen::VectorXd diag_vec = 1e-1 * Eigen::VectorXd::Ones(mat_size);
  const Eigen::MatrixXd preconditioner_matrix = diag_vec.asDiagonal();
  const Eigen::MatrixXd preconditioned_matrix =
    preconditioner_matrix * symmetric_mat * preconditioner_matrix;

  BOOST_CHECK(isSymmetric(symmetric_mat));
  BOOST_CHECK(isSymmetric(preconditioned_matrix));

  DelassusOperatorDense delassus(symmetric_mat);
  PreconditionerDiagonal<Eigen::VectorXd> diag_preconditioner(diag_vec);
  DelassusOperatorPreconditionedTpl<DelassusOperatorDense, PreconditionerDiagonal<Eigen::VectorXd>>
    delassus_preconditioned(delassus, diag_preconditioner);

  Eigen::VectorXd res(mat_size);
  const Eigen::VectorXd rhs = Eigen::VectorXd::Random(mat_size);

  // Checking matrix() method
  BOOST_CHECK(preconditioned_matrix.isApprox(delassus_preconditioned.matrix()));

  // Checking apply on the right
  delassus_preconditioned.applyOnTheRight(rhs, res);
  BOOST_CHECK(res.isApprox((preconditioned_matrix * rhs).eval()));

  // Checking solved
  Eigen::VectorXd damping_vec = 5e-3 * Eigen::VectorXd::Ones(mat_size);
  Eigen::MatrixXd damping = damping_vec.asDiagonal();
  delassus_preconditioned.updateDamping(damping_vec);
  delassus_preconditioned.solve(rhs, res);
  const Eigen::MatrixXd preconditioned_matrix_inv = (preconditioned_matrix + damping).inverse();
  Eigen::VectorXd res_solve = preconditioned_matrix_inv * rhs;
  BOOST_CHECK(res.isApprox(res_solve));

  // Checking solveInPlace
  delassus_preconditioned.solveInPlace(rhs);
  BOOST_CHECK(rhs.isApprox(res_solve));

  // Checking updateDamping
  const double new_damping = 1e-3;
  const Eigen::MatrixXd damped_preconditioned_matrix =
    preconditioned_matrix + new_damping * Eigen::MatrixXd::Identity(mat_size, mat_size);
  delassus_preconditioned.updateDamping(new_damping);
  delassus_preconditioned.applyOnTheRight(rhs, res);
  Eigen::VectorXd res_apply = damped_preconditioned_matrix * rhs;
  BOOST_CHECK(res.isApprox(res_apply));
}

BOOST_AUTO_TEST_SUITE_END()
