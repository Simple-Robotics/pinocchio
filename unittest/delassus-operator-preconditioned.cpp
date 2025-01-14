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

BOOST_AUTO_TEST_CASE(delassus_preconditioned_applyontheright)
{
  const Eigen::DenseIndex mat_size = 50;
  const Eigen::MatrixXd mat_ = Eigen::MatrixXd::Random(mat_size, mat_size);
  const Eigen::MatrixXd symmetric_mat = mat_.transpose() * mat_;
  const Eigen::VectorXd diag_vec = Eigen::VectorXd::Random(mat_size);
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
  delassus_preconditioned.applyOnTheRight(rhs, res);
  BOOST_CHECK(res.isApprox((preconditioned_matrix * rhs).eval()));
}

BOOST_AUTO_TEST_SUITE_END()
