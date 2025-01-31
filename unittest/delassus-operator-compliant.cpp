#define PINOCCHIO_EIGEN_CHECK_MALLOC

#include <pinocchio/fwd.hpp>

#include <boost/variant.hpp> // to avoid C99 warnings

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

#include <pinocchio/algorithm/delassus-operator-dense.hpp>
#include <pinocchio/algorithm/delassus-operator-compliant.hpp>
#include <pinocchio/math/matrix.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

using namespace pinocchio;

BOOST_AUTO_TEST_CASE(delassus_dense_compliant)
{
  const Eigen::DenseIndex mat_size = 50;
  const Eigen::MatrixXd mat = Eigen::MatrixXd::Random(mat_size, mat_size);
  const Eigen::MatrixXd symmetric_mat = mat.transpose() * mat;
  const Eigen::VectorXd diag_vec = 1e-1 * Eigen::VectorXd::Ones(mat_size);
  const Eigen::MatrixXd compliance_matrix = diag_vec.asDiagonal();
  const Eigen::MatrixXd compliant_matrix = symmetric_mat + compliance_matrix;

  BOOST_CHECK(isSymmetric(symmetric_mat));
  BOOST_CHECK(isSymmetric(compliant_matrix));

  DelassusOperatorDense delassus(symmetric_mat);
  Eigen::VectorXd compliance(diag_vec);
  DelassusOperatorCompliantTpl<DelassusOperatorDense, Eigen::VectorXd> delassus_compliant(
    delassus, compliance);

  Eigen::VectorXd res(mat_size);
  const Eigen::VectorXd rhs = Eigen::VectorXd::Random(mat_size);

  // Checking matrix() method
  BOOST_CHECK(compliant_matrix.isApprox(delassus_compliant.matrix()));

  // Checking apply on the right
  delassus_compliant.applyOnTheRight(rhs, res);
  BOOST_CHECK(res.isApprox((compliant_matrix * rhs).eval()));

  // Checking solved
  Eigen::VectorXd damping_vec = 5e-3 * Eigen::VectorXd::Ones(mat_size);
  Eigen::MatrixXd damping = damping_vec.asDiagonal();
  delassus_compliant.updateDamping(damping_vec);
  delassus_compliant.solve(rhs, res);
  const Eigen::MatrixXd compliant_matrix_inv = (compliant_matrix + damping).inverse();
  Eigen::VectorXd res_solve = compliant_matrix_inv * rhs;
  BOOST_CHECK(res.isApprox(res_solve));

  // Checking solveInPlace
  delassus_compliant.solveInPlace(rhs);
  BOOST_CHECK(rhs.isApprox(res_solve));

  // Checking updateDamping
  const double new_damping = 1e-3;
  const Eigen::MatrixXd damped_compliant_matrix =
    compliant_matrix + new_damping * Eigen::MatrixXd::Identity(mat_size, mat_size);
  delassus_compliant.updateDamping(new_damping);
  delassus_compliant.applyOnTheRight(rhs, res);
  Eigen::VectorXd res_apply = damped_compliant_matrix * rhs;
  BOOST_CHECK(res.isApprox(res_apply));
}

BOOST_AUTO_TEST_SUITE_END()
