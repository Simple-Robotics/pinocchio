//
// Copyright (c) 2019-2025 INRIA
//

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/math/matrix.hpp"
#include "pinocchio/math/eigen-helpers.hpp"

#include "pinocchio/utils/std-vector.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_matrix_matrix_product)
{
  using namespace pinocchio;
  using namespace Eigen;
  const Eigen::DenseIndex m = 20, n = 100;
  MatrixXd M1(MatrixXd::Ones(m, n)), M2(MatrixXd::Ones(n, m));
  MatrixMatrixProduct<MatrixXd, MatrixXd>::type res = M1 * M2;
  BOOST_CHECK(!res.eval().isZero());
}

BOOST_AUTO_TEST_CASE(test_scalar_matrix_product)
{
  using namespace pinocchio;
  using namespace Eigen;
  const Eigen::DenseIndex m = 20, n = 100;
  MatrixXd M(MatrixXd::Ones(m, n));
  const double alpha = 0.;
  ScalarMatrixProduct<double, MatrixXd>::type res = alpha * M;
  BOOST_CHECK(res.eval().isZero());
}

BOOST_AUTO_TEST_CASE(test_matrix_scalar_product)
{
  using namespace pinocchio;
  using namespace Eigen;
  const Eigen::DenseIndex m = 20, n = 100;
  MatrixXd M(MatrixXd::Ones(m, n));
  const double alpha = 1.;
  MatrixScalarProduct<MatrixXd, double>::type res = M * alpha;
  BOOST_CHECK(res.eval() == M);
}

BOOST_AUTO_TEST_CASE(test_eigen_helpers)
{
  using namespace pinocchio;
  using namespace Eigen;
  const Eigen::DenseIndex m = 20, n = 100;

  MatrixXd M(MatrixXd::Ones(m, n));

  setZero(M);
  BOOST_CHECK(M.isZero(0));
  setOnes(M);
  BOOST_CHECK(M.isOnes(0));
  setIdentity(M);
  BOOST_CHECK(M.isIdentity(0));
}

BOOST_AUTO_TEST_CASE(test_eigen_helpers_on_std_vector)
{
  using namespace pinocchio;
  using namespace Eigen;
  const Eigen::DenseIndex m = 20, n = 100;

  std::vector<Eigen::MatrixXd> vec(10, MatrixXd::Ones(m, n));

  apply_for_each(vec, setZero<Eigen::MatrixXd>);
  for (const auto & val : vec)
  {
    BOOST_CHECK(val.isZero(0));
  }

  apply_for_each(vec, setOnes<Eigen::MatrixXd>);
  for (const auto & val : vec)
  {
    BOOST_CHECK(val.isOnes(0));
  }

  apply_for_each(vec, setIdentity<Eigen::MatrixXd>);
  for (const auto & val : vec)
  {
    BOOST_CHECK(val.isIdentity(0));
  }
}

BOOST_AUTO_TEST_SUITE_END()
