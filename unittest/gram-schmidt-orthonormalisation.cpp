//
// Copyright (c) 2024-2025 INRIA
//

#include <iostream>

#include <pinocchio/math/gram-schmidt-orthonormalisation.hpp>
#include <Eigen/QR>

#include <boost/variant.hpp> // to avoid C99 warnings

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

using namespace pinocchio;

BOOST_AUTO_TEST_CASE(test_orthogonalization)
{
  for (size_t i = 0; i < 100; ++i)
  {
    const Eigen::DenseIndex size = 20;
    const Eigen::MatrixXd random_mat = Eigen::MatrixXd::Random(size, size);
    const auto qr = random_mat.householderQr();
    const Eigen::MatrixXd basis = qr.householderQ();

    for (size_t k = 0; k < 1000; ++k)
    {
      const Eigen::VectorXd random_vec = Eigen::VectorXd::Random(size);
      orthogonalization(basis.leftCols(10), random_vec);
      BOOST_CHECK((basis.leftCols(10).transpose() * random_vec).isZero());
    }
  }
}

BOOST_AUTO_TEST_CASE(test_orthonormalization)
{
  const size_t num_tests =
#ifdef NDEBUG
    10000
#else
    100
#endif
    ;
  for (size_t i = 0; i < num_tests; ++i)
  {
    const Eigen::DenseIndex size = 100;
    const double prec = size * size * Eigen::NumTraits<double>::dummy_precision();
    const Eigen::MatrixXd random_mat = Eigen::MatrixXd::Random(size, size);
    const Eigen::MatrixXd mat = random_mat + Eigen::MatrixXd::Identity(size, size);

    BOOST_CHECK(!isOrthonormal(mat, prec));
    orthonormalization(mat);
    BOOST_CHECK(isOrthonormal(mat, prec));
  }
}

BOOST_AUTO_TEST_SUITE_END()
