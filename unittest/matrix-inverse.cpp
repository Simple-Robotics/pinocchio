//
// Copyright (c) 2025 INRIA
//

#include <iostream>

#include <pinocchio/math/matrix-inverse.hpp>

#include <boost/variant.hpp> // to avoid C99 warnings

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

using namespace pinocchio;

template<int size>
using MatrixTpl = Eigen::Matrix<double, size, size>;

template<int size>
void test_generated_inverse_impl()
{
  typedef MatrixTpl<4> Matrix;
  Matrix mat = Matrix::Random();
  mat = mat.transpose() * mat;
  make_symmetric(mat);
  BOOST_CHECK(is_symmetric(mat, 0));

  Matrix res = Matrix::Zero();
  matrix_inversion_code_generated(mat, res);
  BOOST_CHECK((res * mat).isIdentity());
  BOOST_CHECK(mat.inverse().isApprox(res));
}

BOOST_AUTO_TEST_CASE(test_generated_inverse)
{
  test_generated_inverse_impl<3>();
  test_generated_inverse_impl<4>();
  test_generated_inverse_impl<5>();
  test_generated_inverse_impl<6>();
  test_generated_inverse_impl<7>();
}

BOOST_AUTO_TEST_SUITE_END()
