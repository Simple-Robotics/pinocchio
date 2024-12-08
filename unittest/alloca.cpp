//
// Copyright (c) 2021-2024 CNRS INRIA
//

#include "pinocchio/fwd.hpp"

#ifdef WIN32
#else
  #include <alloca.h>
#endif

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

template<typename MatrixLikeInput>
typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLikeInput)
  copy(const Eigen::MatrixBase<MatrixLikeInput> & mat)
{
  typedef typename MatrixLikeInput::Scalar Scalar;
  void * alloca_ptr = alloca(size_t(mat.size()) * sizeof(Scalar));

  typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLikeInput) MatrixPlain;
  typedef Eigen::Map<MatrixPlain> MatrixPlainMap;
  MatrixPlainMap alloca_map =
    MatrixPlainMap(static_cast<Scalar *>(alloca_ptr), mat.rows(), mat.cols());

  alloca_map = mat; // copy
  return MatrixPlain(alloca_map);
}

BOOST_AUTO_TEST_CASE(copy_eigen_input)
{
  const int num_tests = 1000;
  const Eigen::DenseIndex rows = 10, cols = 20;
  const Eigen::MatrixXd mat = Eigen::MatrixXd::Random(rows, cols);

  for (int i = 0; i < num_tests; ++i)
  {
    const auto mat_copy = copy(mat);
    BOOST_CHECK(mat_copy == mat);
  }
}

BOOST_AUTO_TEST_SUITE_END()
