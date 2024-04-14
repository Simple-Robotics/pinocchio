//
// Copyright (c) 2024 INRIA
//

#include <pinocchio/math/matrix.hpp>

#include <boost/variant.hpp> // to avoid C99 warnings

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_isSymmetric)
{
  srand(0);

  using namespace pinocchio;
  typedef Eigen::MatrixXd Matrix;
  
#ifdef NDEBUG
  const int max_test = 1e6;
  const int max_size = 1000;
#else
  const int max_test = 1e2;
  const int max_size = 100;
#endif
  for(int i = 0; i < max_test; ++i)
  {
    const Eigen::DenseIndex rows = rand() % max_size + 3; // random row number
    const Eigen::DenseIndex cols = rand() % max_size + 3; // random col number
    
    const Matrix random_matrix = Matrix::Random(rows,cols);
    Matrix symmetric_matrix = random_matrix.transpose() * random_matrix;
    BOOST_CHECK(isSymmetric(symmetric_matrix));
    
    symmetric_matrix.coeffRef(1,0) += 2.; 
    BOOST_CHECK(!isSymmetric(symmetric_matrix));
    
    // Specific check for the Zero matrix
    BOOST_CHECK(isSymmetric(Matrix::Zero(rows,rows)));
    // Specific check for the Identity matrix
    BOOST_CHECK(isSymmetric(Matrix::Identity(rows,rows)));
  }
}

BOOST_AUTO_TEST_SUITE_END()
