//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/container/storage.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
typedef EigenStorageTpl<double> EigenStorage;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(eigen_storage)
{
  const Eigen::DenseIndex rows = 10, cols = 20;

  const Eigen::DenseIndex initial_capacity = rows * cols;
  EigenStorage storage(initial_capacity);

  BOOST_CHECK(storage.capacity() == initial_capacity);

  auto matrix_map = storage.as<Eigen::MatrixXd>(rows, cols);
  BOOST_CHECK(matrix_map.data() == storage.data());

  matrix_map.setIdentity();
  BOOST_CHECK(storage.as<Eigen::MatrixXd>(rows, cols).isIdentity(0.));
  BOOST_CHECK(
    static_cast<const EigenStorage &>(storage).as<Eigen::MatrixXd>(rows, cols).isIdentity(0.));
  matrix_map.setOnes();
  BOOST_CHECK(storage.as<Eigen::VectorXd>(initial_capacity).isOnes(0.));
  BOOST_CHECK(
    static_cast<const EigenStorage &>(storage).as<Eigen::VectorXd>(initial_capacity).isOnes(0.));

  // Check resize
  storage.conservativeResize(2 * storage.capacity());
  BOOST_CHECK(storage.as<Eigen::MatrixXd>(rows, cols).isOnes(0.));
  matrix_map = storage.as<Eigen::MatrixXd>(rows, 2 * cols);
  BOOST_CHECK(matrix_map.leftCols(cols).isOnes(0.));

  // Check throw
  BOOST_CHECK_THROW(storage.as<Eigen::MatrixXd>(2 * rows, 2 * cols), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
