//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/container/storage.hpp"

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
typedef EigenStorageTpl<Eigen::MatrixXd> EigenStorageMatrix;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(eigen_storage_matrix)
{
  const Eigen::DenseIndex rows = 10, cols = 20;

  const Eigen::DenseIndex initial_capacity = rows * cols;
  EigenStorageMatrix storage(rows, cols);

  BOOST_CHECK(storage.capacity() == initial_capacity);
  BOOST_CHECK(storage.rows() == rows);
  BOOST_CHECK(storage.cols() == cols);

  EigenStorageMatrix::RefMapType matrix_map = storage.map();
  BOOST_CHECK(matrix_map.data() == storage.data());

  matrix_map.setIdentity();
  BOOST_CHECK(storage.map().isIdentity(0.));
  BOOST_CHECK(static_cast<const EigenStorageMatrix &>(storage).map().isIdentity(0.));
  matrix_map.setOnes();
  BOOST_CHECK(storage.map().isOnes(0.));
  BOOST_CHECK(static_cast<const EigenStorageMatrix &>(storage).map().isOnes(0.));

  // Check resize
  const Eigen::DenseIndex new_rows = 2 * rows, new_cols = cols;
  storage.conservativeResize(new_rows, new_cols);
  BOOST_CHECK(matrix_map.data() == storage.data());
  BOOST_CHECK(storage.map().topLeftCorner(rows, cols).isOnes(0.));
}

BOOST_AUTO_TEST_SUITE_END()
