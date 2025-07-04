//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/container/storage.hpp"

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
typedef EigenStorageTpl<Eigen::MatrixXd> EigenStorageMatrix;
typedef EigenStorageTpl<PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(Eigen::MatrixXd)>
  EigenStorageRowMatrix;
typedef EigenStorageTpl<Eigen::VectorXd> EigenStorageVector;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(eigen_storage_default)
{
  EigenStorageMatrix storage_matrix;
  BOOST_CHECK(!storage_matrix.isValid());
  EigenStorageVector storage_vector;
  BOOST_CHECK(!storage_vector.isValid());
}

BOOST_AUTO_TEST_CASE(eigen_storage_matrix)
{
  const Eigen::DenseIndex rows = 10, cols = 20;

  const Eigen::DenseIndex initial_capacity = rows * cols;
  EigenStorageMatrix storage(rows, cols);
  BOOST_CHECK(storage.isValid());

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

  // Check copy
  EigenStorageMatrix storage_copy(storage);
  BOOST_CHECK(storage_copy.data() != storage.data());
  BOOST_CHECK(storage_copy.map() == storage.map());
  BOOST_CHECK(storage_copy.capacity() == storage.capacity());
  BOOST_CHECK(storage_copy.storage() == storage.storage());

  // Check resize
  const Eigen::DenseIndex new_rows = 2 * rows, new_cols = cols;
  storage.conservativeResize(new_rows, new_cols);
  BOOST_CHECK(matrix_map.data() == storage.data());
  BOOST_CHECK(storage.map().topLeftCorner(rows, cols).isOnes(0.));
}

BOOST_AUTO_TEST_CASE(cast)
{
  const Eigen::DenseIndex rows = 10, cols = 20;

  EigenStorageMatrix storage(rows, cols);
  storage.map().setConstant(1.895);

  const auto storage_cast_double = storage.cast<double>();
  BOOST_CHECK(storage_cast_double.map() == storage.map());

  const auto storage_cast_long_double = storage.cast<long double>();
  BOOST_CHECK(storage_cast_long_double.cast<double>().map() == storage.map());
}

BOOST_AUTO_TEST_CASE(eigen_storage_row_matrix)
{
  const Eigen::DenseIndex rows = 10, cols = 20;

  const Eigen::DenseIndex initial_capacity = rows * cols;
  EigenStorageRowMatrix storage(rows, cols);
  BOOST_CHECK(storage.isValid());

  BOOST_CHECK(storage.capacity() == initial_capacity);
  BOOST_CHECK(storage.rows() == rows);
  BOOST_CHECK(storage.cols() == cols);

  EigenStorageRowMatrix::RefMapType matrix_map = storage.map();
  BOOST_CHECK(matrix_map.data() == storage.data());

  matrix_map.setIdentity();
  BOOST_CHECK(storage.map().isIdentity(0.));
  BOOST_CHECK(static_cast<const EigenStorageRowMatrix &>(storage).map().isIdentity(0.));
  matrix_map.setOnes();
  BOOST_CHECK(storage.map().isOnes(0.));
  BOOST_CHECK(static_cast<const EigenStorageRowMatrix &>(storage).map().isOnes(0.));

  // Check copy
  EigenStorageRowMatrix storage_copy(storage);
  BOOST_CHECK(storage_copy.data() != storage.data());
  BOOST_CHECK(storage_copy.map() == storage.map());
  BOOST_CHECK(storage_copy.capacity() == storage.capacity());
  BOOST_CHECK(storage_copy.storage() == storage.storage());

  // Check resize
  const Eigen::DenseIndex new_rows = 2 * rows, new_cols = cols;
  storage.conservativeResize(new_rows, new_cols);
  BOOST_CHECK(matrix_map.data() == storage.data());
  BOOST_CHECK(storage.map().topLeftCorner(rows, cols).isOnes(0.));
}

BOOST_AUTO_TEST_CASE(eigen_storage_vector)
{
  const Eigen::DenseIndex size = 100;

  const Eigen::DenseIndex initial_capacity = size;
  EigenStorageVector storage(size);
  BOOST_CHECK(storage.isValid());

  BOOST_CHECK(storage.capacity() == initial_capacity);
  BOOST_CHECK(storage.rows() == size);
  BOOST_CHECK(storage.cols() == 1);

  EigenStorageVector::RefMapType vector_map = storage.map();
  BOOST_CHECK(vector_map.data() == storage.data());

  vector_map.setOnes();
  BOOST_CHECK(storage.map().isOnes(0.));
  BOOST_CHECK(static_cast<const EigenStorageVector &>(storage).map().isOnes(0.));

  // Check copy
  EigenStorageVector storage_copy(storage);
  BOOST_CHECK(storage_copy.data() != storage.data());
  BOOST_CHECK(storage_copy.map() == storage.map());
  BOOST_CHECK(storage_copy.capacity() == storage.capacity());
  BOOST_CHECK(storage_copy.storage() == storage.storage());

  // Check resize
  const Eigen::DenseIndex new_size = 2 * size;
  storage.conservativeResize(new_size);
  BOOST_CHECK(vector_map.data() == storage.data());
  BOOST_CHECK(storage.map().head(size).isOnes(0.));

  storage.conservativeResize(new_size, 1);
  BOOST_CHECK(vector_map.data() == storage.data());
  BOOST_CHECK(storage.map().head(size).isOnes(0.));

  //  storage.conservativeResize(1,new_size); // will not pass
}

BOOST_AUTO_TEST_SUITE_END()
