//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/container/double-entry-container.hpp"

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
typedef Eigen::Matrix<double, 6, 6> Matrix6;
typedef Eigen::aligned_allocator<Matrix6> Allocator;

typedef container::DoubleEntryContainer<Matrix6, Allocator> Container;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_all)
{
  const Eigen::Index nrows = 20, ncols = 20;

  Container container(nrows, ncols);
  BOOST_CHECK(container.size() == 0);
  BOOST_CHECK(container.rows() == nrows);
  BOOST_CHECK(container.cols() == ncols);
  BOOST_CHECK(container.begin() == container.end());

  container.reserve(size_t(nrows));
  BOOST_CHECK(container.values().capacity() == size_t(nrows));

  // Fill diagonal with random values
  for (Eigen::Index id = 0; id < nrows; id++)
  {
    const bool res = container.insert(id, id, Matrix6::Constant(double(id)));
    BOOST_CHECK(res);
    BOOST_CHECK(container.size() == size_t(id + 1));
  }

  for (Eigen::Index id = 0; id < nrows; id++)
  {
    const bool res = container.insert(id, id, Matrix6::Constant(double(id)));
    BOOST_CHECK(!res);
    BOOST_CHECK(container.size() == size_t(nrows));
  }

  {
    const Eigen::VectorXi linear_range = Eigen::VectorXi::LinSpaced(nrows, 0, nrows - 1);
    BOOST_CHECK(container.keys().matrix().diagonal() == linear_range.cast<long>());
  }

  // Set diagonal and check
  container.fill(Matrix6::Identity());
  for (auto val : container)
    BOOST_CHECK(val == Matrix6::Identity());

  // Check method find
  for (Eigen::Index row_id = 0; row_id < nrows; row_id++)
  {
    for (Eigen::Index col_id = 0; col_id < ncols; col_id++)
    {
      if (row_id == col_id)
        BOOST_CHECK(*container.find(row_id, col_id) == Matrix6::Identity());
      else
        BOOST_CHECK(container.find(row_id, col_id) == container.end());
    }
  }

  // Remove elt (4,4)
  Container copy(container);
  BOOST_CHECK(!container.remove(3, 4));
  BOOST_CHECK(container.find(4, 4) != container.end());
  BOOST_CHECK(container.remove(4, 4));
  BOOST_CHECK(container.size() == size_t(nrows - 1));
  BOOST_CHECK(container.find(4, 4) == container.end());

  //
}

BOOST_AUTO_TEST_SUITE_END()
