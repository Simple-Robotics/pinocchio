//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/utils/reference.hpp"

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_get_ref)
{
  using namespace ::pinocchio::helper;

  const double v_const = 10;
  BOOST_CHECK(&v_const == &get_ref(v_const));

  std::reference_wrapper<const double> v_const_ref = v_const;
  BOOST_CHECK(&v_const == &get_ref(v_const_ref));

  std::shared_ptr<const double> v_const_ptr = std::make_shared<const double>(v_const);
  BOOST_CHECK(v_const_ptr.get() != &get_ref(v_const_ref));
  BOOST_CHECK(v_const_ptr.get() == &get_ref(v_const_ptr));
}

BOOST_AUTO_TEST_SUITE_END()
