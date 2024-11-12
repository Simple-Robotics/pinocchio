//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/multibody/joint/joint-basic-visitors.hpp"
#include "pinocchio/multibody/joint/joint-generic.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

namespace bf = boost::fusion;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_check_joint_type)
{
  using namespace pinocchio;

  const JointModel jmodel_rx = JointModelRX();
  const JointModel jmodel_px = JointModelPX();

  typedef boost::mpl::vector<JointModelRX, JointModelRY, JointModelRZ> JointModelSequence;

  BOOST_CHECK(check_joint_type_withiin_sequence<JointModelSequence>(jmodel_rx) == true);
  BOOST_CHECK(check_joint_type_withiin_sequence<JointModelSequence>(jmodel_px) == false);
}

BOOST_AUTO_TEST_SUITE_END()
