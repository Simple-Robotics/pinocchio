//
// Copyright (c) 2023-2024 INRIA
//

#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/multibody/sample-models.hpp"

#include "constraints/init_constraints.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(constraint_variants)
{
  Model model;
  buildModels::humanoidRandom(model, true);

  Data data(model);

  RigidConstraintModel rcm = init_constraint<RigidConstraintModel>(model);
  RigidConstraintData rcd(rcm);

  ConstraintModel::ConstraintModelVariant constraint_model_variant = rcm;
  ConstraintModel constraint_model(rcm);
  ConstraintModel constraint_model_equal = rcm;

  ConstraintData constraint_data = rcm.createData();
}

BOOST_AUTO_TEST_CASE(constraint_visitors)
{
  Model model;
  buildModels::humanoidRandom(model, true);

  Data data(model);

  const SE3 M(SE3::Random());
  RigidConstraintModel rcm(CONTACT_3D, model, 0, M);
  RigidConstraintData rcd(rcm);
  BOOST_CHECK(ConstraintData(rcd) == ConstraintData(rcd));
  BOOST_CHECK(ConstraintData(rcd) == rcd);

  ConstraintModel constraint_model(rcm);

  // Test create data visitor
  {
    RigidConstraintData rcd(rcm);
    ConstraintData constraint_data = visitors::createData(constraint_model);
    constraint_data = rcd;
    BOOST_CHECK(constraint_data == rcd);
  }

  // Test calc visitor
  {
    ConstraintData constraint_data1(rcm.createData());
    visitors::calc(constraint_model, model, data, constraint_data1);
    rcm.calc(model, data, rcd);
    BOOST_CHECK(rcd == constraint_data1);
    ConstraintData constraint_data2(rcm.createData());
    constraint_model.calc(model, data, constraint_data2);
    BOOST_CHECK(rcd == constraint_data2);
  }

  // Test jacobian visitor
  {
    ConstraintData constraint_data(rcm.createData());
    Data::MatrixXs jacobian_matrix1 = Data::Matrix6x::Zero(6, model.nv),
                   jacobian_matrix2 = Data::Matrix6x::Zero(6, model.nv),
                   jacobian_matrix_ref = Data::Matrix6x::Zero(6, model.nv);
    rcm.jacobian(model, data, rcd, jacobian_matrix_ref);
    visitors::jacobian(constraint_model, model, data, constraint_data, jacobian_matrix1);
    BOOST_CHECK(jacobian_matrix1 == jacobian_matrix_ref);
    constraint_model.jacobian(model, data, constraint_data, jacobian_matrix2);
    BOOST_CHECK(jacobian_matrix2 == jacobian_matrix_ref);
  }
}

BOOST_AUTO_TEST_SUITE_END()
