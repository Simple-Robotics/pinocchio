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

  FrictionalPointConstraintModel rcm = init_constraint<FrictionalPointConstraintModel>(model);
  FrictionalPointConstraintData rcd(rcm);

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

  FrictionalPointConstraintModel rcm = init_constraint<FrictionalPointConstraintModel>(model);
  FrictionalPointConstraintData rcd(rcm);
  BOOST_CHECK(ConstraintData(rcd) == ConstraintData(rcd));
  BOOST_CHECK(ConstraintData(rcd) == rcd);

  ConstraintModel constraint_model(rcm);

  // Test size
  {
    BOOST_CHECK(constraint_model.size() == rcm.size());
  }

  // Test create data visitor
  {
    FrictionalPointConstraintData rcd(rcm);
    ConstraintData constraint_data = visitors::createData(constraint_model);
    constraint_data = rcd;
    BOOST_CHECK(constraint_data == rcd);
  }

  // Test calc visitor
  {
    ConstraintData constraint_data1(rcm.createData());
    visitors::calc(constraint_model, model, data, constraint_data1);
    rcm.calc(model, data, rcd);
    //    BOOST_CHECK(rcd == constraint_data1);
    ConstraintData constraint_data2(rcm.createData());
    constraint_model.calc(model, data, constraint_data2);
    //    BOOST_CHECK(rcd == constraint_data2);
  }

  // Test jacobian visitor
  {
    ConstraintData constraint_data(rcm.createData());
    Data::MatrixXs jacobian_matrix1 = Data::MatrixXs::Zero(rcm.size(), model.nv),
                   jacobian_matrix2 = Data::MatrixXs::Zero(rcm.size(), model.nv),
                   jacobian_matrix_ref = Data::MatrixXs::Zero(rcm.size(), model.nv);
    rcm.jacobian(model, data, rcd, jacobian_matrix_ref);
    visitors::jacobian(constraint_model, model, data, constraint_data, jacobian_matrix1);
    BOOST_CHECK(jacobian_matrix1 == jacobian_matrix_ref);
    constraint_model.jacobian(model, data, constraint_data, jacobian_matrix2);
    BOOST_CHECK(jacobian_matrix2 == jacobian_matrix_ref);
  }

  // Test getRowActiveIndexes
  {
    for (Eigen::DenseIndex row_id = 0; row_id < constraint_model.size(); ++row_id)
    {
      BOOST_CHECK(constraint_model.getRowActiveIndexes(row_id) == rcm.getRowActiveIndexes(row_id));
    }
  }

  // Test getRowSparsityPattern
  {
    for (Eigen::DenseIndex row_id = 0; row_id < constraint_model.size(); ++row_id)
    {
      BOOST_CHECK(
        constraint_model.getRowSparsityPattern(row_id) == rcm.getRowSparsityPattern(row_id));
    }
  }

  // Test jacobianMatrixProduct
  {
    const Eigen::Index num_cols = 20;
    ConstraintData constraint_data(rcm.createData());
    const Data::MatrixXs input_matrix = Data::MatrixXs::Random(model.nv, num_cols);
    Data::MatrixXs output_matrix1(rcm.size(), num_cols), output_matrix2(rcm.size(), num_cols),
      output_matrix_ref(rcm.size(), num_cols);
    rcm.jacobianMatrixProduct(model, data, rcd, input_matrix, output_matrix_ref);
    visitors::jacobianMatrixProduct(
      constraint_model, model, data, constraint_data, input_matrix, output_matrix1);
    BOOST_CHECK(output_matrix1 == output_matrix_ref);
    constraint_model.jacobianMatrixProduct(
      model, data, constraint_data, input_matrix, output_matrix2);
    BOOST_CHECK(output_matrix2 == output_matrix_ref);
  }
}

BOOST_AUTO_TEST_SUITE_END()
