//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/constraints/frictional-joint-constraint.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/multibody/sample-models.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
typedef FrictionalJointConstraintModel::EigenIndexVector EigenIndexVector;
typedef FrictionalJointConstraintModel::BooleanVector BooleanVector;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(constraint_empty_constructor)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  const Data data(model);

  const Model::IndexVector empty_active_joint_ids;

  FrictionalJointConstraintModel constraint(model, empty_active_joint_ids);
}

BOOST_AUTO_TEST_CASE(constraint_constructor)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  const Data data(model);
  const auto & parents_fromRow = data.parents_fromRow;

  const std::string RF_name = "rleg6_joint";
  const JointIndex RF_id = model.getJointId(RF_name);

  const Model::IndexVector & RF_support = model.supports[RF_id];
  const Model::IndexVector active_joint_ids(RF_support.begin() + 1, RF_support.end());

  FrictionalJointConstraintModel constraint(model, active_joint_ids);

  // Check size
  {
    int total_size = 0;
    for (const JointIndex joint_id : active_joint_ids)
    {
      total_size += model.joints[joint_id].nv();
    }
    BOOST_CHECK(constraint.size() == total_size);
    BOOST_CHECK(constraint.getActiveDofs().size() == size_t(total_size));
  }

  // Check sparsity pattern
  {
    const EigenIndexVector & active_dofs = constraint.getActiveDofs();
    for (size_t row_id = 0; row_id < size_t(constraint.size()); ++row_id)
    {
      const Eigen::DenseIndex dof_id = active_dofs[row_id];
      const BooleanVector & row_sparsity_pattern =
        constraint.getRowSparsityPattern(Eigen::DenseIndex(row_id));
      const EigenIndexVector & row_active_indexes =
        constraint.getRowActiveIndexes(Eigen::DenseIndex(row_id));

      // Check that the rest of the indexes greater than dof_id are not active.
      BOOST_CHECK((row_sparsity_pattern.tail(model.nv - 1 - dof_id).array() == false).all());

      Eigen::DenseIndex id = dof_id;
      while (parents_fromRow[size_t(id)] > -1)
      {
        BOOST_CHECK(row_sparsity_pattern[id] == true);
        id = parents_fromRow[size_t(id)];
      }

      for (const Eigen::DenseIndex active_id : row_active_indexes)
      {
        BOOST_CHECK(row_sparsity_pattern[active_id] == true);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(constraint_jacobian)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  const Eigen::VectorXd q = neutral(model);

  const Data data(model);

  const std::string RF_name = "rleg6_joint";
  const JointIndex RF_id = model.getJointId(RF_name);

  const Model::IndexVector & RF_support = model.supports[RF_id];
  const Model::IndexVector active_joint_ids(RF_support.begin() + 1, RF_support.end());

  FrictionalJointConstraintModel constraint_model(model, active_joint_ids);
  FrictionalJointConstraintData constraint_data(constraint_model);

  Eigen::MatrixXd jacobian_matrix(constraint_model.size(), model.nv);
  constraint_model.jacobian(model, data, constraint_data, jacobian_matrix);

  const EigenIndexVector & active_dofs = constraint_model.getActiveDofs();
  for (Eigen::DenseIndex row_id = 0; row_id < constraint_model.size(); ++row_id)
  {
    const Eigen::DenseIndex dof_id = active_dofs[size_t(row_id)];
    BOOST_CHECK(jacobian_matrix.row(row_id).sum() == 1.);
    BOOST_CHECK(jacobian_matrix(row_id, dof_id) == 1.);
    BOOST_CHECK(
      (dof_id - 1) > 0 ? (jacobian_matrix.row(row_id).head(dof_id - 1).array() == 0).all() : true);
    BOOST_CHECK(
      (model.nv - dof_id - 1) > 0
        ? (jacobian_matrix.row(row_id).tail(model.nv - dof_id - 1).array() == 0).all()
        : true);
  }
}

BOOST_AUTO_TEST_SUITE_END()
