//
// Copyright (c) 2024-2025 INRIA
//

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/constraints/joint-frictional-constraint.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/multibody/sample-models.hpp"

#include <iostream>

// Helpers
#include "constraints/jacobians-checker.hpp"

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
  //  FrictionalJointConstraintData constraint_data = constraint.createData();

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
        constraint.getRowActivableSparsityPattern(Eigen::DenseIndex(row_id));
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

BOOST_AUTO_TEST_CASE(cast)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  const Eigen::VectorXd q = neutral(model);

  const Data data(model);

  const std::string RF_name = "rleg6_joint";
  const JointIndex RF_id = model.getJointId(RF_name);

  const Model::IndexVector & RF_support = model.supports[RF_id];
  const Model::IndexVector active_joint_ids(RF_support.begin() + 1, RF_support.end());

  FrictionalJointConstraintModel cm(model, active_joint_ids);

  const auto cm_cast_double = cm.cast<double>();
  BOOST_CHECK(cm_cast_double == cm);

  const auto cm_cast_long_double = cm.cast<long double>();
  BOOST_CHECK(cm_cast_long_double.cast<double>() == cm);
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

  check_jacobians_operations(model, data, constraint_model, constraint_data);
}

BOOST_AUTO_TEST_CASE(constraint_coupling_inertia)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  const Eigen::VectorXd q = neutral(model);

  Data data(model);
  computeJointJacobians(model, data, q);

  const std::string RF_name = "rleg6_joint";
  const JointIndex RF_id = model.getJointId(RF_name);

  const Model::IndexVector & RF_support = model.supports[RF_id];
  const Model::IndexVector active_joint_ids(RF_support.begin() + 1, RF_support.end());

  FrictionalJointConstraintModel constraint_model(model, active_joint_ids);
  FrictionalJointConstraintData constraint_data(constraint_model);

  constraint_model.calc(model, data, constraint_data);
  const Eigen::VectorXd diagonal_inertia =
    Eigen::VectorXd::Random(constraint_model.size()).cwiseSquare();
  constraint_model.appendCouplingConstraintInertias(
    model, data, constraint_data, diagonal_inertia, WorldFrameTag());

  Eigen::Index row_id = 0;

  for (const auto joint_id : active_joint_ids)
  {
    //    std::cout << "joint_id: " << joint_id << std::endl;

    const auto & jmodel = model.joints[joint_id];
    const auto & jdata = data.joints[joint_id];
    const auto & oMjoint = data.oMi[joint_id];
    const auto jmodel_nv = jmodel.nv();
    const auto jmodel_idx_v = jmodel.idx_v();

    const auto diagonal_inertia_segment = diagonal_inertia.segment(row_id, jmodel_nv);

    const auto J_cols = data.J.middleCols(jmodel_idx_v, jmodel_nv);
    const Inertia::Matrix6 constraint_world_inertia_ref =
      J_cols * diagonal_inertia_segment.asDiagonal() * J_cols.transpose();
    const Eigen::MatrixXd constraint_world_inertia_ref_projected =
      J_cols.transpose() * constraint_world_inertia_ref * J_cols;

    //    std::cout << "diagonal_inertia_segment: " << diagonal_inertia_segment.transpose() <<
    //    std::endl; std::cout << "constraint_world_inertia_ref_projected:\n" <<
    //    constraint_world_inertia_ref_projected << std::endl;

    const auto S = jdata.S().matrix();
    const Inertia::Matrix6 constraint_local_inertia =
      S * diagonal_inertia_segment.asDiagonal() * S.transpose();
    const Inertia::Matrix6 constraint_world_inertia_ref2 =
      oMjoint.toActionMatrix() * constraint_local_inertia * oMjoint.inverse().toDualActionMatrix();

    const auto & constraint_world_inertia = data.oYaba_augmented[joint_id];
    BOOST_CHECK(constraint_world_inertia.isApprox(constraint_world_inertia_ref));
    BOOST_CHECK(constraint_world_inertia.isApprox(constraint_world_inertia_ref2));

    row_id += jmodel_nv;
    //    std::cout << "----" << std::endl;
  }

  Eigen::MatrixXd jacobian_matrix(constraint_model.size(), model.nv);
  constraint_model.jacobian(model, data, constraint_data, jacobian_matrix);

  Eigen::MatrixXd joint_space_constraint_inertia =
    jacobian_matrix.transpose() * diagonal_inertia.asDiagonal() * jacobian_matrix;

  //  std::cout << "diagonal_inertia: " << diagonal_inertia.transpose() << std::endl;
  //  std::cout << "joint_space_constraint_inertia:\n" << joint_space_constraint_inertia <<
  //  std::endl;

  for (JointIndex joint_id = 1; joint_id < JointIndex(model.njoints); ++joint_id)
  //  for(const auto joint_id: active_joint_ids)
  {
    //    std::cout << "joint_id: " << joint_id << std::endl;
    const auto & jmodel = model.joints[joint_id];
    const auto & jdata = data.joints[joint_id];

    const auto jmodel_nv = jmodel.nv();
    const auto jmodel_idx_v = jmodel.idx_v();

    const auto S = jdata.S().matrix();
    const auto J_cols = data.J.middleCols(jmodel_idx_v, jmodel_nv);
    const auto & oMjoint = data.oMi[joint_id];

    BOOST_CHECK((oMjoint.toActionMatrix() * S).isApprox(J_cols));
    BOOST_CHECK(
      (S.transpose() * oMjoint.inverse().toDualActionMatrix()).isApprox(J_cols.transpose()));
    BOOST_CHECK((oMjoint.toActionMatrixInverse() * J_cols).isApprox(S));

    const auto & constraint_world_inertia = data.oYaba_augmented[joint_id];
    const Inertia::Matrix6 constraint_local_inertia =
      oMjoint.inverse().toActionMatrix() * constraint_world_inertia * oMjoint.toDualActionMatrix();

    const Eigen::MatrixXd projected_constraint_local_inertia =
      S.transpose() * constraint_local_inertia * S;
    if (
      std::find(active_joint_ids.begin(), active_joint_ids.end(), joint_id)
      != active_joint_ids.end())
    {
      const Inertia::Matrix6 res_local = constraint_local_inertia
                                         - (constraint_local_inertia * S)
                                             * projected_constraint_local_inertia.inverse()
                                             * (S.transpose() * constraint_local_inertia);
      BOOST_CHECK(res_local.isZero());

      const Eigen::MatrixXd projected_constraint_world_inertia =
        J_cols.transpose() * constraint_world_inertia * J_cols;
      const Inertia::Matrix6 res_world = constraint_world_inertia
                                         - (constraint_world_inertia * J_cols)
                                             * projected_constraint_world_inertia.inverse()
                                             * (J_cols.transpose() * constraint_world_inertia);
      BOOST_CHECK(res_world.isZero());
    }

    const auto joint_space_constraint_inertia_block =
      joint_space_constraint_inertia.block(jmodel_idx_v, jmodel_idx_v, jmodel_nv, jmodel_nv);

    BOOST_CHECK(projected_constraint_local_inertia.isApprox(joint_space_constraint_inertia_block));
    //    BOOST_CHECK(projected_constraint_world_inertia.isApprox(joint_space_constraint_inertia_block));
    //    std::cout << "----" << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
