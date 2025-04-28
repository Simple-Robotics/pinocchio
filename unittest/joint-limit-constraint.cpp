//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/constraints/joint-limit-constraint.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/multibody/sample-models.hpp"

// Helpers
#include "constraints/jacobians-checker.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
typedef JointLimitConstraintModel::EigenIndexVector EigenIndexVector;
typedef JointLimitConstraintModel::BooleanVector BooleanVector;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(constraint_empty_constructor)
{
  pinocchio::Model model;
  pinocchio::buildModels::manipulator(model);

  const Data data(model);

  const Model::IndexVector empty_activable_joint_ids;

  JointLimitConstraintModel constraint(model, empty_activable_joint_ids);
}

BOOST_AUTO_TEST_CASE(constraint_constructor_with_infinite_bounds)
{
  pinocchio::Model model;
  pinocchio::buildModels::manipulator(model);

  model.lowerPositionLimit.fill(-std::numeric_limits<double>::max());
  model.upperPositionLimit.fill(+std::numeric_limits<double>::max());

  const Data data(model);

  const std::string ee_name = "wrist2_joint";
  const JointIndex ee_id = model.getJointId(ee_name);

  const Model::IndexVector & ee_support = model.supports[ee_id];
  const Model::IndexVector activable_joint_ids(ee_support.begin() + 1, ee_support.end());
  JointLimitConstraintModel constraint(model, activable_joint_ids);

  BOOST_CHECK(constraint.size() == 0);
}

BOOST_AUTO_TEST_CASE(constraint_constructor)
{
  pinocchio::Model model;
  pinocchio::buildModels::manipulator(model);

  Data data(model);
  const auto & parents_fromRow = data.parents_fromRow;

  //  std::cout << "model:\n" << model << std::endl;

  const std::string ee_name = "wrist2_joint";
  const JointIndex ee_id = model.getJointId(ee_name);

  const Model::IndexVector & ee_support = model.supports[ee_id];
  const Model::IndexVector activable_joint_ids(ee_support.begin() + 1, ee_support.end());

  for (const JointIndex joint_id : activable_joint_ids)
  {
    const auto & jmodel = model.joints[joint_id];
    std::cout << "joint type: " << jmodel.shortname() << std::endl;
  }
  JointLimitConstraintModel constraint(model, activable_joint_ids);

  // Check size
  {
    int total_size = 0;
    for (const JointIndex joint_id : activable_joint_ids)
    {
      const auto & jmodel = model.joints[joint_id];
      const int idx_q = jmodel.idx_q();
      const int nq = jmodel.nq();
      for (int k = 0; k < nq; ++k)
      {
        const int index_q = idx_q + k;
        if (model.lowerPositionLimit[index_q] != -std::numeric_limits<double>::max())
          total_size++;
        if (model.upperPositionLimit[index_q] != +std::numeric_limits<double>::max())
          total_size++;
      }
    }
    BOOST_CHECK(constraint.size() == total_size);
  }

  // we set the margin to infinity so all limits are taken into account in what follows.
  model.positionLimitMargin =
    Eigen::VectorXd::Constant(constraint.size(), std::numeric_limits<double>::max());
  constraint = JointLimitConstraintModel(model, activable_joint_ids);
  JointLimitConstraintData constraint_data(constraint);
  const Eigen::VectorXd q0 = randomConfiguration(model);
  data.q_in = q0;
  constraint.resize(model, data, constraint_data);
  constraint.calc(model, data, constraint_data);

  // Check sparsity pattern
  // {
  //   // const EigenIndexVector & active_dofs_lower =
  //   // constraint.getActiveLowerBoundConstraintsTangent(); const EigenIndexVector &
  //   // active_dofs_upper = constraint.getActiveUpperBoundConstraintsTangent();
  //   EigenIndexVector active_dofs(active_dofs_lower);
  //   active_dofs.insert(active_dofs.end(), active_dofs_upper.begin(), active_dofs_upper.end());

  //   for (size_t row_id = 0; row_id < size_t(constraint.activeSize()); ++row_id)
  //   {
  //     const Eigen::DenseIndex dof_id = active_dofs[row_id];
  //     const BooleanVector & row_sparsity_pattern =
  //       constraint.getRowActiveSparsityPattern(Eigen::DenseIndex(row_id));
  //     const EigenIndexVector & row_active_indexes =
  //       constraint.getRowActiveIndexes(Eigen::DenseIndex(row_id));

  //     // Check that the rest of the indexes greater than dof_id are not active.
  //     BOOST_CHECK((row_sparsity_pattern.tail(model.nv - 1 - dof_id).array() == false).all());

  //     Eigen::DenseIndex id = dof_id;
  //     while (parents_fromRow[size_t(id)] > -1)
  //     {
  //       BOOST_CHECK(row_sparsity_pattern[id] == true);
  //       id = parents_fromRow[size_t(id)];
  //     }

  //     for (const Eigen::DenseIndex active_id : row_active_indexes)
  //     {
  //       BOOST_CHECK(row_sparsity_pattern[active_id] == true);
  //     }
  //   }
  // }

  // Check projection on force sets
  {
    const Eigen::DenseIndex nb_lower_active_dofs = Eigen::DenseIndex(constraint.lowerActiveSize());
    const Eigen::DenseIndex nb_upper_active_dofs = Eigen::DenseIndex(constraint.upperActiveSize());

    const Eigen::DenseIndex total_active_dofs = nb_lower_active_dofs + nb_upper_active_dofs;

    const int num_projections = int(1e6);
    for (int k = 0; k < num_projections; ++k)
    {
      const Eigen::VectorXd f = Eigen::VectorXd::Random(total_active_dofs);
      const Eigen::VectorXd f_proj = constraint.set().project(f);

      BOOST_CHECK((f_proj.head(nb_lower_active_dofs).array() <= 0).all());
      BOOST_CHECK((f_proj.tail(nb_lower_active_dofs).array() >= 0).all());
    }
  }
}

BOOST_AUTO_TEST_CASE(cast)
{
  pinocchio::Model model;
  pinocchio::buildModels::manipulator(model);

  const Data data(model);

  const std::string ee_name = "wrist2_joint";
  const JointIndex ee_id = model.getJointId(ee_name);

  const Model::IndexVector & ee_support = model.supports[ee_id];
  const Model::IndexVector active_joint_ids(ee_support.begin() + 1, ee_support.end());
  JointLimitConstraintModel cm(model, active_joint_ids);

  const auto cm_cast_double = cm.cast<double>();
  BOOST_CHECK(cm_cast_double == cm);

  const auto cm_cast_long_double = cm.cast<long double>();
  BOOST_CHECK(cm_cast_long_double.cast<double>() == cm);
}

BOOST_AUTO_TEST_CASE(constraint_jacobian)
{
  pinocchio::Model model;
  pinocchio::buildModels::manipulator(model);

  Data data(model);

  const std::string ee_name = "wrist2_joint";
  const JointIndex ee_id = model.getJointId(ee_name);

  const Model::IndexVector & ee_support = model.supports[ee_id];
  const Model::IndexVector activable_joint_ids(ee_support.begin() + 1, ee_support.end());

  JointLimitConstraintModel constraint_model(model, activable_joint_ids);
  JointLimitConstraintData constraint_data(constraint_model);

  // we set the margin to infinity so all limits are taken into account in what follows.
  model.positionLimitMargin =
    Eigen::VectorXd::Constant(constraint_model.size(), std::numeric_limits<double>::max());

  // Check against finite differences on the drift of the constraint
  const double eps_fd = 1e-8;
  for (int i = 0; i < 1e4; ++i)
  {
    const Eigen::VectorXd q0 = randomConfiguration(model);
    data.q_in = q0;
    constraint_model.resize(model, data, constraint_data);
    constraint_model.calc(model, data, constraint_data);
    Eigen::MatrixXd jacobian_matrix(constraint_model.activeSize(), model.nv);
    constraint_model.jacobian(model, data, constraint_data, jacobian_matrix);
    Data data_fd(model);
    JointLimitConstraintData constraint_data_fd(constraint_model);
    Eigen::MatrixXd jacobian_matrix_fd(constraint_model.activeSize(), model.nv);
    // TODO compute jacobian only for activable constraints
    for (Eigen::DenseIndex k = 0; k < model.nv; ++k)
    {
      Eigen::VectorXd v_eps = Eigen::VectorXd::Zero(model.nv);
      v_eps[k] = eps_fd;
      const Eigen::VectorXd q_plus = integrate(model, q0, v_eps);
      data_fd.q_in = q_plus;
      constraint_model.resize(model, data_fd, constraint_data_fd);
      constraint_model.calc(model, data_fd, constraint_data_fd);

      jacobian_matrix_fd.col(k) =
        (constraint_data_fd.constraint_residual - constraint_data.constraint_residual) / eps_fd;
    }

    BOOST_CHECK(jacobian_matrix.isApprox(jacobian_matrix_fd, math::sqrt(eps_fd)));
  }

  check_jacobians_operations(model, data, constraint_model, constraint_data);
}

BOOST_AUTO_TEST_CASE(dynamic_constraint_residual)
{
  Model model;
  JointIndex joint_id_x = model.addJoint(0, JointModelPX(), SE3::Identity(), "slider_x");
  JointIndex joint_id_y = model.addJoint(joint_id_x, JointModelPY(), SE3::Identity(), "slider_y");
  JointIndex joint_id_z = model.addJoint(joint_id_y, JointModelPZ(), SE3::Identity(), "slider_z");

  const SE3::Vector3 small_box_dims = SE3::Vector3::Ones() * 1e-3;
  const double small_box_mass = 1e-6;
  const Inertia small_box_inertia =
    Inertia::FromBox(small_box_mass, small_box_dims[0], small_box_dims[1], small_box_dims[2]);

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const double box_mass = 10;
  const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);

  model.appendBodyToJoint(joint_id_x, small_box_inertia);
  model.appendBodyToJoint(joint_id_y, small_box_inertia);
  model.appendBodyToJoint(joint_id_z, box_inertia);
  model.gravity.setZero();
  model.lowerPositionLimit[0] = 0.;
  model.lowerPositionLimit[1] = 0.;
  model.lowerPositionLimit[2] = 0.;
  // We deactivate the upper limits
  model.upperPositionLimit[0] = std::numeric_limits<double>::max();
  model.upperPositionLimit[1] = std::numeric_limits<double>::max();
  model.upperPositionLimit[2] = std::numeric_limits<double>::max();

  Data data(model);

  const Model::IndexVector activable_joint_ids = {joint_id_x, joint_id_y, joint_id_z};
  ;

  JointLimitConstraintModel constraint_model(model, activable_joint_ids);
  BOOST_CHECK(constraint_model.size() == model.nq);
  JointLimitConstraintData constraint_data(constraint_model);

  model.positionLimitMargin = Eigen::VectorXd::Constant(constraint_model.size(), 1e-3);

  for (int i = 0; i < 1e4; ++i)
  {
    const Eigen::VectorXd q0 = Eigen::VectorXd::Random(model.nq);
    std::size_t active_size = 0;
    std::vector<std::size_t> active_indexes;
    Eigen::VectorXd activable_residual(constraint_model.size());
    for (int j = 0; j < model.nq; j++)
    {
      activable_residual(j) = -(q0(j) - model.lowerPositionLimit[j]);
      if (-activable_residual(j) < model.positionLimitMargin[j])
      {
        active_size++;
        active_indexes.push_back((std::size_t)j);
      }
    }
    Eigen::VectorXd residual(active_size);
    for (std::size_t j = 0; j < active_size; j++)
    {
      residual((Eigen::Index)j) = activable_residual((Eigen::Index)active_indexes[j]);
    }
    data.q_in = q0;
    constraint_model.resize(model, data, constraint_data);
    constraint_model.calc(model, data, constraint_data);
    BOOST_CHECK((int)active_size == constraint_model.activeSize());
    BOOST_CHECK((int)active_size == constraint_data.constraint_residual.size());
    BOOST_CHECK(constraint_data.constraint_residual.isApprox(residual));
    BOOST_CHECK(constraint_data.activable_constraint_residual.isApprox(activable_residual));
    BOOST_CHECK(active_indexes == constraint_model.getActiveSetIndexes());
  }
}

BOOST_AUTO_TEST_CASE(dynamic_constraint_jacobian)
{
  Model model;
  JointIndex joint_id_x = model.addJoint(0, JointModelPX(), SE3::Identity(), "slider_x");
  JointIndex joint_id_y = model.addJoint(joint_id_x, JointModelPY(), SE3::Identity(), "slider_y");
  JointIndex joint_id_z = model.addJoint(joint_id_y, JointModelPZ(), SE3::Identity(), "slider_z");

  const SE3::Vector3 small_box_dims = SE3::Vector3::Ones() * 1e-3;
  const double small_box_mass = 1e-6;
  const Inertia small_box_inertia =
    Inertia::FromBox(small_box_mass, small_box_dims[0], small_box_dims[1], small_box_dims[2]);

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const double box_mass = 10;
  const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);

  model.appendBodyToJoint(joint_id_x, small_box_inertia);
  model.appendBodyToJoint(joint_id_y, small_box_inertia);
  model.appendBodyToJoint(joint_id_z, box_inertia);
  model.gravity.setZero();
  model.lowerPositionLimit[0] = 0.;
  model.lowerPositionLimit[1] = 0.;
  model.lowerPositionLimit[2] = 0.;
  // We deactivate the upper limits
  model.upperPositionLimit[0] = std::numeric_limits<double>::max();
  model.upperPositionLimit[1] = std::numeric_limits<double>::max();
  model.upperPositionLimit[2] = std::numeric_limits<double>::max();

  Data data(model);

  const Model::IndexVector activable_joint_ids = {joint_id_x, joint_id_y, joint_id_z};
  ;

  JointLimitConstraintModel constraint_model(model, activable_joint_ids);
  JointLimitConstraintData constraint_data(constraint_model);

  model.positionLimitMargin = Eigen::VectorXd::Constant(constraint_model.size(), 1e-3);

  // Check against finite differences on the drift of the constraint
  const double eps_fd = 1e-8;
  for (int i = 0; i < 1e4; ++i)
  {
    const Eigen::VectorXd q0 = Eigen::VectorXd::Random(model.nq);
    int active_size = 0;
    std::vector<Eigen::DenseIndex> active_indexes;
    for (int j = 0; j < model.nq; j++)
    {
      if (q0(j) - model.lowerPositionLimit[j] < model.positionLimitMargin[j])
      {
        active_size++;
        active_indexes.push_back(j);
      }
    }
    data.q_in = q0;
    constraint_model.resize(model, data, constraint_data);
    constraint_model.calc(model, data, constraint_data);
    std::vector<std::size_t> active_set_indexes = constraint_model.getActiveSetIndexes();
    BOOST_CHECK(active_size == constraint_model.activeSize());
    BOOST_CHECK(active_size == constraint_data.constraint_residual.size());
    Eigen::MatrixXd jacobian_matrix(constraint_model.activeSize(), model.nv);
    constraint_model.jacobian(model, data, constraint_data, jacobian_matrix);
    Data data_fd(model);
    JointLimitConstraintData constraint_data_fd(constraint_model);
    Eigen::MatrixXd jacobian_matrix_fd(constraint_model.activeSize(), model.nv);
    // For now, we assume model.nq == model.nv
    assert(model.nq == model.nv);
    for (Eigen::DenseIndex k = 0; k < model.nv; ++k)
    {
      Eigen::VectorXd v_eps = Eigen::VectorXd::Zero(model.nv);
      v_eps[k] = eps_fd;
      const Eigen::VectorXd q_plus = integrate(model, q0, v_eps);
      data_fd.q_in = q_plus;
      constraint_model.resize(model, data_fd, constraint_data_fd);
      constraint_model.calc(model, data_fd, constraint_data_fd);
      bool same_active_set = active_set_indexes == constraint_model.getActiveSetIndexes();
      // if the active set is identical we can check the jacobian
      if (!same_active_set)
        continue;
      jacobian_matrix_fd.col(k) =
        (constraint_data_fd.constraint_residual - constraint_data.constraint_residual) / eps_fd;
    }

    BOOST_CHECK(jacobian_matrix.isApprox(jacobian_matrix_fd, math::sqrt(eps_fd)));
  }
}

BOOST_AUTO_TEST_SUITE_END()
