//
// Copyright (c) 2024 INRIA
//

#include "utils/model-generator.hpp"
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
  Model model;
  buildAllJointsModel(model);

  const Data data(model);

  const Model::IndexVector empty_activable_joint_ids;

  JointLimitConstraintModel constraint(model, empty_activable_joint_ids);
}

BOOST_AUTO_TEST_CASE(constraint_constructor_with_infinite_bounds)
{
  Model model;
  buildAllJointsModel(model);

  model.lowerPositionLimit.fill(-std::numeric_limits<double>::max());
  model.upperPositionLimit.fill(+std::numeric_limits<double>::max());

  const Data data(model);

  const std::string ee_name = "translation_joint";
  const JointIndex ee_id = model.getJointId(ee_name);

  const Model::IndexVector & ee_support = model.supports[ee_id];
  const Model::IndexVector activable_joint_ids(ee_support.begin() + 1, ee_support.end());
  JointLimitConstraintModel constraint(model, activable_joint_ids);

  BOOST_CHECK(constraint.size() == 0);
}

BOOST_AUTO_TEST_CASE(constraint_constructor)
{
  Model model;
  buildAllJointsModel(model);

  Data data(model);
  const auto & parents_fromRow = data.parents_fromRow;

  const std::string ee_name = "translation_joint";
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
      const auto has_configuration_limit = jmodel.hasConfigurationLimit();
      for (int k = 0; k < nq; ++k)
      {
        const int index_q = idx_q + k;
        if (has_configuration_limit[k])
        {
          if (model.lowerPositionLimit[index_q] != -std::numeric_limits<double>::max())
            total_size++;
          if (model.upperPositionLimit[index_q] != +std::numeric_limits<double>::max())
            total_size++;
        }
      }
    }
    BOOST_CHECK(constraint.size() == total_size);
  }

  // we set the margin to infinity so all limits are taken into account in what follows.
  model.positionLimitMargin =
    Eigen::VectorXd::Constant(model.nq, std::numeric_limits<double>::max());
  constraint = JointLimitConstraintModel(model, activable_joint_ids);
  JointLimitConstraintData constraint_data(constraint);
  const Eigen::VectorXd q0 = randomConfiguration(model);
  data.q_in = q0;
  constraint.resize(model, data, constraint_data);
  constraint.calc(model, data, constraint_data);

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
  Model model;
  buildAllJointsModel(model);

  const Data data(model);

  const std::string ee_name = "translation_joint";
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
  Model model;
  buildAllJointsModel(model);

  Data data(model);

  const std::string ee_name = "translation_joint";
  const JointIndex ee_id = model.getJointId(ee_name);

  const Model::IndexVector & ee_support = model.supports[ee_id];
  const Model::IndexVector activable_joint_ids(ee_support.begin() + 1, ee_support.end());

  // we set the margin to infinity so all limits are taken into account in what follows.
  model.positionLimitMargin =
    Eigen::VectorXd::Constant(model.nq, std::numeric_limits<double>::max());

  JointLimitConstraintModel constraint_model(model, activable_joint_ids);
  JointLimitConstraintData constraint_data(constraint_model);

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
  buildAllJointsModel(model);

  model.gravity.setZero();
  model.lowerPositionLimit.fill(0.);
  const Eigen::VectorXd qmin(Eigen::VectorXd::Zero(model.nq));
  const Eigen::VectorXd qmax(Eigen::VectorXd::Ones(model.nq));
  // We deactivate the upper limits
  model.upperPositionLimit.fill(std::numeric_limits<double>::max());
  model.positionLimitMargin = Eigen::VectorXd::Constant(model.nq, .5);

  Data data(model);

  int total_size = 0;
  Model::IndexVector activable_joint_ids;
  for (Model::JointIndex i = 1; i < (Model::JointIndex)model.njoints; ++i)
  {
    activable_joint_ids.push_back(i);

    const auto & jmodel = model.joints[i];
    const int nq = jmodel.nq();
    const auto has_configuration_limit = jmodel.hasConfigurationLimit();
    for (int k = 0; k < nq; ++k)
    {
      if (has_configuration_limit[size_t(k)])
        total_size++;
    }
  }

  JointLimitConstraintModel constraint_model(model, activable_joint_ids);
  BOOST_CHECK(constraint_model.size() == total_size);
  JointLimitConstraintData constraint_data(constraint_model);

  for (int i = 0; i < 1e4; ++i)
  {
    const Eigen::VectorXd q0 = randomConfiguration(model, qmin, qmax);
    std::size_t active_size = 0;
    std::vector<std::size_t> active_indexes;
    Eigen::VectorXd activable_residual(constraint_model.size());

    int set_idx = 0;
    for (Model::JointIndex i = 1; i < (Model::JointIndex)model.njoints; ++i)
    {
      const auto & jmodel = model.joints[i];
      const int idx_q = jmodel.idx_q();
      const int nq = jmodel.nq();

      const auto has_configuration_limit = jmodel.hasConfigurationLimit();
      for (int k = 0; k < nq; ++k)
      {
        if (!has_configuration_limit[size_t(k)])
          continue;
        const int q_index = idx_q + k;
        activable_residual(set_idx) = -(q0(q_index) - model.lowerPositionLimit[q_index]);
        if (-activable_residual(set_idx) < model.positionLimitMargin[q_index])
        {
          active_size++;
          active_indexes.push_back((std::size_t)set_idx);
        }
        set_idx++;
      }
    }
    BOOST_CHECK(constraint_model.size() == set_idx);

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
  buildAllJointsModel(model);

  model.gravity.setZero();
  model.lowerPositionLimit.fill(0.);
  const Eigen::VectorXd qmin(Eigen::VectorXd::Zero(model.nq));
  const Eigen::VectorXd qmax(Eigen::VectorXd::Ones(model.nq));
  // We deactivate the upper limits
  model.upperPositionLimit.fill(std::numeric_limits<double>::max());
  model.positionLimitMargin = Eigen::VectorXd::Constant(model.nq, .5);

  Data data(model);

  int total_size = 0;
  Model::IndexVector activable_joint_ids;
  for (Model::JointIndex i = 1; i < (Model::JointIndex)model.njoints; ++i)
  {
    activable_joint_ids.push_back(i);

    const auto & jmodel = model.joints[i];
    const int nq = jmodel.nq();
    const auto has_configuration_limit = jmodel.hasConfigurationLimit();
    for (int k = 0; k < nq; ++k)
    {
      if (has_configuration_limit[size_t(k)])
        total_size++;
    }
  }

  JointLimitConstraintModel constraint_model(model, activable_joint_ids);
  BOOST_CHECK(constraint_model.size() == total_size);
  JointLimitConstraintData constraint_data(constraint_model);

  const double eps_fd = 1e-8;
  for (int i = 0; i < 1e4; ++i)
  {
    const Eigen::VectorXd q0 = randomConfiguration(model, qmin, qmax);
    int active_size = 0;
    std::vector<std::size_t> active_indexes;
    Eigen::VectorXd activable_residual(constraint_model.size());

    int set_idx = 0;
    for (Model::JointIndex i = 1; i < (Model::JointIndex)model.njoints; ++i)
    {
      const auto & jmodel = model.joints[i];
      const int idx_q = jmodel.idx_q();
      const int nq = jmodel.nq();

      const auto has_configuration_limit = jmodel.hasConfigurationLimit();
      for (int k = 0; k < nq; ++k)
      {
        if (!has_configuration_limit[size_t(k)])
          continue;
        const int q_index = idx_q + k;
        activable_residual(set_idx) = -(q0(q_index) - model.lowerPositionLimit[q_index]);
        if (-activable_residual(set_idx) < model.positionLimitMargin[q_index])
        {
          active_size++;
          active_indexes.push_back((std::size_t)set_idx);
        }
        set_idx++;
      }
    }
    BOOST_CHECK(constraint_model.size() == set_idx);

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

    bool ok_to_check = true;
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
      {
        ok_to_check = false;
        continue;
      }
      jacobian_matrix_fd.col(k) =
        (constraint_data_fd.constraint_residual - constraint_data.constraint_residual) / eps_fd;
      BOOST_CHECK(jacobian_matrix.col(k).isApprox(jacobian_matrix_fd.col(k), math::sqrt(eps_fd)));
    }

    if (ok_to_check)
      BOOST_CHECK(jacobian_matrix.isApprox(jacobian_matrix_fd, math::sqrt(eps_fd)));
  }
}

BOOST_AUTO_TEST_SUITE_END()
