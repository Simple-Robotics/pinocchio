//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/constraints/joint-limit-constraint.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/multibody/sample-models.hpp"

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

  const Model::IndexVector empty_active_joint_ids;

  JointLimitConstraintModel constraint(model, empty_active_joint_ids);
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
  const Model::IndexVector active_joint_ids(ee_support.begin() + 1, ee_support.end());
  JointLimitConstraintModel constraint(model, active_joint_ids);

  BOOST_CHECK(constraint.size() == 0);
}

BOOST_AUTO_TEST_CASE(constraint_constructor)
{
  pinocchio::Model model;
  pinocchio::buildModels::manipulator(model);

  const Data data(model);
  const auto & parents_fromRow = data.parents_fromRow;

  //  std::cout << "model:\n" << model << std::endl;

  const std::string ee_name = "wrist2_joint";
  const JointIndex ee_id = model.getJointId(ee_name);

  const Model::IndexVector & ee_support = model.supports[ee_id];
  const Model::IndexVector active_joint_ids(ee_support.begin() + 1, ee_support.end());

  for (const JointIndex joint_id : active_joint_ids)
  {
    const auto & jmodel = model.joints[joint_id];
    std::cout << "joint type: " << jmodel.shortname() << std::endl;
  }
  JointLimitConstraintModel constraint(model, active_joint_ids);

  // Check size
  {
    int total_size = 0;
    for (const JointIndex joint_id : active_joint_ids)
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

  // Check sparsity pattern
  {
    const EigenIndexVector & active_dofs_lower = constraint.getActiveLowerBoundConstraintsTangent();
    const EigenIndexVector & active_dofs_upper = constraint.getActiveUpperBoundConstraintsTangent();
    EigenIndexVector active_dofs(active_dofs_lower);
    active_dofs.insert(active_dofs.end(), active_dofs_upper.begin(), active_dofs_upper.end());

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

  // Check projection on force sets
  {
    const EigenIndexVector & active_dofs_lower = constraint.getActiveLowerBoundConstraintsTangent();
    const EigenIndexVector & active_dofs_upper = constraint.getActiveUpperBoundConstraintsTangent();

    const Eigen::DenseIndex nb_lower_active_dofs = Eigen::DenseIndex(active_dofs_lower.size());
    const Eigen::DenseIndex nb_upper_active_dofs = Eigen::DenseIndex(active_dofs_upper.size());

    const Eigen::DenseIndex total_active_dofs = nb_lower_active_dofs + nb_upper_active_dofs;

    const int num_projections = int(1e6);
    for (int k = 0; k < num_projections; ++k)
    {
      const Eigen::VectorXd f = Eigen::VectorXd::Random(total_active_dofs);
      const Eigen::VectorXd f_proj = constraint.set().project(f);

      BOOST_CHECK((f_proj.head(nb_lower_active_dofs).array() >= 0).all());
      BOOST_CHECK((f_proj.tail(nb_lower_active_dofs).array() <= 0).all());
    }
  }
}

BOOST_AUTO_TEST_CASE(constraint_jacobian)
{
  pinocchio::Model model;
  pinocchio::buildModels::manipulator(model);

  Data data(model);

  const std::string ee_name = "wrist2_joint";
  const JointIndex ee_id = model.getJointId(ee_name);

  const Model::IndexVector & ee_support = model.supports[ee_id];
  const Model::IndexVector active_joint_ids(ee_support.begin() + 1, ee_support.end());

  JointLimitConstraintModel constraint_model(model, active_joint_ids);
  JointLimitConstraintData constraint_data(constraint_model);

  Eigen::MatrixXd jacobian_matrix(constraint_model.size(), model.nv);

  // Check against finite differences on the drift of the constraint
  const double eps_fd = 1e-8;
  for (int i = 0; i < 1e4; ++i)
  {
    const Eigen::VectorXd q0 = randomConfiguration(model);
    data.q_in = q0;
    constraint_model.calc(model, data, constraint_data);
    constraint_model.jacobian(model, data, constraint_data, jacobian_matrix);
    Data data_fd(model);
    JointLimitConstraintData constraint_data_fd(constraint_model);
    Eigen::MatrixXd jacobian_matrix_fd(constraint_model.size(), model.nv);
    for (Eigen::DenseIndex k = 0; k < model.nv; ++k)
    {
      Eigen::VectorXd v_eps = Eigen::VectorXd::Zero(model.nv);
      v_eps[k] = eps_fd;
      const Eigen::VectorXd q_plus = integrate(model, q0, v_eps);
      data_fd.q_in = q_plus;
      constraint_model.calc(model, data_fd, constraint_data_fd);

      jacobian_matrix_fd.col(k) =
        (constraint_data_fd.constraint_residual - constraint_data.constraint_residual) / eps_fd;
    }

    BOOST_CHECK(jacobian_matrix.isApprox(jacobian_matrix_fd, math::sqrt(eps_fd)));
  }
}

BOOST_AUTO_TEST_SUITE_END()
