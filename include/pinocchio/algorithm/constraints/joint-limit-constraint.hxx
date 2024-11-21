//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__

#include "pinocchio/multibody/joint/joint-basic-visitors.hpp"

namespace pinocchio
{
  //  template<typename Scalar, int Options>
  //  template<template<typename, int> class JointCollectionTpl>
  //  int JointLimitConstraintModelTpl<Scalar, Options>::check_active_joints(
  //    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
  //    const JointIndexVector & active_joints)
  //  {
  //    for (const JointIndex joint_id : active_joints)
  //    {
  //      const JointModel & jmodel = model.joints[joint_id];
  //
  //      if (!check_joint_type_within_sequence<ValidJointTypes>(jmodel))
  //        return int(joint_id);
  //    }
  //
  //    return -1;
  //  }

  template<typename Scalar, int Options>
  template<
    template<typename, int> class JointCollectionTpl,
    typename VectorLowerConfiguration,
    typename VectorUpperConfiguration>
  void JointLimitConstraintModelTpl<Scalar, Options>::init(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const JointIndexVector & active_joints,
    const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
    const Eigen::MatrixBase<VectorUpperConfiguration> & ub)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointModel JointModel;

    PINOCCHIO_CHECK_ARGUMENT_SIZE(lb.size(), model.nq);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(ub.size(), model.nq);

    // Check validity of active_joints input
    for (const JointIndex joint_id : active_joints)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        joint_id < model.joints.size(),
        "joint_id is larger than the total number of joints contained in the model.");
    }
    //    PINOCCHIO_CHECK_INPUT_ARGUMENT(
    //      check_active_joints(model, active_joints) == -1,
    //      "One of the joint is not supported by JointLimitConstraintModelTpl.")

    // Collect all potentially active bounds
    BooleanVector is_lower_bound_constraint_active = BooleanVector::Constant(model.nq, false),
                  is_upper_bound_constraint_active = BooleanVector::Constant(model.nq, false);

    active_configuration_components.reserve(size_t(model.nq));
    active_lower_bound_constraints.reserve(size_t(model.nq));
    active_lower_bound_constraints_tangent.reserve(size_t(model.nv));
    active_upper_bound_constraints.reserve(size_t(model.nq));
    active_upper_bound_constraints_tangent.reserve(size_t(model.nv));
    for (const JointIndex joint_id : active_joints)
    {
      const JointModel & jmodel = model.joints[joint_id];

      const int idx_q = jmodel.idx_q();
      const int idx_v = jmodel.idx_v();
      const int nq = jmodel.nq();

      assert(nq == jmodel.nv() && "joint nv and nq dimensions should be equal.");
      const auto has_configuration_limit = jmodel.hasConfigurationLimit();
      for (int k = 0; k < nq; ++k)
      {
        const int q_index = idx_q + k;
        if (!has_configuration_limit[size_t(k)])
          continue;

        const int v_index = idx_v + k;
        if (lb[q_index] != -std::numeric_limits<Scalar>::max())
        {
          is_lower_bound_constraint_active[q_index] = true;
          active_lower_bound_constraints.push_back(q_index);
          active_lower_bound_constraints_tangent.push_back(v_index);
        }
        if (ub[q_index] != +std::numeric_limits<Scalar>::max())
        {
          is_upper_bound_constraint_active[q_index] = true;
          active_upper_bound_constraints.push_back(q_index);
          active_upper_bound_constraints_tangent.push_back(v_index);
        }

        if (is_lower_bound_constraint_active[q_index] || is_upper_bound_constraint_active[q_index])
        {
          active_configuration_components.push_back(q_index);
        }
      }
    }

    // Fill lower and upper bounds for active components of the configuration vector.
    {
      VectorXs active_lower_bound_limit =
                 VectorXs::Zero(Eigen::DenseIndex(active_configuration_components.size())),
               active_upper_bound_limit =
                 VectorXs::Zero(Eigen::DenseIndex(active_configuration_components.size()));
      Eigen::Index row_id = 0;
      for (const auto active_configuration_index : active_configuration_components)
      {
        active_lower_bound_limit[row_id] = lb[active_configuration_index];
        active_upper_bound_limit[row_id] = ub[active_configuration_index];
        row_id++;
      }

      active_configuration_limits = BoxSet(active_lower_bound_limit, active_upper_bound_limit);
    }

    // Fill constraint sparsity pattern
    VectofOfEigenIndexVector row_active_indexes_upper;
    VectofOfEigenIndexVector & row_active_indexes_lower = row_active_indexes;

    EigenIndexVector extended_support;
    extended_support.reserve(size_t(model.nv));
    for (const JointIndex joint_id : active_joints)
    {
      const JointModel & jmodel = model.joints[joint_id];
      const auto & jsupport = model.supports[joint_id];

      const int nq = jmodel.nq();
      const int idx_q = jmodel.idx_q();

      const int nv = jmodel.nv();
      const int idx_v = jmodel.idx_v();

      PINOCCHIO_USED_VARIABLE_FOR_DEBUG_ONLY(nq);
      assert(nq == nv && "joint nv and nq dimensions should be equal.");

      extended_support.clear();
      for (size_t j = 1; j < jsupport.size() - 1; ++j)
      {
        const JointIndex jsupport_id = jsupport[j];
        const JointModel & jsupport = model.joints[jsupport_id];

        const int jsupport_nv = jsupport.nv();
        const int jsupport_idx_v = jsupport.idx_v();

        for (int k = 0; k < jsupport_nv; ++k)
        {
          const int extended_row_id = jsupport_idx_v + k;
          extended_support.push_back(extended_row_id);
        }
      }

      for (int k = 0; k < nv; ++k)
      {
        // TODO(jcarpent): potential issue for mapping row_id_v and row_id_q together for joints
        // with nq != nv.
        const int row_id_v = idx_v + k;
        const int row_id_q = idx_q + k;

        extended_support.push_back(row_id_v);
        if (is_lower_bound_constraint_active[row_id_q])
          row_active_indexes_lower.push_back(extended_support);
        if (is_upper_bound_constraint_active[row_id_q])
          row_active_indexes_upper.push_back(extended_support);
      }
    }

    // append row_active_indexes_upper to row_active_indexes_lower
    row_active_indexes_lower.insert(
      row_active_indexes_lower.end(), row_active_indexes_upper.begin(),
      row_active_indexes_upper.end());

    const size_t total_size =
      active_lower_bound_constraints.size() + active_upper_bound_constraints.size();
    assert(row_active_indexes.size() == total_size);

    // Fill force limits
    m_set.resize(Eigen::DenseIndex(total_size));
    auto & set_lb = m_set.lb();
    set_lb.setZero();
    auto & set_ub = m_set.ub();
    set_ub.setZero();

    // Compare to the reference document, we need to reverse the box bounds, as -force \in BoxSet
    set_ub.head(Eigen::DenseIndex(active_lower_bound_constraints.size()))
      .fill(+std::numeric_limits<Scalar>::max());

    set_lb.tail(Eigen::DenseIndex(active_upper_bound_constraints.size()))
      .fill(-std::numeric_limits<Scalar>::max());

    // Fill row_sparsity_pattern from row_active_indexes content
    row_sparsity_pattern.resize(total_size, BooleanVector::Zero(model.nv));
    for (size_t row_id = 0; row_id < total_size; ++row_id)
    {
      auto & sparsity_pattern = row_sparsity_pattern[row_id];
      const auto & extended_support = row_active_indexes[row_id];

      for (const auto val : extended_support)
        sparsity_pattern[val] = true;
    }

    m_compliance = ComplianceVectorType::Zero(size());
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl>
  void JointLimitConstraintModelTpl<Scalar, Options>::calc(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    ConstraintData & cdata) const
  {
    // Compute notably the constraint constraint_residual
    // TODO(jcarpent): change model.lowerPositionLimit[q_index] and
    // model.upperPositionLimit[q_index] for internal limit values stored in the constraint model.
    auto & constraint_residual = cdata.constraint_residual;
    Eigen::DenseIndex row_index = 0;

    for (const auto q_index : active_lower_bound_constraints)
    {
      constraint_residual[row_index++] = data.q_in[q_index] - model.lowerPositionLimit[q_index];
    }

    for (const auto q_index : active_upper_bound_constraints)
    {
      constraint_residual[row_index++] = data.q_in[q_index] - model.upperPositionLimit[q_index];
    }

    assert(row_index == this->size());
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
  void JointLimitConstraintModelTpl<Scalar, Options>::jacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & /*data*/,
    ConstraintData & /*cdata*/,
    const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const
  {
    JacobianMatrix & jacobian_matrix = _jacobian_matrix.const_cast_derived();

    const JointLimitConstraintModelTpl & cmodel = *this;
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian_matrix.rows(), cmodel.size(),
      "The input/output Jacobian matrix does not have the right number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian_matrix.cols(), model.nv,
      "The input/output Jacobian matrix does not have the right number of cols.");

    jacobian_matrix.setZero();
    Eigen::DenseIndex row_id = 0;
    for (size_t constraint_id = 0; constraint_id < active_lower_bound_constraints_tangent.size();
         ++constraint_id, ++row_id)
    {
      const auto col_id = active_lower_bound_constraints_tangent[constraint_id];
      jacobian_matrix(row_id, col_id) = Scalar(1);
    }
    for (size_t constraint_id = 0; constraint_id < active_upper_bound_constraints_tangent.size();
         ++constraint_id, ++row_id)
    {
      const auto col_id = active_upper_bound_constraints_tangent[constraint_id];
      jacobian_matrix(row_id, col_id) = Scalar(1);
    }
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__
