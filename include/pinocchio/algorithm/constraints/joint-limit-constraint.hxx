//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__

#include "pinocchio/multibody/joint/joint-basic-visitors.hpp"
#include <iostream>

namespace pinocchio
{
  //  template<typename Scalar, int Options>
  //  template<template<typename, int> class JointCollectionTpl>
  //  int JointLimitConstraintModelTpl<Scalar, Options>::check_activable_joints(
  //    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
  //    const JointIndexVector & activable_joints)
  //  {
  //    for (const JointIndex joint_id : activable_joints)
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
    const JointIndexVector & activable_joints,
    const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
    const Eigen::MatrixBase<VectorUpperConfiguration> & ub)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointModel JointModel;

    PINOCCHIO_CHECK_ARGUMENT_SIZE(lb.size(), model.nq);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(ub.size(), model.nq);

    // Check validity of activable_joints input
    for (const JointIndex joint_id : activable_joints)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        joint_id < model.joints.size(),
        "joint_id is larger than the total number of joints contained in the model.");
    }
    //    PINOCCHIO_CHECK_INPUT_ARGUMENT(
    //      check_activable_joints(model, activable_joints) == -1,
    //      "One of the joint is not supported by JointLimitConstraintModelTpl.")

    // Collect all potentially active bounds
    BooleanVector is_lower_bound_constraint_active = BooleanVector::Constant(model.nq, false),
                  is_upper_bound_constraint_active = BooleanVector::Constant(model.nq, false);

    activable_configuration_components.reserve(size_t(model.nq));
    activable_lower_bound_constraints.reserve(size_t(model.nq));
    activable_lower_bound_constraints_tangent.reserve(size_t(model.nv));
    activable_upper_bound_constraints.reserve(size_t(model.nq));
    activable_upper_bound_constraints_tangent.reserve(size_t(model.nv));
    for (const JointIndex joint_id : activable_joints)
    {
      const JointModel & jmodel = model.joints[joint_id];

      const int idx_q = jmodel.idx_q();
      const int idx_v = jmodel.idx_v();
      const int nq = jmodel.nq();

      const auto has_configuration_limit = jmodel.hasConfigurationLimit();
      for (int k = 0; k < nq; ++k)
      {
        const int q_index = idx_q + k;
        if (!has_configuration_limit[size_t(k)])
          continue;

        const int v_index = idx_v + k;
        if (!(lb[q_index] == -std::numeric_limits<Scalar>::max()
              || lb[q_index] == -std::numeric_limits<Scalar>::infinity()))
        {
          is_lower_bound_constraint_active[q_index] = true;
          activable_lower_bound_constraints.push_back(q_index);
          activable_lower_bound_constraints_tangent.push_back(v_index);
        }
        if (!(ub[q_index] == +std::numeric_limits<Scalar>::max()
              || ub[q_index] == +std::numeric_limits<Scalar>::infinity()))
        {
          is_upper_bound_constraint_active[q_index] = true;
          activable_upper_bound_constraints.push_back(q_index);
          activable_upper_bound_constraints_tangent.push_back(v_index);
        }

        if (is_lower_bound_constraint_active[q_index] || is_upper_bound_constraint_active[q_index])
        {
          activable_configuration_components.push_back(q_index);
        }
      }
    }

    // Fill lower and upper bounds for active components of the configuration vector.
    {
      // TODO: this code should be removed ? activable_configuration_limits is not used anymore
      VectorXs active_lower_bound_limit =
                 VectorXs::Zero(Eigen::DenseIndex(activable_configuration_components.size())),
               active_upper_bound_limit =
                 VectorXs::Zero(Eigen::DenseIndex(activable_configuration_components.size()));
      Eigen::Index row_id = 0;
      for (const auto active_configuration_index : activable_configuration_components)
      {
        active_lower_bound_limit[row_id] = lb[active_configuration_index];
        active_upper_bound_limit[row_id] = ub[active_configuration_index];
        row_id++;
      }

      activable_configuration_limits = BoxSet(active_lower_bound_limit, active_upper_bound_limit);
    }

    // Fill constraint sparsity pattern
    VectofOfEigenIndexVector row_activable_indexes_upper;
    VectofOfEigenIndexVector & row_activable_indexes_lower = row_activable_indexes;

    EigenIndexVector extended_support;
    extended_support.reserve(size_t(model.nv));
    for (const JointIndex joint_id : activable_joints)
    {
      const JointModel & jmodel = model.joints[joint_id];
      const auto & jsupport = model.supports[joint_id];

      const int idx_q = jmodel.idx_q();

      const int nv = jmodel.nv();
      const int idx_v = jmodel.idx_v();

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
          row_activable_indexes_lower.push_back(extended_support);
        if (is_upper_bound_constraint_active[row_id_q])
          row_activable_indexes_upper.push_back(extended_support);
      }
    }

    // append row_activable_indexes_upper to row_activable_indexes_lower
    row_activable_indexes_lower.insert(
      row_activable_indexes_lower.end(), row_activable_indexes_upper.begin(),
      row_activable_indexes_upper.end());

    const size_t total_size =
      activable_lower_bound_constraints.size() + activable_upper_bound_constraints.size();
    assert(row_activable_indexes.size() == total_size);

    // Fill row_activable_sparsity_pattern from row_activable_indexes content
    row_activable_sparsity_pattern.resize(total_size, BooleanVector::Zero(model.nv));
    for (size_t row_id = 0; row_id < total_size; ++row_id)
    {
      auto & sparsity_pattern = row_activable_sparsity_pattern[row_id];
      const auto & extended_support = row_activable_indexes[row_id];

      for (const auto val : extended_support)
        sparsity_pattern[val] = true;
    }

    m_compliance = ComplianceVectorType::Zero(size());
    m_baumgarte_parameters = BaumgarteCorrectorParameters();

    // Allocate the maximum size for the dynamic limits
    active_lower_bound_constraints.reserve(this->getActivableLowerBoundConstraints().size());
    active_lower_bound_constraints_tangent.reserve(
      this->getActivableLowerBoundConstraintsTangent().size());
    active_upper_bound_constraints.reserve(this->getActivableUpperBoundConstraints().size());
    active_upper_bound_constraints_tangent.reserve(
      this->getActivableUpperBoundConstraintsTangent().size());
    active_set_indexes.reserve(
      active_upper_bound_constraints_tangent.capacity()
      + active_lower_bound_constraints_tangent.capacity());
    active_compliance_storage.resize(0);
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl>
  void JointLimitConstraintModelTpl<Scalar, Options>::resize(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    ConstraintData & cdata)
  {
    // Compute notably the constraint constraint_residual
    // This allows to compute which limits are active in the current configuration (data.q_in) which
    // corresponds to the current active set.
    // TODO(jcarpent): change model.lowerPositionLimit[q_index] and
    // model.upperPositionLimit[q_index] for internal limit values stored in the constraint model.
    auto & activable_constraint_residual = cdata.activable_constraint_residual;
    active_lower_bound_constraints.clear();
    active_upper_bound_constraints.clear();
    active_lower_bound_constraints_tangent.clear();
    active_upper_bound_constraints_tangent.clear();
    active_set_indexes.clear();

    std::size_t row_index = 0;

    // Fill the constraint residual for all activable constraints and detect the active ones.
    for (std::size_t i = 0; i < activable_lower_bound_constraints.size(); i++)
    {
      const auto q_index = activable_lower_bound_constraints[i];
      activable_constraint_residual[int(row_index)] =
        -(data.q_in[q_index] - model.lowerPositionLimit[q_index]);
      assert(model.positionLimitMargin[q_index] >= 0);
      if (activable_constraint_residual[int(row_index)] >= -model.positionLimitMargin[q_index])
      {
        const auto v_index = activable_lower_bound_constraints_tangent[i];
        active_lower_bound_constraints.push_back(q_index);
        active_lower_bound_constraints_tangent.push_back(v_index);
        active_set_indexes.push_back(row_index);
      }
      row_index++;
    }

    for (std::size_t i = 0; i < activable_upper_bound_constraints.size(); i++)
    {
      const auto q_index = activable_upper_bound_constraints[i];
      activable_constraint_residual[int(row_index)] =
        -(data.q_in[q_index] - model.upperPositionLimit[q_index]);
      assert(model.positionLimitMargin[q_index] >= 0);
      if (activable_constraint_residual[int(row_index)] <= model.positionLimitMargin[q_index])
      {
        const auto v_index = activable_upper_bound_constraints_tangent[i];
        active_upper_bound_constraints.push_back(q_index);
        active_upper_bound_constraints_tangent.push_back(v_index);
        active_set_indexes.push_back(row_index);
      }
      row_index++;
    }

    assert(row_index == (std::size_t)this->size());

    // Resize the constraint residual/compliance storage to the active set size.
    std::size_t active_size = active_set_indexes.size();
    cdata.constraint_residual_storage.resize(int(active_size));

    // Update the active compliance
    active_compliance_storage.resize(int(active_size));
    for (std::size_t active_row_index = 0; active_row_index < active_size; active_row_index++)
    {
      active_compliance[int(active_row_index)] =
        m_compliance[int(active_set_indexes[active_row_index])];
    }

    // Resize the constraint set so it corresponds to the active set.
    m_set.resize(
      Eigen::DenseIndex(active_lower_bound_constraints.size()),
      Eigen::DenseIndex(active_upper_bound_constraints.size()));
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl>
  void JointLimitConstraintModelTpl<Scalar, Options>::calc(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    ConstraintData & cdata) const
  {
    PINOCCHIO_UNUSED_VARIABLE(model);
    PINOCCHIO_UNUSED_VARIABLE(data);
    // Fill the constraint residual for all active constraints.
    std::size_t active_size = std::size_t(activeSize());
    auto & activable_constraint_residual = cdata.activable_constraint_residual;
    auto & constraint_residual = cdata.constraint_residual;

    for (std::size_t active_row_index = 0; active_row_index < active_size; active_row_index++)
    {
      constraint_residual[int(active_row_index)] =
        activable_constraint_residual[int(active_set_indexes[active_row_index])];
    }
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

    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian_matrix.rows(), this->activeSize(),
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
      jacobian_matrix(row_id, col_id) = -Scalar(1);
    }
    for (size_t constraint_id = 0; constraint_id < active_upper_bound_constraints_tangent.size();
         ++constraint_id, ++row_id)
    {
      const auto col_id = active_upper_bound_constraints_tangent[constraint_id];
      jacobian_matrix(row_id, col_id) = -Scalar(1);
    }
  }

  template<typename Scalar, int Options>
  template<
    typename InputMatrix,
    typename OutputMatrix,
    template<typename, int> class JointCollectionTpl,
    AssignmentOperatorType op>
  void JointLimitConstraintModelTpl<Scalar, Options>::jacobianMatrixProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<InputMatrix> & mat,
    const Eigen::MatrixBase<OutputMatrix> & _res,
    AssignmentOperatorTag<op> aot) const
  {
    OutputMatrix & res = _res.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), model.nv);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.cols(), res.cols());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), this->activeSize());
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(cdata);
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    Eigen::DenseIndex row_id = 0;
    for (size_t constraint_id = 0; constraint_id < active_lower_bound_constraints_tangent.size();
         ++constraint_id, ++row_id)
    {
      const auto col_id = active_lower_bound_constraints_tangent[constraint_id];

      if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(row_id) -= -mat.row(col_id);
      else
        res.row(row_id) += -mat.row(col_id);
    }
    for (size_t constraint_id = 0; constraint_id < active_upper_bound_constraints_tangent.size();
         ++constraint_id, ++row_id)
    {
      const auto col_id = active_upper_bound_constraints_tangent[constraint_id];

      if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(row_id) -= -mat.row(col_id);
      else
        res.row(row_id) += -mat.row(col_id);
    }
  }

  template<typename Scalar, int Options>
  template<
    typename InputMatrix,
    typename OutputMatrix,
    template<typename, int> class JointCollectionTpl,
    AssignmentOperatorType op>
  void JointLimitConstraintModelTpl<Scalar, Options>::jacobianTransposeMatrixProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<InputMatrix> & mat,
    const Eigen::MatrixBase<OutputMatrix> & _res,
    AssignmentOperatorTag<op> aot) const
  {
    OutputMatrix & res = _res.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), this->activeSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.cols(), mat.cols());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(cdata);
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    Eigen::DenseIndex row_id = 0;
    for (size_t constraint_id = 0; constraint_id < active_lower_bound_constraints_tangent.size();
         ++constraint_id, ++row_id)
    {
      const auto col_id = active_lower_bound_constraints_tangent[constraint_id];

      if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(col_id) -= -mat.row(row_id);
      else
        res.row(col_id) += -mat.row(row_id);
    }
    for (size_t constraint_id = 0; constraint_id < active_upper_bound_constraints_tangent.size();
         ++constraint_id, ++row_id)
    {
      const auto col_id = active_upper_bound_constraints_tangent[constraint_id];

      if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(col_id) -= -mat.row(row_id);
      else
        res.row(col_id) += -mat.row(row_id);
    }
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__
