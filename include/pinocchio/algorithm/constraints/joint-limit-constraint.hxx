//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__

#include <iostream>
#include <algorithm>

#include "pinocchio/multibody/joint/joint-basic-visitors.hpp"

namespace pinocchio
{

  template<typename Scalar, int Options>
  template<
    template<typename, int> class JointCollectionTpl,
    typename VectorLowerConfiguration,
    typename VectorUpperConfiguration,
    typename VectorMarginConfiguration>
  void JointLimitConstraintModelTpl<Scalar, Options>::init(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const JointIndexVector & _activable_joints,
    const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
    const Eigen::MatrixBase<VectorUpperConfiguration> & ub,
    const Eigen::MatrixBase<VectorMarginConfiguration> & marg)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointModel JointModel;

    PINOCCHIO_CHECK_ARGUMENT_SIZE(lb.size(), model.nq);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(ub.size(), model.nq);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(marg.size(), model.nq);

    // Check validity of _activable_joints input
    for (const JointIndex joint_id : _activable_joints)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        joint_id < model.joints.size(),
        "joint_id is larger than the total number of joints contained in the model.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(joint_id > 0, "joint_id is not valid.");
    }

    //    PINOCCHIO_CHECK_INPUT_ARGUMENT(
    //      check__activable_joints(model, _activable_joints) == -1,
    //      "One of the joint is not supported by JointLimitConstraintModelTpl.")

    // TODO: Should we reserve some activable quantities ?

    // Loop on all q components of activable jointds to identify activable lower and upper
    // constraints, and for each track row_id of related activable joint, idx_q in the configuration
    // and idx_q_reduce in the subpart of q due to activable joints
    VectorOfSize & activable_idx_rows_lower = activable_idx_rows;
    VectorOfSize activable_idx_rows_upper;

    EigenIndexVector & activable_idx_qs_reduce_lower = activable_idx_qs_reduce;
    EigenIndexVector activable_idx_qs_reduce_upper;

    EigenIndexVector & activable_idx_qs_lower = activable_idx_qs;
    EigenIndexVector activable_idx_qs_upper;

    // Prepare the structure to compute sparsity pattern
    EigenIndexVector extended_support;
    extended_support.reserve(size_t(model.nv));

    size_t idx_row = 0;
    nq_reduce = 0;
    for (const JointIndex joint_id : _activable_joints)
    {
      const JointModel & jmodel = model.joints[joint_id];

      const int idx_q = jmodel.idx_q();
      const int idx_v = jmodel.idx_v();
      const int nq = jmodel.nq();
      const int nv = jmodel.nv();
      const auto has_configuration_limit = jmodel.hasConfigurationLimit();

      bool is_joint_really_active = false;
      for (int j_qi = 0; j_qi < nq; ++j_qi)
      {
        if (!has_configuration_limit[size_t(j_qi)])
          continue;

        const int q_index = idx_q + j_qi;
        const int q_reduce_index = nq_reduce + j_qi;

        if (!(lb[q_index] == -std::numeric_limits<Scalar>::max()
              || lb[q_index] == -std::numeric_limits<Scalar>::infinity()))
        {
          activable_idx_rows_lower.push_back(idx_row);
          activable_idx_qs_lower.push_back(q_index);
          activable_idx_qs_reduce_lower.push_back(q_reduce_index);
          is_joint_really_active = true;
        }
        if (!(ub[q_index] == +std::numeric_limits<Scalar>::max()
              || ub[q_index] == +std::numeric_limits<Scalar>::infinity()))
        {
          activable_idx_rows_upper.push_back(idx_row);
          activable_idx_qs_upper.push_back(q_index);
          activable_idx_qs_reduce_upper.push_back(q_reduce_index);
          is_joint_really_active = true;
        }
      }

      // At least one lower or upper constraint for a component of the joint is active so update the
      // quantity
      if (is_joint_really_active)
      {
        activable_joints.push_back(joint_id);
        idx_row += 1;
        nq_reduce += nq;

        // Compute the sparsity pattern of the joint
        const auto & jsupport = model.supports[joint_id];
        extended_support.clear();
        for (size_t i = 1; i < jsupport.size() - 1; ++i)
        {
          const JointIndex jsupport_id = jsupport[i];
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
          const int extended_row_id = idx_v + k;
          extended_support.push_back(extended_row_id);
        }
        row_indexes.push_back(extended_support);
      }
    }

    // Recover sizes of constraints
    lower_activable_size = static_cast<int>(activable_idx_rows_lower.size());
    int upper_activable_size = static_cast<int>(activable_idx_rows_upper.size());
    int activable_size = lower_activable_size + upper_activable_size;
    PINOCCHIO_UNUSED_VARIABLE(activable_size);

    // Recompose one vectors for all constraint with convention lower | upper
    activable_idx_rows.insert(
      activable_idx_rows.end(), activable_idx_rows_upper.begin(), activable_idx_rows_upper.end());
    activable_idx_qs_reduce.insert(
      activable_idx_qs_reduce.end(), activable_idx_qs_reduce_upper.begin(),
      activable_idx_qs_reduce_upper.end());
    activable_idx_qs.insert(
      activable_idx_qs.end(), activable_idx_qs_upper.begin(), activable_idx_qs_upper.end());
    assert(size() == activable_size);

    // Fill bound limit and margin for lower and upper constraint
    // Another strategy could be to query the model again but it is not coherent with the existing
    // constructors.
    bound_position_limit = VectorXs::Zero(Eigen::DenseIndex(size()));
    bound_position_margin = VectorXs::Zero(Eigen::DenseIndex(size()));
    Eigen::DenseIndex bound_row_id = 0;
    for (std::size_t i = 0; i < static_cast<std::size_t>(lowerSize()); i++)
    {
      const auto activable_idx_q = activable_idx_qs[i];
      bound_position_limit[bound_row_id] = lb[activable_idx_q];
      assert(marg[activable_idx_q] >= 0);
      bound_position_margin[bound_row_id] = marg[activable_idx_q];
      bound_row_id++;
    }
    for (std::size_t i = static_cast<std::size_t>(lowerSize());
         i < static_cast<std::size_t>(size()); i++)
    {
      const auto activable_idx_q = activable_idx_qs[i];
      bound_position_limit[bound_row_id] = ub[activable_idx_q];
      assert(marg[activable_idx_q] >= 0);
      bound_position_margin[bound_row_id] = marg[activable_idx_q];
      bound_row_id++;
    }
    assert(bound_row_id == static_cast<Eigen::DenseIndex>(size()));

    // Get nvs and idx_vs of all actibale joints to compute nv_max_atom
    // and activable_nvs, activable_idx_vs
    std::vector<int> reduce_nvs, reduce_idx_vs;
    reduce_nvs.reserve(static_cast<std::size_t>(nq_reduce));
    reduce_idx_vs.reserve(static_cast<std::size_t>(nq_reduce));

    pinocchio::indexvInfo(model, activable_joints, reduce_nvs, reduce_idx_vs);
    assert(nq_reduce == static_cast<int>(reduce_nvs.size()));
    assert(nq_reduce == static_cast<int>(reduce_idx_vs.size()));
    auto nv_max_atom_iter = std::max_element(reduce_nvs.begin(), reduce_nvs.end());
    nv_max_atom = nv_max_atom_iter != reduce_nvs.end() ? *nv_max_atom_iter : 1;
    assert(nv_max_atom <= MAX_JOINT_NV);

    std::size_t r_size = static_cast<std::size_t>(size());
    activable_nvs.reserve(r_size);
    activable_idx_vs.reserve(r_size);
    for (const auto activable_idx_q_reduce : activable_idx_qs_reduce)
    {
      std::size_t idx_query = static_cast<std::size_t>(activable_idx_q_reduce);
      activable_nvs.push_back(static_cast<Eigen::DenseIndex>(reduce_nvs[idx_query]));
      activable_idx_vs.push_back(static_cast<Eigen::DenseIndex>(reduce_idx_vs[idx_query]));
    }
    assert(r_size == activable_nvs.size());
    assert(r_size == activable_idx_vs.size());

    // Fill row_activable_sparsity_pattern from row_activable_indexes content
    row_sparsity_pattern.resize(row_indexes.size(), BooleanVector::Zero(model.nv));
    for (size_t joint_id = 0; joint_id < row_indexes.size(); ++joint_id)
    {
      auto & sparsity_pattern = row_sparsity_pattern[joint_id];
      const auto & extended_support = row_indexes[joint_id];
      for (const auto val : extended_support)
        sparsity_pattern[val] = true;
    }

    m_compliance = ComplianceVectorType::Zero(size());
    m_baumgarte_parameters = BaumgarteCorrectorParameters();

    // Allocate the maximum size for the dynamic quantity
    lower_active_size = 0;
    active_set_indexes.reserve(r_size);
    active_idx_rows.reserve(r_size);
    active_idx_qs_reduce.reserve(r_size);
    active_nvs.reserve(r_size);
    active_idx_vs.reserve(r_size);
    // active_compliance_storage.resize(size()) would allocate double size...
    // active_compliance_storage.resize(0);  and we resize again
    active_compliance_storage.reserve(size());
    assert(activeSize() == lowerActiveSize() == upperActiveSize() == 0);
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl>
  void JointLimitConstraintModelTpl<Scalar, Options>::resize_impl(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & /* model */,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    ConstraintData & cdata)
  {
    // Compute notably the constraint constraint_residual
    // This allows to compute which limits are active in the current configuration (data.q_in) which
    // corresponds to the current active set.
    auto & activable_constraint_residual = cdata.activable_constraint_residual;

    active_set_indexes.clear();
    active_idx_rows.clear();
    active_idx_qs_reduce.clear();
    active_nvs.clear();
    active_idx_vs.clear();
    lower_active_size = 0;

    // Fill the constraint residual for all activable constraints and detect the active ones.
    // The convention is lower | upper, and negative | positive constraint so:
    // q_l <= q + TMv <= q_up
    // -TMv + (q_l - q) <= 0
    // -TMv + (q_u - q) >= 0

    // We compute all active quanties
    // active_[idx_rows|idx_qs_reduce|nvs|idx_vs] are store
    // but they are not necessary as they are recoverable from active_set_indexes
    // However it imply double referencing in all jacobian methods
    // And for one call of resize/calc, their can be multiple call to jacobian methods !

    // Lower
    for (std::size_t i = 0; i < static_cast<std::size_t>(lowerSize()); i++)
    {
      const Eigen::DenseIndex ie = static_cast<Eigen::DenseIndex>(i);
      const Eigen::DenseIndex idx_q = activable_idx_qs[i];
      activable_constraint_residual[ie] = bound_position_limit[ie] - data.q_in[idx_q];
      if (activable_constraint_residual[ie] >= -bound_position_margin[ie])
      {
        active_set_indexes.push_back(i);
        active_idx_rows.push_back(activable_idx_rows[i]);
        active_idx_qs_reduce.push_back(activable_idx_qs_reduce[i]);
        active_nvs.push_back(activable_nvs[i]);
        active_idx_vs.push_back(activable_idx_vs[i]);
        lower_active_size += 1;
      }
    }
    // Upper
    for (std::size_t i = static_cast<std::size_t>(lowerSize());
         i < static_cast<std::size_t>(size()); i++)
    {
      const Eigen::DenseIndex ie = static_cast<Eigen::DenseIndex>(i);
      const Eigen::DenseIndex idx_q = activable_idx_qs[i];
      activable_constraint_residual[ie] = bound_position_limit[ie] - data.q_in[idx_q];
      if (activable_constraint_residual[ie] <= bound_position_margin[ie])
      {
        active_set_indexes.push_back(i);
        active_idx_rows.push_back(activable_idx_rows[i]);
        active_idx_qs_reduce.push_back(activable_idx_qs_reduce[i]);
        active_nvs.push_back(activable_nvs[i]);
        active_idx_vs.push_back(activable_idx_vs[i]);
      }
    }

    // Resize the constraint residual/compliance storage to the active set size.
    const int active_size = activeSize();
    cdata.constraint_residual_storage.resize(active_size);

    // Update the active compliance
    active_compliance_storage.resize(active_size);
    for (int active_row_index = 0; active_row_index < active_size; active_row_index++)
    {
      active_compliance[active_row_index] = m_compliance[static_cast<int>(
        active_set_indexes[static_cast<std::size_t>(active_row_index)])];
    }

    // Resize the constraint set so it corresponds to the active set.
    m_set.resize(Eigen::DenseIndex(lowerActiveSize()), Eigen::DenseIndex(upperActiveSize()));
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

    const std::size_t active_size = static_cast<std::size_t>(this->activeSize());
    auto & activable_constraint_residual = cdata.activable_constraint_residual;
    auto & constraint_residual = cdata.constraint_residual;

    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      constraint_residual.size(), this->activeSize(),
      "The active constraint_residual size in constraint data is different from the constraint "
      "model active size. You should probably use cmodel.resize(model, data, cdata) first.");

    // Fill the constraint residual for all active constraints.
    for (std::size_t active_row_index = 0; active_row_index < active_size; active_row_index++)
    {
      constraint_residual[int(active_row_index)] =
        activable_constraint_residual[int(active_set_indexes[active_row_index])];
    }

    // Fill the compact tangent map
    pinocchio::compactTangentMap(model, activable_joints, data.q_in, cdata.compact_tangent_map);
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
  void JointLimitConstraintModelTpl<Scalar, Options>::jacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & /*data*/,
    ConstraintData & cdata,
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
    const CompactTangentMap_t & TMc = cdata.compact_tangent_map;
    Eigen::DenseIndex row_id = 0;
    for (size_t constraint_id = 0; constraint_id < static_cast<std::size_t>(activeSize());
         ++constraint_id, ++row_id)
    {
      jacobian_matrix.block(row_id, active_idx_vs[constraint_id], 1, active_nvs[constraint_id]) =
        -TMc.block(active_idx_qs_reduce[constraint_id], 0, 1, active_nvs[constraint_id]);
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
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    const CompactTangentMap_t & TMc = cdata.compact_tangent_map;
    Eigen::DenseIndex row_id = 0;
    for (size_t constraint_id = 0; constraint_id < static_cast<std::size_t>(activeSize());
         ++constraint_id, ++row_id)
    {
      if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(row_id) -=
          -TMc.block(active_idx_qs_reduce[constraint_id], 0, 1, active_nvs[constraint_id])
          * mat.middleRows(active_idx_vs[constraint_id], active_nvs[constraint_id]);
      else
        res.row(row_id) +=
          -TMc.block(active_idx_qs_reduce[constraint_id], 0, 1, active_nvs[constraint_id])
          * mat.middleRows(active_idx_vs[constraint_id], active_nvs[constraint_id]);
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
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    const CompactTangentMap_t & TMc = cdata.compact_tangent_map;
    Eigen::DenseIndex row_id = 0;
    for (size_t constraint_id = 0; constraint_id < static_cast<std::size_t>(activeSize());
         ++constraint_id, ++row_id)
    {
      if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.middleRows(active_idx_vs[constraint_id], active_nvs[constraint_id]) -=
          -TMc.block(active_idx_qs_reduce[constraint_id], 0, 1, active_nvs[constraint_id])
             .transpose()
          * mat.row(row_id);
      else
        res.middleRows(active_idx_vs[constraint_id], active_nvs[constraint_id]) +=
          -TMc.block(active_idx_qs_reduce[constraint_id], 0, 1, active_nvs[constraint_id])
             .transpose()
          * mat.row(row_id);
    }
  }

  template<typename Scalar, int Options>
  template<
    template<typename, int> class JointCollectionTpl,
    typename VectorNLike,
    ReferenceFrame rf>
  void JointLimitConstraintModelTpl<Scalar, Options>::appendCouplingConstraintInertias(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
    const ReferenceFrameTag<rf> reference_frame) const
  {
    PINOCCHIO_UNUSED_VARIABLE(model);
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(cdata);
    PINOCCHIO_UNUSED_VARIABLE(diagonal_constraint_inertia);
    PINOCCHIO_UNUSED_VARIABLE(reference_frame);
    // TODO(jcarpent)
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__
