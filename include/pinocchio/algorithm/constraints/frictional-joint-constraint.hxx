//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_frictional_joint_constraint_hxx__
#define __pinocchio_algorithm_constraints_frictional_joint_constraint_hxx__

namespace pinocchio
{
  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl>
  void FrictionalJointConstraintModelTpl<Scalar, Options>::init(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const JointIndexVector & active_joints)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointModel JointModel;
    active_dofs.reserve(size_t(model.nv));
    for (const JointIndex joint_id : active_joints)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        joint_id < model.joints.size(),
        "joint_id is larger than the total number of joints contained in the model.");
      const JointModel & jmodel = model.joints[joint_id];
      const auto & jsupport = model.supports[joint_id];

      const int nv = jmodel.nv();
      const int idx_v = jmodel.idx_v();

      for (int k = 0; k < nv; ++k)
      {
        const int row_id = idx_v + k;
        active_dofs.push_back(row_id);
      }

      EigenIndexVector extended_support;
      extended_support.reserve(size_t(model.nv));
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
        const int row_id = idx_v + k;
        extended_support.push_back(row_id);
        row_active_indexes.push_back(extended_support);
      }
    }

    const size_t total_size = active_dofs.size();
    assert(
      row_active_indexes.size() == total_size && "The two vectors should be of the same size.");

    // Fill row_sparsity_pattern from row_active_indexes content
    row_sparsity_pattern.resize(total_size, BooleanVector::Zero(model.nv));
    for (size_t row_id = 0; row_id < total_size; ++row_id)
    {
      auto & sparsity_pattern = row_sparsity_pattern[row_id];
      const auto & extended_support = row_active_indexes[row_id];

      for (const auto val : extended_support)
        sparsity_pattern[val] = true;
    }
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
  void FrictionalJointConstraintModelTpl<Scalar, Options>::jacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & /*data*/,
    ConstraintData & /*cdata*/,
    const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const
  {
    JacobianMatrix & jacobian_matrix = _jacobian_matrix.const_cast_derived();

    const FrictionalJointConstraintModelTpl & cmodel = *this;
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian_matrix.rows(), cmodel.size(),
      "The input/output Jacobian matrix does not have the right number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian_matrix.cols(), model.nv,
      "The input/output Jacobian matrix does not have the right number of cols.");

    jacobian_matrix.setZero();
    for (size_t row_id = 0; row_id < active_dofs.size(); ++row_id)
    {
      const auto col_id = active_dofs[row_id];
      jacobian_matrix(Eigen::DenseIndex(row_id), col_id) = Scalar(1);
    }
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_frictional_joint_constraint_hxx__
