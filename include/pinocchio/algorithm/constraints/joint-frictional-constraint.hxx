//
// Copyright (c) 2024-2025 INRIA
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
    active_dofs.reserve(size_t(model.nv));
    for (const JointIndex joint_id : active_joints)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        joint_id < model.joints.size(),
        "joint_id is larger than the total number of joints contained in the model.");
      const auto & jsupport = model.supports[joint_id];

      const auto nv = model.nvs[joint_id];
      const auto idx_v = model.idx_vs[joint_id];

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

        const int jsupport_nv = model.nvs[jsupport_id];
        const int jsupport_idx_v = model.idx_vs[jsupport_id];

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

    m_compliance = ComplianceVectorType::Zero(activeSize());
    m_set = ConstraintSet(activeSize());
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
      jacobian_matrix.rows(), cmodel.activeSize(),
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

  template<typename Scalar, int Options>
  template<
    typename InputMatrix,
    typename OutputMatrix,
    template<typename, int> class JointCollectionTpl,
    AssignmentOperatorType op>
  void FrictionalJointConstraintModelTpl<Scalar, Options>::jacobianMatrixProduct(
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
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), activeSize());
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(cdata);
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    for (size_t row_id = 0; row_id < active_dofs.size(); ++row_id)
    {
      const auto col_id = active_dofs[row_id];

      if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(Eigen::DenseIndex(row_id)) -= mat.row(col_id);
      else
        res.row(Eigen::DenseIndex(row_id)) += mat.row(col_id);
    }
  }

  template<typename Scalar, int Options>
  template<
    typename InputMatrix,
    typename OutputMatrix,
    template<typename, int> class JointCollectionTpl,
    AssignmentOperatorType op>
  void FrictionalJointConstraintModelTpl<Scalar, Options>::jacobianTransposeMatrixProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<InputMatrix> & mat,
    const Eigen::MatrixBase<OutputMatrix> & _res,
    AssignmentOperatorTag<op> aot) const
  {
    OutputMatrix & res = _res.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), activeSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.cols(), mat.cols());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(cdata);
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    for (size_t row_id = 0; row_id < active_dofs.size(); ++row_id)
    {
      const auto col_id = active_dofs[row_id];

      if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(col_id) -= mat.row(Eigen::DenseIndex(row_id));
      else
        res.row(col_id) += mat.row(Eigen::DenseIndex(row_id));
    }
  }

  template<typename Scalar, int Options>
  template<
    template<typename, int> class JointCollectionTpl,
    typename VectorNLike,
    ReferenceFrame rf>
  void FrictionalJointConstraintModelTpl<Scalar, Options>::appendCouplingConstraintInertias(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
    const ReferenceFrameTag<rf> reference_frame) const
  {
    PINOCCHIO_UNUSED_VARIABLE(cdata);
    PINOCCHIO_UNUSED_VARIABLE(reference_frame);

    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      diagonal_constraint_inertia.size(), activeSize(),
      "The diagonal_constraint_inertia is of wrong size.");

    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
    using Matrix6 = typename Data::Inertia::Matrix6;

    Eigen::DenseIndex row_id = 0;
    Matrix6 SI;
    for (const JointIndex joint_id : active_joints)
    {
      const auto joint_nv = model.nvs[joint_id];
      const auto joint_idx_v = model.idx_vs[joint_id];
      const auto joint_diagonal_constraint_inertia =
        diagonal_constraint_inertia.segment(row_id, joint_nv);

      data.joint_apparent_inertia.segment(joint_idx_v, joint_nv) +=
        joint_diagonal_constraint_inertia;

      row_id += joint_nv;
    }
  }

  template<typename Scalar, int Options>
  template<
    template<typename, int> class JointCollectionTpl,
    typename ConstraintForceLike,
    typename JointTorqueLike>
  void FrictionalJointConstraintModelTpl<Scalar, Options>::mapConstraintForceToJointTorques(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<ConstraintForceLike> & constraint_forces,
    const Eigen::MatrixBase<JointTorqueLike> & joint_torques_) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), activeSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_torques_.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(cdata);

    auto & joint_torques = joint_torques_.const_cast_derived();

    for (size_t dof_id = 0; dof_id < active_dofs.size(); ++dof_id)
    {
      const auto row_id = active_dofs[dof_id];

      joint_torques.row(row_id) += constraint_forces.row(Eigen::DenseIndex(dof_id));
    }
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_frictional_joint_constraint_hxx__
