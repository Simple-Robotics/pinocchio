//
// Copyright (c) 2022-2024 INRIA
//

#include "pinocchio/spatial/fwd.hpp"

#ifndef PINOCCHIO_SKIP_ALGORITHM_CONTACT_JACOBIAN

  #include "pinocchio/algorithm/contact-jacobian.hpp"

namespace pinocchio
{
  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void getConstraintJacobian<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::RigidConstraintModel,
    context::RigidConstraintData,
    context::MatrixXs>(
    const context::Model &,
    const context::Data &,
    const ConstraintModelBase<context::RigidConstraintModel> &,
    ConstraintDataBase<context::RigidConstraintData> &,
    const Eigen::MatrixBase<context::MatrixXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void getConstraintsJacobian<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::MatrixXs,
    context::RigidConstraintModel,
    typename context::RigidConstraintModelVector::allocator_type,
    context::RigidConstraintData,
    typename context::RigidConstraintDataVector::allocator_type>(
    const context::Model &,
    const context::Data &,
    const context::RigidConstraintModelVector &,
    context::RigidConstraintDataVector &,
    const Eigen::MatrixBase<context::MatrixXs> &);

} // namespace pinocchio

#endif // PINOCCHIO_SKIP_ALGORITHM_CONTACT_JACOBIAN
