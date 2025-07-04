//
// Copyright (c) 2022 INRIA
//

#include "pinocchio/algorithm/aba.hpp"

namespace pinocchio
{
  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI const context::VectorXs & aba<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    Eigen::Ref<const context::VectorXs>,
    Eigen::Ref<const context::VectorXs>,
    Eigen::Ref<const context::VectorXs>>(
    const context::Model &,
    context::Data &,
    const Eigen::MatrixBase<Eigen::Ref<const context::VectorXs>> &,
    const Eigen::MatrixBase<Eigen::Ref<const context::VectorXs>> &,
    const Eigen::MatrixBase<Eigen::Ref<const context::VectorXs>> &,
    const Convention);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI const context::VectorXs & aba<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    Eigen::Ref<const context::VectorXs>,
    Eigen::Ref<const context::VectorXs>,
    Eigen::Ref<const context::VectorXs>,
    context::Force,
    Eigen::aligned_allocator<context::Force>>(
    const context::Model &,
    context::Data &,
    const Eigen::MatrixBase<Eigen::Ref<const context::VectorXs>> &,
    const Eigen::MatrixBase<Eigen::Ref<const context::VectorXs>> &,
    const Eigen::MatrixBase<Eigen::Ref<const context::VectorXs>> &,
    const std::vector<context::Force, Eigen::aligned_allocator<context::Force>> &,
    const Convention);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI const context::RowMatrixXs &
  computeMinverse<context::Scalar, context::Options, JointCollectionDefaultTpl>(
    const context::Model &, context::Data &);
} // namespace pinocchio
