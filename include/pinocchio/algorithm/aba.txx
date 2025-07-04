//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_aba_txx__
#define __pinocchio_algorithm_aba_txx__

namespace pinocchio
{
  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI const context::VectorXs & aba<
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

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI const context::VectorXs & aba<
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

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI const context::RowMatrixXs &
  computeMinverse<context::Scalar, context::Options, JointCollectionDefaultTpl>(
    const context::Model &, context::Data &);
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_aba_txx__
