//
// Copyright (c) 2022 INRIA
//

#ifndef __pinocchio_algorithm_aba_txx__
#define __pinocchio_algorithm_aba_txx__

namespace pinocchio
{
  namespace impl
  {
    extern template PINOCCHIO_DLLAPI const context::VectorXs & aba<
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
      const Eigen::MatrixBase<Eigen::Ref<const context::VectorXs>> &);

    extern template PINOCCHIO_DLLAPI const context::VectorXs & aba<
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
      const container::aligned_vector<ForceTpl<context::Scalar, context::Options>> &);

    namespace minimal
    {
      extern template PINOCCHIO_DLLAPI const context::VectorXs & aba<
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
        const Eigen::MatrixBase<Eigen::Ref<const context::VectorXs>> &);

      extern template PINOCCHIO_DLLAPI const context::VectorXs & aba<
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
        const container::aligned_vector<ForceTpl<context::Scalar, context::Options>> &);
    } // namespace minimal

    extern template PINOCCHIO_DLLAPI const context::RowMatrixXs & computeMinverse<
      context::Scalar,
      context::Options,
      JointCollectionDefaultTpl,
      Eigen::Ref<const context::VectorXs>>(
      const context::Model &,
      context::Data &,
      const Eigen::MatrixBase<Eigen::Ref<const context::VectorXs>> &);
  } // namespace impl
  extern template PINOCCHIO_DLLAPI const context::RowMatrixXs &
  computeMinverse<context::Scalar, context::Options, JointCollectionDefaultTpl>(
    const context::Model &, context::Data &);
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_aba_txx__
