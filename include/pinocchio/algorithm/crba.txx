//
// Copyright (c) 2022 INRIA
//

#ifndef __pinocchio_algorithm_crba_txx__
#define __pinocchio_algorithm_crba_txx__

namespace pinocchio
{
  namespace impl
  {
    namespace minimal
    {

      extern template PINOCCHIO_DLLAPI const context::MatrixXs & crba<
        context::Scalar,
        context::Options,
        JointCollectionDefaultTpl,
        Eigen::Ref<const context::VectorXs>>(
        const context::Model &,
        context::Data &,
        const Eigen::MatrixBase<Eigen::Ref<const context::VectorXs>> &);

    } // namespace minimal

    extern template PINOCCHIO_DLLAPI const context::MatrixXs & crba<
      context::Scalar,
      context::Options,
      JointCollectionDefaultTpl,
      Eigen::Ref<const context::VectorXs>>(
      const context::Model &,
      context::Data &,
      const Eigen::MatrixBase<Eigen::Ref<const context::VectorXs>> &);
  } // namespace impl
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_crba_txx__
