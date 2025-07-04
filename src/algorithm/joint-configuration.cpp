//
// Copyright (c) 2022 INRIA
//

#include "pinocchio/algorithm/joint-configuration.hpp"

namespace pinocchio
{

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void integrate<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void integrate<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void interpolate<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const context::Scalar &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void interpolate<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const context::Scalar &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void difference<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void difference<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void squaredDistance<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void squaredDistance<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void randomConfiguration<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void randomConfiguration<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void neutral<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs>(const context::Model &, const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void
  neutral<context::Scalar, context::Options, JointCollectionDefaultTpl, context::VectorXs>(
    const context::Model &, const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void dIntegrate<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const ArgumentPosition,
    const AssignmentOperatorType);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void dIntegrate<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const ArgumentPosition);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void dIntegrate<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const ArgumentPosition,
    const AssignmentOperatorType);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void tangentMap<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const AssignmentOperatorType);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void tangentMap<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const AssignmentOperatorType);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void tangentMapProduct<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::MatrixXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const AssignmentOperatorType);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void tangentMapProduct<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::MatrixXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const AssignmentOperatorType);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void tangentMapTransposeProduct<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::MatrixXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const AssignmentOperatorType);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void tangentMapTransposeProduct<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::MatrixXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const AssignmentOperatorType);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void dIntegrateTransport<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::MatrixXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const ArgumentPosition);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void dIntegrateTransport<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::MatrixXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const ArgumentPosition);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void dIntegrateTransport<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const ArgumentPosition);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void dDifference<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const ArgumentPosition);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void dDifference<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &,
    const ArgumentPosition);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::Scalar squaredDistanceSum<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::Scalar squaredDistanceSum<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::Scalar distance<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::Scalar distance<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void normalize<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs>(const context::Model &, const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void
  normalize<context::Scalar, context::Options, JointCollectionDefaultTpl, context::VectorXs>(
    const context::Model &, const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void
  lieGroup<LieGroupMap, context::Scalar, context::Options, JointCollectionDefaultTpl>(
    const context::Model &,
    typename LieGroupMap::template operationProduct<context::Scalar, context::Options>::type &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void
  lieGroup<context::Scalar, context::Options, JointCollectionDefaultTpl>(
    const context::Model &,
    typename LieGroupMap::template operationProduct<context::Scalar, context::Options>::type &);

#ifndef PINOCCHIO_SKIP_CASADI_UNSUPPORTED

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI bool isNormalized<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs>(
    const context::Model &, const Eigen::MatrixBase<context::VectorXs> &, const context::Scalar &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI bool
  isNormalized<context::Scalar, context::Options, JointCollectionDefaultTpl, context::VectorXs>(
    const context::Model &, const Eigen::MatrixBase<context::VectorXs> &, const context::Scalar &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI bool isSameConfiguration<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const context::Scalar &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI bool isSameConfiguration<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const context::Scalar &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void integrateCoeffWiseJacobian<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI void integrateCoeffWiseJacobian<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::MatrixXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::MatrixXs> &);

#endif // PINOCCHIO_SKIP_CASADI_UNSUPPORTED

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs integrate<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs integrate<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs interpolate<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const context::Scalar &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs interpolate<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const context::Scalar &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs difference<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs difference<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs squaredDistance<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs squaredDistance<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs randomConfiguration<
    LieGroupMap,
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs randomConfiguration<
    context::Scalar,
    context::Options,
    JointCollectionDefaultTpl,
    context::VectorXs,
    context::VectorXs>(
    const context::Model &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const Eigen::MatrixBase<context::VectorXs> &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs
  randomConfiguration<LieGroupMap, context::Scalar, context::Options, JointCollectionDefaultTpl>(
    const context::Model &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs
  randomConfiguration<context::Scalar, context::Options, JointCollectionDefaultTpl>(
    const context::Model &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs
  neutral<LieGroupMap, context::Scalar, context::Options, JointCollectionDefaultTpl>(
    const context::Model &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI context::VectorXs
  neutral<context::Scalar, context::Options, JointCollectionDefaultTpl>(const context::Model &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI
    typename LieGroupMap::template operationProduct<context::Scalar, context::Options>::type
    lieGroup<LieGroupMap, context::Scalar, context::Options, JointCollectionDefaultTpl>(
      const context::Model &);

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI
    typename LieGroupMap::template operationProduct<context::Scalar, context::Options>::type
    lieGroup<context::Scalar, context::Options, JointCollectionDefaultTpl>(const context::Model &);
} // namespace pinocchio
