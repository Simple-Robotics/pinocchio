//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_fwd_hpp__
#define __pinocchio_algorithm_constraints_fwd_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include <boost/variant.hpp>

namespace pinocchio
{
  // Constraints
  template<typename Scalar, int Options = 0>
  struct RigidConstraintModelTpl;
  template<typename Scalar, int Options = 0>
  struct RigidConstraintDataTpl;

  template<typename Scalar, int Options = 0>
  struct FictiousConstraintModelTpl;
  template<typename Scalar, int Options = 0>
  struct FictiousConstraintDataTpl;

  template<typename Scalar, int Options = 0>
  struct FrictionalJointConstraintModelTpl;
  typedef FrictionalJointConstraintModelTpl<context::Scalar> FrictionalJointConstraintModel;

  template<typename Scalar, int Options = 0>
  struct FrictionalJointConstraintDataTpl;
  typedef FrictionalJointConstraintDataTpl<context::Scalar> FrictionalJointConstraintData;

  template<typename Scalar, int Options = 0>
  struct JointLimitConstraintModelTpl;
  typedef JointLimitConstraintModelTpl<context::Scalar> JointLimitConstraintModel;

  template<typename Scalar, int Options = 0>
  struct JointLimitConstraintDataTpl;
  typedef JointLimitConstraintDataTpl<context::Scalar> JointLimitConstraintData;

  template<typename Scalar, int Options = 0>
  struct BilateralPointConstraintModelTpl;
  typedef BilateralPointConstraintModelTpl<context::Scalar> BilateralPointConstraintModel;
  template<typename Scalar, int Options = 0>
  struct BilateralPointConstraintDataTpl;
  typedef BilateralPointConstraintDataTpl<context::Scalar> BilateralPointConstraintData;

  template<typename Scalar, int Options = 0>
  struct FrictionalPointConstraintModelTpl;
  typedef FrictionalPointConstraintModelTpl<context::Scalar> FrictionalPointConstraintModel;
  template<typename Scalar, int Options = 0>
  struct FrictionalPointConstraintDataTpl;
  typedef FrictionalPointConstraintDataTpl<context::Scalar> FrictionalPointConstraintData;

  template<typename Scalar, int Options = 0>
  struct WeldConstraintModelTpl;
  typedef WeldConstraintModelTpl<context::Scalar> WeldConstraintModel;
  template<typename Scalar, int Options = 0>
  struct WeldConstraintDataTpl;
  typedef WeldConstraintDataTpl<context::Scalar> WeldConstraintData;

  template<typename Scalar, int Options = 0>
  struct ConstraintCollectionDefaultTpl;

  typedef ConstraintCollectionDefaultTpl<context::Scalar, context::Options>
    ConstraintCollectionDefault;

  template<
    typename Scalar,
    int _Options,
    template<typename S, int O> class ConstraintCollectionTpl = ConstraintCollectionDefaultTpl>
  struct ConstraintModelTpl;
  typedef ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>
    ConstraintModel;

  template<
    typename Scalar,
    int _Options,
    template<typename S, int O> class ConstraintCollectionTpl = ConstraintCollectionDefaultTpl>
  struct ConstraintDataTpl;
  typedef ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>
    ConstraintData;

  // Sets
  template<typename Scalar, int Options = 0>
  struct BoxSetTpl;
  typedef BoxSetTpl<context::Scalar> BoxSet;

  template<typename Scalar, int Options = 0>
  struct UnboundedSetTpl;
  typedef UnboundedSetTpl<context::Scalar> UnboundedSet;

  // Convex sets
  template<typename Scalar>
  struct CoulombFrictionConeTpl;
  typedef CoulombFrictionConeTpl<context::Scalar> CoulombFrictionCone;

  template<typename Scalar>
  struct DualCoulombFrictionConeTpl;
  typedef DualCoulombFrictionConeTpl<context::Scalar> DualCoulombFrictionCone;

  template<typename Scalar>
  struct PositiveOrthantConeTpl;
  typedef PositiveOrthantConeTpl<context::Scalar> PositiveOrthantCone;

  template<typename Scalar>
  struct NegativeOrthantConeTpl;
  typedef NegativeOrthantConeTpl<context::Scalar> NegativeOrthantCone;

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_fwd_hpp__
