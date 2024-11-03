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
  struct FrictionalJointConstraintModelTpl;
  typedef FrictionalJointConstraintModelTpl<context::Scalar> FrictionalJointConstraintModel;

  template<typename Scalar, int Options = 0>
  struct FrictionalJointConstraintDataTpl;
  typedef FrictionalJointConstraintDataTpl<context::Scalar> FrictionalJointConstraintData;

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

  template<typename _Scalar, int _Options>
  struct ConstraintCollectionDefaultTpl
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef RigidConstraintModelTpl<Scalar, Options> RigidConstraintModel;
    typedef RigidConstraintDataTpl<Scalar, Options> RigidConstraintData;

    typedef FrictionalJointConstraintModelTpl<Scalar, Options> FrictionalJointConstraintModel;
    typedef FrictionalJointConstraintDataTpl<Scalar, Options> FrictionalJointConstraintData;

    typedef boost::variant<RigidConstraintModel, FrictionalJointConstraintModel>
      ConstraintModelVariant;
    typedef boost::variant<RigidConstraintData, FrictionalJointConstraintData>
      ConstraintDataVariant;
  }; // struct ConstraintCollectionDefaultTpl

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
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_fwd_hpp__
