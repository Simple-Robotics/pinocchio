//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_collection_default_hpp__
#define __pinocchio_algorithm_constraints_constraint_collection_default_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include <boost/variant.hpp>

namespace pinocchio
{
  template<typename _Scalar, int _Options>
  struct ConstraintCollectionDefaultTpl
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    //    typedef FictiousConstraintModelTpl<Scalar, Options> FictiousConstraintModel;
    //    typedef FictiousConstraintDataTpl<Scalar, Options> FictiousConstraintData;

    typedef BilateralPointConstraintModelTpl<Scalar, Options> BilateralPointConstraintModel;
    typedef BilateralPointConstraintDataTpl<Scalar, Options> BilateralPointConstraintData;

    typedef FrictionalPointConstraintModelTpl<Scalar, Options> FrictionalPointConstraintModel;
    typedef FrictionalPointConstraintDataTpl<Scalar, Options> FrictionalPointConstraintData;

    typedef FrictionalJointConstraintModelTpl<Scalar, Options> FrictionalJointConstraintModel;
    typedef FrictionalJointConstraintDataTpl<Scalar, Options> FrictionalJointConstraintData;

    typedef JointLimitConstraintModelTpl<Scalar, Options> JointLimitConstraintModel;
    typedef JointLimitConstraintDataTpl<Scalar, Options> JointLimitConstraintData;

    typedef WeldConstraintModelTpl<Scalar, Options> WeldConstraintModel;
    typedef WeldConstraintDataTpl<Scalar, Options> WeldConstraintData;

    typedef boost::variant<
      boost::blank,
      //      FictiousConstraintModel,
      BilateralPointConstraintModel,
      FrictionalPointConstraintModel,
      FrictionalJointConstraintModel,
      JointLimitConstraintModel,
      WeldConstraintModel>
      ConstraintModelVariant;

    typedef boost::variant<
      boost::blank,
      //      FictiousConstraintData,
      BilateralPointConstraintData,
      FrictionalPointConstraintData,
      FrictionalJointConstraintData,
      JointLimitConstraintData,
      WeldConstraintData>
      ConstraintDataVariant;
  }; // struct ConstraintCollectionDefaultTpl

  typedef ConstraintCollectionDefault::ConstraintModelVariant ConstraintModelVariant;
  typedef ConstraintCollectionDefault::ConstraintDataVariant ConstraintDataVariant;

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_constraint_collection_default_hpp__
