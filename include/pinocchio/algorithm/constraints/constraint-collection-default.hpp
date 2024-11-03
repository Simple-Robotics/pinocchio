//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_collection_default_hpp__
#define __pinocchio_algorithm_constraints_constraint_collection_default_hpp__

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

    typedef BilateralPointConstraintModelTpl<Scalar, Options> BilateralPointConstraintModel;
    typedef BilateralPointConstraintDataTpl<Scalar, Options> BilateralPointConstraintData;

    typedef FrictionalPointConstraintModelTpl<Scalar, Options> FrictionalPointConstraintModel;
    typedef FrictionalPointConstraintDataTpl<Scalar, Options> FrictionalPointConstraintData;

    typedef FrictionalJointConstraintModelTpl<Scalar, Options> FrictionalJointConstraintModel;
    typedef FrictionalJointConstraintDataTpl<Scalar, Options> FrictionalJointConstraintData;

    typedef boost::variant<
      BilateralPointConstraintModel,
      FrictionalPointConstraintModel,
      FrictionalJointConstraintModel>
      ConstraintModelVariant;

    typedef boost::variant<
      BilateralPointConstraintData,
      FrictionalPointConstraintData,
      FrictionalJointConstraintData>
      ConstraintDataVariant;
  }; // struct ConstraintCollectionDefaultTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_constraint_collection_default_hpp__
