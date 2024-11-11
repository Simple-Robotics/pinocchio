//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_fictious_constraint_hpp__
#define __pinocchio_algorithm_constraints_fictious_constraint_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, FictiousConstraintModelTpl<Scalar, Options>>
  {
    typedef FictiousConstraintModelTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<FictiousConstraintModelTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef FictiousConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef FictiousConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef void ConstraintSet;
  };

  template<typename _Scalar, int _Options>
  struct traits<FictiousConstraintDataTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef FictiousConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef FictiousConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef void ConstraintSet;
  };

  /// \brief Fictious constraint model used for variant definition
  template<typename Scalar, int Options>
  struct FictiousConstraintModelTpl
  : ConstraintModelBase<FictiousConstraintModelTpl<Scalar, Options>>
  {
    FictiousConstraintModelTpl()
    {
    }

    bool operator==(const FictiousConstraintModelTpl &) const
    {
      return true;
    }

    inline FictiousConstraintDataTpl<Scalar, Options> createData() const;
  };

  /// \brief Fictious constraint data used for variant definition
  template<typename Scalar, int Options>
  struct FictiousConstraintDataTpl : ConstraintDataBase<FictiousConstraintDataTpl<Scalar, Options>>
  {
    FictiousConstraintDataTpl()
    {
    }

    bool operator==(const FictiousConstraintDataTpl &) const
    {
      return true;
    }
  };

  template<typename Scalar, int Options>
  inline FictiousConstraintDataTpl<Scalar, Options>
  FictiousConstraintModelTpl<Scalar, Options>::createData() const
  {
    return FictiousConstraintDataTpl<Scalar, Options>();
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_fictious_constraint_hpp__
