//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_fictious_constraint_hpp__
#define __pinocchio_algorithm_constraints_fictious_constraint_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"

namespace pinocchio
{

  template<typename NewScalar>
  struct CastType<NewScalar, FictiousConstraintModel>
  {
    typedef FictiousConstraintModel type;
  };

  /// \brief Fictious constraint model used for variant definition
  struct FictiousConstraintModel
  {
    FictiousConstraintModel()
    {
    }
  };

  /// \brief Fictious constraint data used for variant definition
  struct FictiousConstraintData
  {
    FictiousConstraintData()
    {
    }
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_fictious_constraint_hpp__
