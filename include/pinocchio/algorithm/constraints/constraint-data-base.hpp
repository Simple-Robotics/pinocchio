//
// Copyright (c) 2023-2024 INRIA
//

#ifndef __pinocchio_algorithm_constraint_data_base_hpp__
#define __pinocchio_algorithm_constraint_data_base_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/common/data-entity.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct ConstraintDataBase
  : NumericalBase<Derived>
  , DataEntity<Derived>
  {
    typedef typename traits<Derived>::Scalar Scalar;
    typedef typename traits<Derived>::ConstraintModel ConstraintModel;
    typedef DataEntity<Derived> Base;

    Derived & derived()
    {
      return static_cast<Derived &>(*this);
    }
    const Derived & derived() const
    {
      return static_cast<const Derived &>(*this);
    }
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraint_data_base_hpp__
