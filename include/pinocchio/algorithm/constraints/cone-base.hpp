//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_cone_base_hpp__
#define __pinocchio_algorithm_constraints_cone_base_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/set-base.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct ConeBase : SetBase<Derived>
  {
    typedef typename traits<Derived>::Scalar Scalar;
    typedef typename traits<Derived>::DualCone DualCone;

    typedef SetBase<Derived> Base;

    using Base::derived;
    using Base::dim;
    using Base::isInside;
    using Base::project;

    DualCone dual() const
    {
      return derived().dual();
    }

  }; // struct ConeBase

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_cone_base_hpp__
