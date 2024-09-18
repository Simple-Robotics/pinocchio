//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_set_base_hpp__
#define __pinocchio_algorithm_constraints_set_base_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct SetBase
  {
    typedef typename traits<Derived>::Scalar Scalar;

    Derived & derived()
    {
      return static_cast<Derived &>(*this);
    }

    const Derived & derived() const
    {
      return static_cast<const Derived &>(*this);
    }

    int dim() const
    {
      return derived().dim();
    }

    template<typename Vector>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(Vector) project(const Eigen::MatrixBase<Vector> & x) const
    {
      return derived().project(x);
    }

    template<typename Vector>
    bool isInside(const Eigen::MatrixBase<Vector> & x, Scalar prec = Scalar(0)) const
    {
      return derived().isInside(x, prec);
    }
  }; // struct SetBase

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_set_base_hpp__
