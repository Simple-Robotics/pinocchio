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

    template<typename VectorLikeIn, typename VectorLikeIn2, typename VectorLikeOut>
    void scaledProject_impl(
      const Eigen::MatrixBase<VectorLikeIn> & x,
      const Eigen::MatrixBase<VectorLikeIn2> & scale,
      const Eigen::MatrixBase<VectorLikeOut> & x_proj) const
    {
      assert(x.size() == scale.size() && " x and scale should have the same size.");
      assert(
        scale.isApprox(scale(0) * VectorLikeIn2::Ones(scale.size()))
        && "Only scalar scaling are supported.");
      PINOCCHIO_UNUSED_VARIABLE(scale); // the cone is preserved when scaled by a scalar
      return project(x, x_proj);
    }

  }; // struct ConeBase

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_cone_base_hpp__
