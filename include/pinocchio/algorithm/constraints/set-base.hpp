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

    Eigen::DenseIndex dim() const
    {
      return derived().dim();
    }

    template<typename VectorLike>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike)
      project(const Eigen::MatrixBase<VectorLike> & x) const
    {
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike) ReturnType;
      ReturnType res(x.size());
      derived().project(x, res);
      return res;
    }

    template<typename VectorLikeIn, typename VectorLikeOut>
    void project(
      const Eigen::MatrixBase<VectorLikeIn> & x,
      const Eigen::MatrixBase<VectorLikeOut> & x_proj) const
    {
      return derived().projec(x.derived(), x_proj.const_cast_derived());
    }

    template<typename Vector>
    bool isInside(const Eigen::MatrixBase<Vector> & x, Scalar prec = Scalar(0)) const
    {
      return derived().isInside(x, prec);
    }
  }; // struct SetBase

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_set_base_hpp__
