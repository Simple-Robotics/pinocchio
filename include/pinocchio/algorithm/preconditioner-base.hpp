//
// Copyright (c) 2025 INRIA
//
#ifndef __pinocchio_preconditioner_base_hpp__
#define __pinocchio_preconditioner_base_hpp__

#include "pinocchio/algorithm/fwd.hpp"

namespace pinocchio
{

  template<typename PreconditionerDerived>
  struct PreconditionerBase
  {

    PreconditionerDerived & derived()
    {
      return static_cast<PreconditionerDerived &>(*this);
    }
    const PreconditionerDerived & derived() const
    {
      return static_cast<const PreconditionerDerived &>(*this);
    }

    explicit PreconditionerBase()
    {
    }

    template<typename MatrixIn, typename MatrixOut>
    void
    scale(const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & res) const
    {
      derived().scale(x.derived(), res.const_cast_derived());
    }

    template<typename MatrixIn>
    void scaleInPlace(const Eigen::MatrixBase<MatrixIn> & x) const
    {
      derived().scaleInPlace(x.derived());
    }

    template<typename MatrixIn, typename MatrixOut>
    void
    unscale(const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & res) const
    {
      derived().unscale(x.derived(), res.const_cast_derived());
    }

    template<typename MatrixIn>
    void unscaleInPlace(const Eigen::MatrixBase<MatrixIn> & x) const
    {
      derived().unscaleInPlace(x.derived());
    }

    Eigen::DenseIndex rows() const
    {
      return derived().rows();
    }
    Eigen::DenseIndex cols() const
    {
      return derived().cols();
    }

  }; // struct PreconditionerBase

} // namespace pinocchio

#endif // #ifndef __pinocchio_preconditioner_base_hpp__
