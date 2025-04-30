//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_unary_constraint_base_hpp__
#define __pinocchio_algorithm_constraints_unary_constraint_base_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct UnaryConstraintModelBase : ConstraintModelBase<Derived>
  {
    typedef typename traits<Derived>::Scalar Scalar;
    typedef ConstraintModelBase<Derived> Base;

    using Base::derived;
    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    template<int Options, template<typename, int> class JointCollectionTpl>
    UnaryConstraintModelBase(const ModelTpl<Scalar, Options, JointCollectionTpl> & model)
    : Base(model)
    {
    }

    template<typename OtherDerived>
    bool operator==(const UnaryConstraintModelBase<OtherDerived> & other) const
    {
      return base() == other.base();
    }

    template<typename OtherDerived>
    bool operator!=(const UnaryConstraintModelBase<OtherDerived> & other) const
    {
      return !(*this == other);
    }

  protected:
    /// \brief Default constructor
    UnaryConstraintModelBase()
    {
    }
  }; // struct UnaryConstraintModelBase
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_unary_constraint_base_hpp__
