//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_binary_constraint_base_hpp__
#define __pinocchio_algorithm_constraints_binary_constraint_base_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct BinaryConstraintModelBase : ConstraintModelBase<Derived>
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
    BinaryConstraintModelBase(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const JointIndex joint2_id)
    : Base(model)
    , joint1_id(joint1_id)
    , joint2_id(joint2_id)
    {
    }

    /// \brief Index of the first joint in the model tree
    JointIndex joint1_id;

    /// \brief Index of the second joint in the model tree
    JointIndex joint2_id;

    template<typename OtherDerived>
    bool operator==(const BinaryConstraintModelBase<OtherDerived> & other) const
    {
      return base() == other.base() && joint1_id == other.joint1_id && joint2_id == other.joint2_id;
    }

    template<typename OtherDerived>
    bool operator!=(const BinaryConstraintModelBase<OtherDerived> & other) const
    {
      return !(*this == other);
    }

  protected:
    /// \brief Default constructor
    BinaryConstraintModelBase()
    {
    }
  }; // struct BinaryConstraintModelBase
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_binary_constraint_base_hpp__
