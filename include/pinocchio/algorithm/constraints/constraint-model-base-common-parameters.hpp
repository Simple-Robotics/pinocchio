//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraint_model_common_parameters_hpp__
#define __pinocchio_algorithm_constraint_model_common_parameters_hpp__

#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include <iostream>

namespace pinocchio
{

  template<typename Derived>
  struct ConstraintModelBaseCommonParameters : ConstraintModelBase<Derived>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef typename traits<Derived>::Scalar Scalar;
    typedef ConstraintModelBaseCommonParameters<Derived> Self;
    typedef ConstraintModelBase<Derived> Base;

    typedef typename traits<Derived>::ComplianceVectorType ComplianceVectorType;
    typedef typename traits<Derived>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Derived>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;

    using Base::size;

    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    template<typename OtherDerived>
    friend struct ConstraintModelBaseCommonParameters;

    /// \brief Cast to NewScalar
    template<typename NewScalar, typename OtherDerived>
    void cast(ConstraintModelBaseCommonParameters<OtherDerived> & other) const
    {
      Base::cast(other);
      other.m_compliance = m_compliance.template cast<NewScalar>();
    }

    /// \brief Comparison operator
    bool operator==(const Self & other) const
    {
      return base() == other.base() && m_compliance == other.m_compliance;
    }

    /// \brief Comparison operator
    bool operator!=(const Self & other) const
    {
      return !(*this == other);
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ComplianceVectorTypeConstRef compliance() const
    {
      return m_compliance;
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ComplianceVectorTypeRef compliance()
    {
      return m_compliance;
    }

  protected:
    template<int Options, template<typename, int> class JointCollectionTpl>
    explicit ConstraintModelBaseCommonParameters(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & /*model*/)
    {
    }

    /// \brief Default constructor - protected so that the class cannot be instanciated on its own.
    ConstraintModelBaseCommonParameters()
    {
    }

    ComplianceVectorType m_compliance;
  };

} // namespace pinocchio

#endif // __pinocchio_algorithm_constraint_model_common_parameters_hpp__
