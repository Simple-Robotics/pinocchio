//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_model_common_parameters_hpp__
#define __pinocchio_algorithm_constraints_constraint_model_common_parameters_hpp__

namespace pinocchio
{

  template<typename _BaumgarteVector>
  struct BaumgarteCorrectorParametersTpl;

  template<typename Derived>
  struct ConstraintModelCommonParameters
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef ConstraintModelCommonParameters<Derived> Self;

    typedef typename traits<Derived>::ComplianceVectorType ComplianceVectorType;
    typedef typename traits<Derived>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Derived>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;
    typedef typename traits<Derived>::BaumgarteVectorType BaumgarteVectorType;
    typedef typename traits<Derived>::BaumgarteCorrectorParameters BaumgarteCorrectorParameters;
    typedef
      typename traits<Derived>::BaumgarteCorrectorParametersRef BaumgarteCorrectorParametersRef;
    typedef typename traits<Derived>::BaumgarteCorrectorParametersConstRef
      BaumgarteCorrectorParametersConstRef;

    template<typename OtherDerived>
    friend struct ConstraintModelCommonParameters;

    /// \brief Cast to NewScalar
    template<typename NewScalar, typename OtherDerived>
    void cast(ConstraintModelCommonParameters<OtherDerived> & other) const
    {
      other.m_compliance = m_compliance.template cast<NewScalar>();
      other.m_baumgarte_parameters = m_baumgarte_parameters.template cast<NewScalar>();
    }

    /// \brief Comparison operator
    bool operator==(const Self & other) const
    {
      return m_compliance == other.m_compliance
             && m_baumgarte_parameters == other.m_baumgarte_parameters;
    }

    /// \brief Comparison operator
    bool operator!=(const Self & other) const
    {
      return !(*this == other);
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ComplianceVectorTypeConstRef compliance_impl() const
    {
      return m_compliance;
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ComplianceVectorTypeRef compliance_impl()
    {
      return m_compliance;
    }

    /// \brief Returns the baumgarte parameters internally stored in the constraint model
    BaumgarteCorrectorParametersConstRef baumgarte_corrector_parameters_impl() const
    {
      return m_baumgarte_parameters;
    }

    /// \brief Returns the baumgarte parameters internally stored in the constraint model
    BaumgarteCorrectorParametersRef baumgarte_corrector_parameters_impl()
    {
      return m_baumgarte_parameters;
    }

  protected:
    /// \brief Default constructor - protected so that the class cannot be instanciated on its own.
    ConstraintModelCommonParameters()
    {
    }

    ComplianceVectorType m_compliance;
    BaumgarteCorrectorParameters m_baumgarte_parameters;
  };

} // namespace pinocchio

#endif // __pinocchio_algorithm_constraints_constraint_model_common_parameters_hpp__
