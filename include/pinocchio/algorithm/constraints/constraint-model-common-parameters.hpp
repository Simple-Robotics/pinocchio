//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_model_common_parameters_hpp__
#define __pinocchio_algorithm_constraints_constraint_model_common_parameters_hpp__

namespace pinocchio
{

  template<typename _BaumgarteVector>
  struct BaumgarteCorrectorVectorParametersTpl;

  template<typename _Scalar>
  struct BaumgarteCorrectorParametersTpl;

  template<typename Derived>
  struct ConstraintModelCommonParameters
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    template<typename OtherDerived>
    friend struct ConstraintModelCommonParameters;

    typedef ConstraintModelCommonParameters<Derived> Self;
    typedef typename traits<Derived>::Scalar Scalar;
    typedef typename traits<Derived>::ComplianceVectorType ComplianceVectorType;
    typedef typename traits<Derived>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Derived>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;
    typedef typename traits<Derived>::BaumgarteVectorType BaumgarteVectorType;
    typedef typename traits<Derived>::BaumgarteCorrectorVectorParameters
      BaumgarteCorrectorVectorParameters;
    typedef typename traits<Derived>::BaumgarteCorrectorVectorParametersRef
      BaumgarteCorrectorVectorParametersRef;
    typedef typename traits<Derived>::BaumgarteCorrectorVectorParametersConstRef
      BaumgarteCorrectorVectorParametersConstRef;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    /// \brief Cast to NewScalar
    template<typename NewScalar, typename OtherDerived>
    void cast(ConstraintModelCommonParameters<OtherDerived> & other) const
    {
      other.m_compliance = m_compliance.template cast<NewScalar>();

      // CHOICE: right now we use the scalar Baumgarte
      // other.m_baumgarte_vector_parameters = m_baumgarte_vector_parameters.template
      // cast<NewScalar>();
      other.m_baumgarte_parameters = m_baumgarte_parameters.template cast<NewScalar>();
    }

    /// \brief Comparison operator
    bool operator==(const Self & other) const
    {
      // CHOICE: right now we use the scalar Baumgarte
      // return m_compliance == other.m_compliance
      //        && m_baumgarte_vector_parameters == other.m_baumgarte_vector_parameters;
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

    // CHOICE: right now we use the scalar Baumgarte
    // /// \brief Returns the Baumgarte vector parameters internally stored in the constraint model
    // BaumgarteCorrectorVectorParametersConstRef baumgarte_corrector_vector_parameters_impl() const
    // {
    //   return m_baumgarte_vector_parameters;
    // }

    // /// \brief Returns the Baumgarte vector parameters internally stored in the constraint model
    // BaumgarteCorrectorVectorParametersRef baumgarte_corrector_vector_parameters_impl()
    // {
    //   return m_baumgarte_vector_parameters;
    // }

    /// \brief Returns the Baumgarte parameters internally stored in the constraint model
    const BaumgarteCorrectorParameters & baumgarte_corrector_parameters_impl() const
    {
      return m_baumgarte_parameters;
    }

    /// \brief Returns the Baumgarte parameters internally stored in the constraint model
    BaumgarteCorrectorParameters & baumgarte_corrector_parameters_impl()
    {
      return m_baumgarte_parameters;
    }

  protected:
    /// \brief Default constructor - protected so that the class cannot be instanciated on its own.
    ConstraintModelCommonParameters()
    {
    }

    ComplianceVectorType m_compliance;
    // CHOICE: right now we use the scalar Baumgarte
    // BaumgarteCorrectorVectorParameters m_baumgarte_vector_parameters;
    BaumgarteCorrectorParameters m_baumgarte_parameters;
  };

} // namespace pinocchio

#endif // __pinocchio_algorithm_constraints_constraint_model_common_parameters_hpp__
