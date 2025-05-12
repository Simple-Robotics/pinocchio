//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraint_model_generic_hpp__
#define __pinocchio_algorithm_constraint_model_generic_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-generic.hpp"
#include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-vector-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"

namespace pinocchio
{

  template<
    typename _Scalar,
    int _Options,
    template<typename S, int O> class ConstraintCollectionTpl>
  struct traits<ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    static constexpr bool has_baumgarte_corrector = true;
    static constexpr bool has_baumgarte_corrector_vector = true;

    typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;
    typedef ConstraintData Data;
    typedef boost::blank ConstraintSet;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> JacobianMatrixType;
    typedef VectorXs VectorConstraintSize;

    typedef VectorXs ComplianceVectorType;
    typedef Eigen::Ref<ComplianceVectorType> ComplianceVectorTypeRef;
    typedef Eigen::Ref<const ComplianceVectorType> ComplianceVectorTypeConstRef;

    typedef ComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef ComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;

    typedef VectorXs BaumgarteVectorType;
    typedef Eigen::Ref<VectorXs> BaumgarteVectorTypeRef;
    typedef BaumgarteCorrectorVectorParametersTpl<BaumgarteVectorTypeRef>
      BaumgarteCorrectorVectorParameters;
    typedef BaumgarteCorrectorVectorParameters BaumgarteCorrectorVectorParametersRef;
    typedef Eigen::Ref<const VectorXs> BaumgarteVectorTypeConstRef;
    typedef BaumgarteCorrectorVectorParametersTpl<BaumgarteVectorTypeConstRef>
      BaumgarteCorrectorVectorParametersConstRef;

    template<typename InputMatrix>
    struct JacobianMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::
        Matrix<Scalar, Eigen::Dynamic, InputMatrix::ColsAtCompileTime, InputMatrixPlain::Options>
          type;
    };

    template<typename InputMatrix>
    struct JacobianTransposeMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::Matrix<
        Scalar,
        Eigen::Dynamic,
        InputMatrixPlain::ColsAtCompileTime,
        InputMatrixPlain::Options>
        type;
    };
  };

  template<
    typename _Scalar,
    int _Options,
    template<typename S, int O> class ConstraintCollectionTpl>
  struct ConstraintModelTpl
  : ConstraintModelBase<ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl>>
  , ConstraintCollectionTpl<_Scalar, _Options>::ConstraintModelVariant
  , serialization::Serializable<ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef ConstraintModelTpl Self;
    typedef ConstraintModelBase<Self> Base;
    typedef ConstraintModelBase<Self> RootBase;

    typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;
    typedef ConstraintCollectionTpl<Scalar, Options> ConstraintCollection;
    typedef typename ConstraintCollection::ConstraintDataVariant ConstraintDataVariant;
    typedef typename ConstraintCollection::ConstraintModelVariant ConstraintModelVariant;
    typedef typename traits<Self>::ComplianceVectorType ComplianceVectorType;
    typedef typename traits<Self>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Self>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;
    typedef typename traits<Self>::ActiveComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef
      typename traits<Self>::ActiveComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;
    typedef
      typename traits<Self>::BaumgarteCorrectorVectorParameters BaumgarteCorrectorVectorParameters;
    typedef typename traits<Self>::BaumgarteCorrectorVectorParametersRef
      BaumgarteCorrectorVectorParametersRef;
    typedef typename traits<Self>::BaumgarteCorrectorVectorParametersConstRef
      BaumgarteCorrectorVectorParametersConstRef;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    using RootBase::jacobian;
    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    ConstraintModelTpl()
    : ConstraintModelVariant()
    {
    }

    ConstraintModelTpl(const ConstraintModelVariant & cmodel_variant)
    : ConstraintModelVariant(cmodel_variant)
    {
    }

    template<typename ContraintModelDerived>
    ConstraintModelTpl(const ConstraintModelBase<ContraintModelDerived> & cmodel)
    : ConstraintModelVariant((ConstraintModelVariant)cmodel.derived())
    {
      BOOST_MPL_ASSERT(
        (boost::mpl::contains<typename ConstraintModelVariant::types, ContraintModelDerived>));
    }

    ConstraintData createData() const
    {
      return ::pinocchio::visitors::createData<Scalar, Options, ConstraintCollectionTpl>(*this);
    }

    template<template<typename, int> class JointCollectionTpl>
    void calc(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      ::pinocchio::visitors::calc(*this, model, data, cdata);
    }

    template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobian(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & jacobian_matrix) const
    {
      ::pinocchio::visitors::jacobian(
        *this, model, data, cdata, jacobian_matrix.const_cast_derived());
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    const EigenIndexVector & getRowActivableIndexes(const Eigen::DenseIndex row_id) const
    {
      return ::pinocchio::visitors::getRowActivableIndexes(*this, row_id);
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    const EigenIndexVector & getRowActiveIndexes(const Eigen::DenseIndex row_id) const
    {
      return ::pinocchio::visitors::getRowActiveIndexes(*this, row_id);
    }

    /// \brief Returns the sparsity pattern associated with a given row
    const BooleanVector & getRowActivableSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      return ::pinocchio::visitors::getRowActivableSparsityPattern(*this, row_id);
    }

    /// \brief Returns the sparsity pattern associated with a given row
    const BooleanVector & getRowActiveSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      return ::pinocchio::visitors::getRowActiveSparsityPattern(*this, row_id);
    }

    /// \brief Returns the compliance associated to the current active set
    ActiveComplianceVectorTypeConstRef getActiveCompliance_impl() const
    {
      return ::pinocchio::visitors::getActiveCompliance(*this);
    }

    /// \brief Returns the compliance associated to the current active set
    ActiveComplianceVectorTypeRef getActiveCompliance_impl()
    {
      return ::pinocchio::visitors::getActiveCompliance(*this);
    }

    /// \brief Runs the underlying jacobian multiplication with a matrix.
    template<
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix>
    typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix) const
    {
      typedef typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
        ReturnType;
      ReturnType res(size(), input_matrix.cols());
      jacobianMatrixProduct(model, data, cdata, input_matrix.derived(), res);
      return res;
    }

    /// \brief Runs the underlying jacobian multiplication with a matrix.
    template<
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op = SETTO>
    void jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix,
      const Eigen::MatrixBase<OutputMatrix> & result_matrix,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      ::pinocchio::visitors::jacobianMatrixProduct(
        *this, model, data, cdata, input_matrix.derived(), result_matrix.const_cast_derived(), aot);
    }

    /// \brief Runs the underlying jacobian transpose multiplication with a matrix.
    template<
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix>
    typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix) const
    {
      typedef
        typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
          ReturnType;
      ReturnType res(model.nv, input_matrix.cols());
      jacobianTransposeMatrixProduct(*this, model, data, cdata, input_matrix.derived(), res);
      return res;
    }

    template<
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op = SETTO>
    void jacobianTransposeMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix,
      const Eigen::MatrixBase<OutputMatrix> & result_matrix,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      ::pinocchio::visitors::jacobianTransposeMatrixProduct(
        *this, model, data, cdata, input_matrix.derived(), result_matrix.const_cast_derived(), aot);
    }

    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorNLike,
      ReferenceFrame rf>
    void appendCouplingConstraintInertias(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame) const
    {
      ::pinocchio::visitors::appendCouplingConstraintInertias(
        *this, model, data, cdata, diagonal_constraint_inertia.derived(), reference_frame);
    }

    template<
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForceLike,
      typename ForceAllocator,
      typename JointTorquesLike,
      ReferenceFrame rf>
    void mapConstraintForceToJointSpace(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques,
      ReferenceFrameTag<rf> reference_frame) const
    {
      ::pinocchio::visitors::mapConstraintForceToJointSpace(
        *this, model, data, cdata, constraint_forces, joint_forces,
        joint_torques.const_cast_derived(), reference_frame);
    }

    template<
      template<typename, int> class JointCollectionTpl,
      typename MotionAllocator,
      typename JointMotionsLike,
      typename VectorLike,
      ReferenceFrame rf>
    void mapJointSpaceToConstraintMotion(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
      const Eigen::MatrixBase<JointMotionsLike> & joint_generalized_velocity,
      const Eigen::MatrixBase<VectorLike> & constraint_motions,
      ReferenceFrameTag<rf> reference_frame) const
    {
      ::pinocchio::visitors::mapJointSpaceToConstraintMotion(
        *this, model, data, cdata, joint_motions, joint_generalized_velocity,
        constraint_motions.const_cast_derived(), reference_frame);
    }

    static std::string classname()
    {
      return "ConstraintModel";
    }
    std::string shortname() const
    {
      return ::pinocchio::visitors::shortname(*this);
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ComplianceVectorTypeConstRef compliance_impl() const
    {
      return ::pinocchio::visitors::compliance(*this);
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ComplianceVectorTypeRef compliance_impl()
    {
      return ::pinocchio::visitors::compliance(*this);
    }

    // CHOICE: right now we use the scalar Baumgarte
    // /// \brief Returns the Baumgarte vector parameters internally stored in the constraint model
    // BaumgarteCorrectorVectorParametersConstRef baumgarte_corrector_vector_parameters_impl() const
    // {
    //   return ::pinocchio::visitors::getBaumgarteCorrectorVectorParameters(*this);
    // }

    // /// \brief Returns the Baumgarte vector parameters internally stored in the constraint model
    // BaumgarteCorrectorVectorParametersRef baumgarte_corrector_vector_parameters_impl()
    // {
    //   return ::pinocchio::visitors::getBaumgarteCorrectorVectorParameters(*this);
    // }

    /// \brief Returns the Baumgarte parameters internally stored in the constraint model
    const BaumgarteCorrectorParameters & baumgarte_corrector_parameters_impl() const
    {
      return ::pinocchio::visitors::getBaumgarteCorrectorParameters(*this);
    }

    /// \brief Returns the Baumgarte parameters internally stored in the constraint model
    BaumgarteCorrectorParameters & baumgarte_corrector_parameters_impl()
    {
      return ::pinocchio::visitors::getBaumgarteCorrectorParameters(*this);
    }

    /// \brief Returns the size of the constraint
    int size() const
    {
      return ::pinocchio::visitors::size(*this);
    }

    /// \brief Returns the size of the active constraints
    int activeSize() const
    {
      return ::pinocchio::visitors::activeSize(*this);
    }

    boost::blank & set()
    {
      static boost::blank val;
      PINOCCHIO_THROW_PRETTY(
        std::runtime_error, "Set method is not accessible for ConstraintModelTpl.");
      return val;
    }

    const boost::blank & set() const
    {
      static boost::blank val;
      PINOCCHIO_THROW_PRETTY(
        std::runtime_error, "Set method is not accessible for ConstraintModelTpl.");
      return val;
    }

    ConstraintModelVariant & toVariant()
    {
      return static_cast<ConstraintModelVariant &>(*this);
    }

    const ConstraintModelVariant & toVariant() const
    {
      return static_cast<const ConstraintModelVariant &>(*this);
    }

    template<typename ConstraintModelDerived>
    bool isEqual(const ConstraintModelBase<ConstraintModelDerived> & other) const
    {
      return ::pinocchio::isEqual(*this, other.derived());
    }

    bool isEqual(const ConstraintModelTpl & other) const
    {
      return toVariant() == other.toVariant();
    }

    bool operator==(const ConstraintModelTpl & other) const
    {
      return isEqual(other);
    }

    bool operator!=(const ConstraintModelTpl & other) const
    {
      return !(*this == other);
    }
  }; // struct ConstraintModelTpl
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraint_model_generic_hpp__
