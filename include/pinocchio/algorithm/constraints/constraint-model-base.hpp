//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_model_base_hpp__
#define __pinocchio_algorithm_constraints_constraint_model_base_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/common/model-entity.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-vector-parameters.hpp"

template<typename Scalar>
struct BaumgarteCorrectorParametersTpl;

namespace pinocchio
{

  enum struct ConstraintFormulationLevel
  {
    POSITION_LEVEL,    // scaling dt^2
    VELOCITY_LEVEL,    // scaling dt
    ACCELERATION_LEVEL // scaling 1
  };

  template<class Derived>
  struct ConstraintModelBase
  : NumericalBase<Derived>
  , ModelEntity<Derived>
  {
    typedef typename traits<Derived>::Scalar Scalar;
    enum
    {
      Options = traits<Derived>::Options
    };

    typedef ModelEntity<Derived> Base;

    typedef typename traits<Derived>::ConstraintData ConstraintData;
    typedef typename traits<Derived>::ConstraintSet ConstraintSet;
    typedef typename traits<Derived>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Derived>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;
    typedef typename traits<Derived>::ActiveComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef typename traits<Derived>::ActiveComplianceVectorTypeConstRef
      ActiveComplianceVectorTypeConstRef;
    typedef typename traits<Derived>::BaumgarteCorrectorVectorParametersRef
      BaumgarteCorrectorVectorParametersRef;
    typedef typename traits<Derived>::BaumgarteCorrectorVectorParametersConstRef
      BaumgarteCorrectorVectorParametersConstRef;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> BooleanVector;
    typedef std::vector<Eigen::DenseIndex> EigenIndexVector;

    using Base::createData;

    Derived & derived()
    {
      return static_cast<Derived &>(*this);
    }
    const Derived & derived() const
    {
      return static_cast<const Derived &>(*this);
    }

    template<typename NewScalar>
    typename CastType<NewScalar, Derived>::type cast() const
    {
      return derived().template cast<NewScalar>();
    }

    template<typename OtherDerived>
    void cast(ConstraintModelBase<OtherDerived> & other) const
    {
      other.name = name;
    }

    /// \brief Resize the constraint if needed at the current state given by data and store the
    /// results in cdata.
    template<template<typename, int> class JointCollectionTpl>
    void resize(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata)
    {
      derived().resize_impl(model, data, cdata);
    }

    /// \brief Evaluate the constraint values at the current state given by data and store the
    /// results in cdata.
    template<int Options, template<typename, int> class JointCollectionTpl>
    void calc(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      derived().calc(model, data, cdata);
    }

    template<int Options, template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobian(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & jacobian_matrix) const
    {
      derived().jacobian(model, data, cdata, jacobian_matrix.const_cast_derived());
    }

    template<int Options, template<typename, int> class JointCollectionTpl>
    typename traits<Derived>::JacobianMatrixType jacobian(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      typedef typename traits<Derived>::JacobianMatrixType ReturnType;
      ReturnType res = ReturnType::Zero(activeSize(), model.nv);

      jacobian(model, data, cdata, res);

      return res;
    }

    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Derived>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      return derived().jacobianMatrixProduct(model, data, cdata, mat.derived());
    }

    template<
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op = SETTO>
    void jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & res,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      derived().jacobianMatrixProduct(
        model, data, cdata, mat.derived(), res.const_cast_derived(), aot);
    }

    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Derived>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      return derived().jacobianTransposeMatrixProduct(model, data, cdata, mat.derived());
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
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & res,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      derived().jacobianTransposeMatrixProduct(
        model, data, cdata, mat.derived(), res.const_cast_derived(), aot);
    }

    ///
    /// \brief Map the constraint forces (aka constraint Lagrange multipliers) to joint space (e.g.,
    /// joint forces, joint torque vector).
    ///
    /// \param[in] model The model of the rigid body system.
    /// \param[in] data The data associated with model.
    /// \param[in] cdata The constraint data associated with the constraint model.
    /// \param[in] constraint_forces Input constraint forces (Lagrange multipliers) associated with
    /// the constraint.
    /// \param[out] joint_forces Output joint forces associated with each joint of the model.
    /// \param[out] joint_torques Output joint torques associated with the model.
    /// \param[in] reference_frame Input reference frame in which the forces are expressed.
    ///
    /// \note The results will be added to the joint_forces and joint_torques ouput argument.
    ///
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
      derived().mapConstraintForceToJointSpace(
        model, data, cdata, constraint_forces, joint_forces, joint_torques.const_cast_derived(),
        reference_frame);
    }

    ///
    /// \brief Map the joint space quantities (e.g.,
    /// joint motions, joint motion vector) to the constraint motions.
    ///
    /// \param[in] model The model of the rigid body system.
    /// \param[in] data The data associated with model.
    /// \param[in] cdata The constraint data associated with the constraint model.
    /// \param[in] joint_motions Input joint motions associated with the model.
    /// \param[in] joint_generalized_velocity Input joint motions associated with the model.
    /// \param[out] constraint_motions Output constraint motions.
    /// \param[in] reference_frame Input reference frame in which the joint motions are expressed.
    ///
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
      derived().mapJointSpaceToConstraintMotion(
        model, data, cdata, joint_motions, joint_generalized_velocity, constraint_motions,
        reference_frame);
    }

    // Attributes common to all constraints

    /// \brief Name of the constraint
    std::string name;

    template<typename OtherDerived>
    bool operator==(const ConstraintModelBase<OtherDerived> & other) const
    {
      return name == other.name;
    }

    template<typename OtherDerived>
    bool operator!=(const ConstraintModelBase<OtherDerived> & other) const
    {
      return !(*this == other);
    }

    template<typename OtherDerived>
    ConstraintModelBase & operator=(const ConstraintModelBase<OtherDerived> & other)
    {
      name = other.name;

      return *this;
    }

    ConstraintData createData() const
    {
      return derived().createData();
    }

    /// \brief Returns the colwise sparsity associated with a given row
    const BooleanVector & getRowActivableSparsityPattern(const Eigen::Index row_id) const
    {
      return derived().getRowActivableSparsityPattern(row_id);
    }

    /// \brief Returns the colwise sparsity associated with a given row of the active set of
    /// cosntraints
    const BooleanVector & getRowActiveSparsityPattern(const Eigen::Index row_id) const
    {
      return derived().getRowActiveSparsityPattern(row_id);
    }

    /// \brief Returns the vector of the activable indexes associated with a given row
    const EigenIndexVector & getRowActivableIndexes(const Eigen::DenseIndex row_id) const
    {
      return derived().getRowActivableIndexes(row_id);
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    const EigenIndexVector & getRowActiveIndexes(const Eigen::DenseIndex row_id) const
    {
      return derived().getRowActiveIndexes(row_id);
    }

    /// \brief Returns the active compliance internally stored in the constraint and corresponding
    /// to the active set contained in cdata
    ActiveComplianceVectorTypeConstRef getActiveCompliance() const
    {
      return derived().getActiveCompliance_impl();
    }

    /// \brief Returns the active compliance internally stored in the constraint and corresponding
    /// to the active set contained in cdata
    ActiveComplianceVectorTypeRef getActiveCompliance()
    {
      return derived().getActiveCompliance_impl();
    }

    int size() const
    {
      return derived().size();
    }

    int activeSize() const
    {
      return derived().activeSize();
    }

    ConstraintSet & set()
    {
      return derived().set();
    }
    const ConstraintSet & set() const
    {
      return derived().set();
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
      derived().appendCouplingConstraintInertias(
        model, data, cdata, diagonal_constraint_inertia.derived(), reference_frame);
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ComplianceVectorTypeConstRef compliance() const
    {
      return derived().compliance_impl();
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ComplianceVectorTypeRef compliance()
    {
      return derived().compliance_impl();
    }

    // CHOICE: right now we use the scalar Baumgarte
    // /// \brief Returns the Baumgarte vector parameters internally stored in the constraint model
    // BaumgarteCorrectorVectorParametersConstRef baumgarte_corrector_vector_parameters() const
    // {
    //   return derived().baumgarte_corrector_vector_parameters_impl();
    // }

    // /// \brief Returns the Baumgarte vector parameters internally stored in the constraint model
    // BaumgarteCorrectorVectorParametersRef baumgarte_corrector_vector_parameters()
    // {
    //   return derived().baumgarte_corrector_vector_parameters_impl();
    // }

    /// \brief Returns the Baumgarte parameters internally stored in the constraint model
    const BaumgarteCorrectorParameters & baumgarte_corrector_parameters() const
    {
      return derived().baumgarte_corrector_parameters_impl();
    }

    /// \brief Returns the Baumgarte parameters internally stored in the constraint model
    BaumgarteCorrectorParameters & baumgarte_corrector_parameters()
    {
      return derived().baumgarte_corrector_parameters_impl();
    }

    /// \brief Default implementation: do nothing
    template<template<typename, int> class JointCollectionTpl>
    void resize_impl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata)
    {
      PINOCCHIO_UNUSED_VARIABLE(model);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
    }

    ConstraintModelBase & base()
    {
      return *this;
    }

    const ConstraintModelBase & base() const
    {
      return *this;
    }

    std::string shortname() const
    {
      return derived().shortname();
    }
    static std::string classname()
    {
      return Derived::classname();
    }

    void disp(std::ostream & os) const
    {
      using namespace std;
      os << shortname() << endl;
    }

    friend std::ostream &
    operator<<(std::ostream & os, const ConstraintModelBase<Derived> & constraint)
    {
      constraint.disp(os);
      return os;
    }

  protected:
    template<int Options, template<typename, int> class JointCollectionTpl>
    explicit ConstraintModelBase(const ModelTpl<Scalar, Options, JointCollectionTpl> & /*model*/)
    {
    }

    /// \brief Default constructor
    ConstraintModelBase()
    {
    }
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_constraint_model_base_hpp__
