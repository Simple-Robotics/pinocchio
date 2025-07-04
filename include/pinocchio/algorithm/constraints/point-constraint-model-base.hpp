//
// Copyright (c) 2019-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_point_constraint_model_base_hpp__
#define __pinocchio_algorithm_constraints_point_constraint_model_base_hpp__

#include <algorithm>

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/spatial/skew.hpp"
#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/kinematics-constraint-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-common-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-vector-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct PointConstraintModelBase;

  template<typename Derived>
  struct traits<PointConstraintModelBase<Derived>>
  {
    static constexpr ConstraintFormulationLevel constraint_formulation_level =
      ConstraintFormulationLevel::VELOCITY_LEVEL;
    static constexpr bool has_baumgarte_corrector = true;
    static constexpr bool has_baumgarte_corrector_vector = true;

    template<typename InputMatrix>
    struct JacobianMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::Matrix<Scalar, 3, InputMatrix::ColsAtCompileTime, InputMatrixPlain::Options>
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

  ///
  ///  \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename Derived>
  struct PointConstraintModelBase
  : KinematicsConstraintModelBase<Derived>
  , ConstraintModelCommonParameters<Derived>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef typename traits<Derived>::Scalar Scalar;
    enum
    {
      Options = traits<Derived>::Options
    };

    typedef KinematicsConstraintModelBase<Derived> Base;
    typedef ConstraintModelCommonParameters<Derived> BaseCommonParameters;
    typedef ConstraintModelBase<Derived> RootBase;

    template<typename OtherDerived>
    friend struct PointConstraintModelBase;

    static const ConstraintFormulationLevel constraint_formulation_level =
      traits<PointConstraintModelBase>::constraint_formulation_level;
    typedef typename traits<Derived>::ConstraintData ConstraintData;
    typedef typename traits<Derived>::ComplianceVectorType ComplianceVectorType;
    typedef typename traits<Derived>::ActiveComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef typename traits<Derived>::ActiveComplianceVectorTypeConstRef
      ActiveComplianceVectorTypeConstRef;
    typedef typename traits<Derived>::BaumgarteCorrectorVectorParameters
      BaumgarteCorrectorVectorParameters;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    using Base::derived;
    using Base::joint1_id;
    using Base::joint2_id;
    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    typedef SE3Tpl<Scalar, Options> SE3;
    typedef MotionTpl<Scalar, Options> Motion;
    typedef ForceTpl<Scalar, Options> Force;
    typedef Eigen::Matrix<Scalar, 3, 6, Options> Matrix36;
    typedef Eigen::Matrix<Scalar, 6, 6, Options> Matrix6;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef Eigen::Matrix<Scalar, 6, 1, Options> Vector6;
    typedef Vector3 VectorConstraintSize;

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    BaseCommonParameters & base_common_parameters()
    {
      return static_cast<BaseCommonParameters &>(*this);
    }
    const BaseCommonParameters & base_common_parameters() const
    {
      return static_cast<const BaseCommonParameters &>(*this);
    }

    /// \brief Position of attached point with respect to the frame of joint1.
    SE3 joint1_placement;

    /// \brief Position of attached point with respect to the frame of joint2.
    SE3 joint2_placement;

    /// \brief Desired constraint shift at position level
    Vector3 desired_constraint_offset;

    /// \brief Desired constraint velocity at velocity level
    Vector3 desired_constraint_velocity;

    /// \brief Desired constraint velocity at acceleration level
    Vector3 desired_constraint_acceleration;

    /// \brief Colwise sparsity pattern associated with joint 1.
    BooleanVector colwise_joint1_sparsity;

    /// \brief Colwise sparsity pattern associated with joint 2.
    BooleanVector colwise_joint2_sparsity;

    /// \brief Jointwise span indexes associated with joint 1.
    EigenIndexVector joint1_span_indexes;

    /// \brief Jointwise span indexes associated with joint 2.
    EigenIndexVector joint2_span_indexes;

    EigenIndexVector loop_span_indexes;

    /// \brief Sparsity pattern associated to the constraint;
    BooleanVector colwise_sparsity;

    /// \brief Indexes of the columns spanned by the constraints.
    EigenIndexVector colwise_span_indexes;

    /// \brief Dimensions of the model
    int nv;

    ///  \brief Depth of the kinematic tree for joint1 and joint2
    size_t depth_joint1, depth_joint2;

  protected:
    using BaseCommonParameters::m_compliance;
    // CHOICE: right now we use the scalar Baumgarte
    // using BaseCommonParameters::m_baumgarte_vector_parameters;
    using BaseCommonParameters::m_baumgarte_parameters;

  public:
    ///
    ///  \brief Default constructor
    ///
    PointConstraintModelBase()
    : nv(-1)
    , depth_joint1(0)
    , depth_joint2(0)
    {
    }

    //  public:
    ///
    ///  \brief Contructor with from a given type, joint indexes and placements.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] model Model associated to the constraint.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint2_id Index of the joint 2 in the model tree.
    /// \param[in] joint1_placement Placement of the constraint w.r.t the frame of joint1.
    /// \param[in] joint2_placement Placement of the constraint w.r.t the frame of joint2.
    /// \param[in] reference_frame Reference frame in which the constraints quantities are
    /// expressed.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement,
      const JointIndex joint2_id,
      const SE3 & joint2_placement)
    : Base(model, joint1_id, joint2_id)
    , joint1_placement(joint1_placement)
    , joint2_placement(joint2_placement)
    , desired_constraint_offset(Vector3::Zero())
    , desired_constraint_velocity(Vector3::Zero())
    , desired_constraint_acceleration(Vector3::Zero())
    , colwise_joint1_sparsity(model.nv)
    , colwise_joint2_sparsity(model.nv)
    , loop_span_indexes((size_t)model.nv)
    {
      init(model);
    }

    ///
    ///  \brief Contructor with from a given type, joint1_id and placement.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint1_placement Placement of the constraint w.r.t the frame of joint1.
    /// \param[in] reference_frame Reference frame in which the constraints quantities are
    /// expressed.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement)
    : Base(model, joint1_id, 0)
    , joint1_placement(joint1_placement)
    , joint2_placement(SE3::Identity())
    , desired_constraint_offset(Vector3::Zero())
    , desired_constraint_velocity(Vector3::Zero())
    , desired_constraint_acceleration(Vector3::Zero())
    , colwise_joint1_sparsity(model.nv)
    , colwise_joint2_sparsity(model.nv)
    , loop_span_indexes((size_t)model.nv)
    {
      init(model);
    }

    ///
    ///  \brief Contructor with from a given type and the joint ids.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint2_id Index of the joint 2 in the model tree.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const JointIndex joint2_id)
    : Base(model, joint1_id, joint2_id)
    , joint1_placement(SE3::Identity())
    , joint2_placement(SE3::Identity())
    , desired_constraint_offset(Vector3::Zero())
    , desired_constraint_velocity(Vector3::Zero())
    , desired_constraint_acceleration(Vector3::Zero())
    , colwise_joint1_sparsity(model.nv)
    , colwise_joint2_sparsity(model.nv)
    , loop_span_indexes((size_t)model.nv)
    {
      init(model);
    }

    ///
    ///  \brief Contructor with from a given type and .
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    ///
    /// \remarks The second joint id (joint2_id) is set to be 0 (corresponding to the index of the
    /// universe).
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model, const JointIndex joint1_id)
    : Base(model, joint1_id, 0)
    , joint1_placement(SE3::Identity())
    , joint2_placement(SE3::Identity())
    , desired_constraint_offset(Vector3::Zero())
    , desired_constraint_velocity(Vector3::Zero())
    , desired_constraint_acceleration(Vector3::Zero())
    , colwise_joint1_sparsity(model.nv)
    , colwise_joint2_sparsity(model.nv)
    , loop_span_indexes((size_t)model.nv)
    {
      init(model);
    }

    ///
    /// \brief Create data storage associated to the constraint
    ///
    ConstraintData createData() const
    {
      return ConstraintData(*this);
    }

    /// \brief Returns the colwise sparsity associated with a given row
    const BooleanVector & getRowActivableSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return colwise_sparsity;
    }

    /// \brief Returns the sparsity associated with a given row
    const BooleanVector & getRowActiveSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      return getRowActivableSparsityPattern(row_id);
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    const EigenIndexVector & getRowActivableIndexes(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return colwise_span_indexes;
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    const EigenIndexVector & getRowActiveIndexes(const Eigen::DenseIndex row_id) const
    {
      return getRowActivableIndexes(row_id);
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ActiveComplianceVectorTypeConstRef getActiveCompliance_impl() const
    {
      return this->compliance();
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ActiveComplianceVectorTypeRef getActiveCompliance_impl()
    {
      return this->compliance();
    }

    ///
    ///  \brief Comparison operator
    ///
    /// \param[in] other Other PointConstraintModelBase to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributs
    /// must be the same).
    ///
    bool operator==(const PointConstraintModelBase & other) const
    {
      if (this == &other)
        return true;
      return base() == other.base() && base_common_parameters() == other.base_common_parameters()
             && joint1_id == other.joint1_id && joint2_id == other.joint2_id
             && joint1_placement == other.joint1_placement
             && joint2_placement == other.joint2_placement && nv == other.nv
             && desired_constraint_offset == other.desired_constraint_offset
             && desired_constraint_velocity == other.desired_constraint_velocity
             && desired_constraint_acceleration == other.desired_constraint_acceleration
             && colwise_joint1_sparsity == other.colwise_joint1_sparsity
             && colwise_joint2_sparsity == other.colwise_joint2_sparsity
             && joint1_span_indexes == other.joint1_span_indexes
             && joint2_span_indexes == other.joint2_span_indexes && nv == other.nv
             && depth_joint1 == other.depth_joint1 && depth_joint2 == other.depth_joint2
             && colwise_sparsity == other.colwise_sparsity
             && colwise_span_indexes == other.colwise_span_indexes
             && loop_span_indexes == other.loop_span_indexes;
    }

    ///
    ///  \brief Oposite of the comparison operator.
    ///
    /// \param[in] other Other PointConstraintModelBase to compare with.
    ///
    /// \returns false if the two *this is not equal to other (at least type, joint1_id or placement
    /// attributs is different).
    ///
    bool operator!=(const PointConstraintModelBase & other) const
    {
      return !(*this == other);
    }

    /// \brief Evaluate the constraint values at the current state given by data and store the
    /// results in cdata.
    template<template<typename, int> class JointCollectionTpl>
    void calc(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(model);

      if (joint1_id > 0)
        cdata.oMc1 = data.oMi[joint1_id] * joint1_placement;
      else
        cdata.oMc1 = joint1_placement;

      if (joint2_id > 0)
        cdata.oMc2 = data.oMi[joint2_id] * joint2_placement;
      else
        cdata.oMc2 = joint2_placement;

      // Compute relative placement
      cdata.c1Mc2 = cdata.oMc1.actInv(cdata.oMc2);
      const auto & _1R2_ = cdata.c1Mc2.rotation();

      // Compute errors
      auto & position_error = cdata.constraint_position_error;
      position_error.noalias() = cdata.c1Mc2.translation();
      //      cdata.constraint_position_error = cdata.oMc1.inverse().translation();

      const auto vf1 = joint1_placement.actInv(data.v[this->joint1_id]);
      const auto vf2 = joint2_placement.actInv(data.v[this->joint2_id]);

      auto & velocity_error = cdata.constraint_velocity_error;
      const Vector3 velocity_error_component1 = _1R2_ * vf2.linear() - vf1.linear();
      velocity_error.noalias() = velocity_error_component1 - vf1.angular().cross(position_error);

      const auto af1 = joint1_placement.actInv(data.a[this->joint1_id]);
      const auto af2 = joint2_placement.actInv(data.a[this->joint2_id]);
      auto & acceleration_error = cdata.constraint_acceleration_error;
      acceleration_error.noalias() = _1R2_ * (af2.linear() + vf2.angular().cross(vf2.linear()))
                                     - (af1.linear() + vf1.angular().cross(vf1.linear()));
      acceleration_error.noalias() -= af1.angular().cross(position_error);
      acceleration_error.noalias() += vf1.angular().cross(vf1.angular().cross(position_error));
      acceleration_error.noalias() -= 2 * vf1.angular().cross(velocity_error_component1);

      cdata.A1_world = this->getA1(cdata, WorldFrameTag());
      cdata.A2_world = this->getA2(cdata, WorldFrameTag());
      cdata.A_world = cdata.A1_world + cdata.A2_world;

      cdata.A1_local = this->getA1(cdata, LocalFrameTag());
      cdata.A2_local = this->getA2(cdata, LocalFrameTag());
      cdata.A_local = cdata.A1_local + cdata.A2_local;
    }

    /// \brief Returns the constraint projector associated with joint 1.
    /// This matrix transforms a spatial velocity expressed at the origin to the first component of
    /// the constraint associated with joint 1.
    template<ReferenceFrame rf>
    Matrix36 getA1(const ConstraintData & cdata, ReferenceFrameTag<rf>) const
    {
      Matrix36 res;

      if (std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value)
      {
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = oM1.rotation().transpose() * v_tmp;

        const SE3 & oM1 = cdata.oMc1;
        Vector3 v_tmp;
        res.template leftCols<3>() = -oM1.rotation().transpose();
        INTERNAL_LOOP(0, -oM1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, -oM1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, -oM1.translation(), res.template rightCols<3>());

        for (int i = 0; i < 3; ++i)
        {
          res.template rightCols<3>().col(i) +=
            cdata.constraint_position_error.cross(oM1.rotation().transpose().col(i));
        }

#undef INTERNAL_LOOP
      }
      else if (std::is_same<ReferenceFrameTag<rf>, LocalFrameTag>::value)
      {
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = iM1.rotation().transpose() * v_tmp;

        const SE3 & iM1 = this->joint1_placement;
        Vector3 v_tmp;
        res.template leftCols<3>() = -iM1.rotation().transpose();
        INTERNAL_LOOP(0, -iM1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, -iM1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, -iM1.translation(), res.template rightCols<3>());

        for (int i = 0; i < 3; ++i)
        {
          res.template rightCols<3>().col(i) +=
            cdata.constraint_position_error.cross(iM1.rotation().transpose().col(i));
        }
      }
      else
      {
        assert(false && "Should never happened");
      }

#undef INTERNAL_LOOP

      return res;
    }

    /// \brief Returns the constraint projector associated with joint 2.
    /// This matrix transforms a spatial velocity expressed at the origin to the first component of
    /// the constraint associated with joint 2.
    template<ReferenceFrame rf>
    Matrix36 getA2(const ConstraintData & cdata, ReferenceFrameTag<rf>) const
    {
      Matrix36 res;
      typedef typename SE3::Vector3 Vector3;

      if (std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value)
      {
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = oM1.rotation().transpose() * v_tmp;

        const SE3 & oM1 = cdata.oMc1;
        const SE3 & oM2 = cdata.oMc2;
        res.template leftCols<3>() = oM1.rotation().transpose();
        Vector3 v_tmp;
        INTERNAL_LOOP(0, oM2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, oM2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, oM2.translation(), res.template rightCols<3>());

#undef INTERNAL_LOOP
      }
      else if (std::is_same<ReferenceFrameTag<rf>, LocalFrameTag>::value)
      {
        const SE3 & j2Mc2 = this->joint2_placement;
        const SE3 & c1Mc2 = cdata.c1Mc2;
        const typename SE3::Matrix3 c1Rj2 = c1Mc2.rotation() * j2Mc2.rotation().transpose();
        res.template leftCols<3>() = c1Rj2;
        Vector3 v_tmp;
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = c1Rj2 * v_tmp;

        INTERNAL_LOOP(0, j2Mc2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, j2Mc2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, j2Mc2.translation(), res.template rightCols<3>());

#undef INTERNAL_LOOP
      }
      else
      {
        assert(false && "Should never happened");
      }

      return res;
    }

    ///
    /// @brief This function computes the spatial inertia associated with the constraint.
    /// This function is useful to express the constraint inertia associated with the constraint for
    /// AL-based approaches.
    ///
    template<typename Vector3Like>
    Matrix6 computeConstraintSpatialInertia(
      const SE3Tpl<Scalar, Options> & placement,
      const Eigen::MatrixBase<Vector3Like> & diagonal_constraint_inertia) const
    {
      EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Vector3Like, Vector3);
      Matrix6 res;

      const auto & R = placement.rotation();
      const auto & t = placement.translation();

      typedef Eigen::Matrix<Scalar, 3, 3, Options> Matrix3;
      const Matrix3 R_Sigma = R * diagonal_constraint_inertia.asDiagonal();
      const Matrix3 t_skew = skew(t);

      auto block_LL = res.template block<3, 3>(SE3::LINEAR, SE3::LINEAR);
      auto block_LA = res.template block<3, 3>(SE3::LINEAR, SE3::ANGULAR);
      auto block_AL = res.template block<3, 3>(SE3::ANGULAR, SE3::LINEAR);
      auto block_AA = res.template block<3, 3>(SE3::ANGULAR, SE3::ANGULAR);

      block_LL.noalias() = R_Sigma * R.transpose();
      block_LA.noalias() = -block_LL * t_skew;
      block_AL.noalias() = block_LA.transpose();
      block_AA.noalias() = t_skew * block_LA;

      return res;
    }

    template<
      typename Matrix6LikeOut1,
      typename Matrix6LikeOut2,
      typename Matrix6LikeOut3,
      ReferenceFrame rf>
    void computeConstraintInertias(
      const ConstraintData & cdata,
      const Scalar & constraint_inertia_value,
      const Eigen::MatrixBase<Matrix6LikeOut1> & I11,
      const Eigen::MatrixBase<Matrix6LikeOut2> & I12,
      const Eigen::MatrixBase<Matrix6LikeOut3> & I22,
      const ReferenceFrameTag<rf> reference_frame) const
    {
      computeConstraintInertias(
        cdata, Vector3::Constant(constraint_inertia_value), I11.const_cast_derived(),
        I12.const_cast_derived(), I22.const_cast_derived(), reference_frame);
    }

    template<
      typename Vector3Like,
      typename Matrix6LikeOut1,
      typename Matrix6LikeOut2,
      typename Matrix6LikeOut3,
      ReferenceFrame rf>
    void computeConstraintInertias(
      const ConstraintData & cdata,
      const Eigen::MatrixBase<Vector3Like> & diagonal_constraint_inertia,
      const Eigen::MatrixBase<Matrix6LikeOut1> & I11,
      const Eigen::MatrixBase<Matrix6LikeOut2> & I12,
      const Eigen::MatrixBase<Matrix6LikeOut3> & I22,
      const ReferenceFrameTag<rf> reference_frame) const
    {
      EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Vector3Like, Vector3);
      EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Matrix6LikeOut1, Matrix6);
      EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Matrix6LikeOut2, Matrix6);
      EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Matrix6LikeOut3, Matrix6);
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);
      //      assert((check_expression_if_real<Scalar,
      //      true>(diagonal_constraint_inertia.isZero(Scalar(0)))));

      const auto & A1 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A1_world : cdata.A1_local;
      const auto & A2 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A2_world : cdata.A2_local;

      Matrix36 diagonal_constraint_inertia_time_A;

      if (joint1_id > 0)
      {
        diagonal_constraint_inertia_time_A.noalias() =
          diagonal_constraint_inertia.asDiagonal() * A1;
        I11.const_cast_derived().noalias() = A1.transpose() * diagonal_constraint_inertia_time_A;
      }
      else
        I11.const_cast_derived().setZero();

      if (joint2_id > 0)
      {
        diagonal_constraint_inertia_time_A.noalias() =
          diagonal_constraint_inertia.asDiagonal() * A2;
        I22.const_cast_derived().noalias() = A2.transpose() * diagonal_constraint_inertia_time_A;
      }
      else
        I22.const_cast_derived().setZero();

      // Compute the cross coupling term
      if (joint1_id > 0 && joint2_id > 0)
      {
        I12.const_cast_derived().noalias() = A1.transpose() * diagonal_constraint_inertia_time_A;
      }
      else
        I12.const_cast_derived().setZero();
    }

    template<
      template<typename, int> class JointCollectionTpl,
      typename Vector3Like,
      ReferenceFrame rf>
    void appendCouplingConstraintInertias(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<Vector3Like> & diagonal_constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_UNUSED_VARIABLE(model);

      Matrix6 I11, I12, I22;
      computeConstraintInertias(cdata, diagonal_constraint_inertia, I11, I12, I22, reference_frame);
      assert(
        (std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value
         || std::is_same<ReferenceFrameTag<rf>, LocalFrameTag>::value)
        && "must never happened");

      Matrix6 & Y1 = std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value
                       ? data.oYaba_augmented[joint1_id]
                       : data.oYaba_augmented[joint1_id];

      if (joint1_id > 0)
        Y1 += I11;

      Matrix6 & Y2 = std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value
                       ? data.oYaba_augmented[joint2_id]
                       : data.oYaba_augmented[joint2_id];

      if (joint2_id > 0)
        Y2 += I22;

      if (joint1_id > 0 && joint2_id > 0)
      {
        if (joint1_id < joint2_id)
        {
          data.joint_cross_coupling.get({joint1_id, joint2_id}) += I12;
        }
        else
        {
          data.joint_cross_coupling.get({joint2_id, joint1_id}) += I12.transpose();
        }
      }
    }

    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Derived>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef typename traits<Derived>::template JacobianMatrixProductReturnType<InputMatrix>::type
        ReturnType;
      ReturnType res(3, mat.cols());
      jacobianMatrixProduct(model, data, cdata, mat.derived(), res);
      return res;
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
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef typename Data::Vector3 Vector3;
      OutputMatrix & res = _res.const_cast_derived();

      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), model.nv);
      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.cols(), res.cols());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), size());
      PINOCCHIO_UNUSED_VARIABLE(aot);

      if (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
        res.setZero();

      //      const Eigen::DenseIndex constraint_dim = size();
      //
      //      const Eigen::DenseIndex
      //      complexity_strategy_1 = 6 * res.cols() * 36 + constraint_dim * 36 * res.cols(),
      //      complexity_strategy_2 = 36 * constraint_dim * 6 + constraint_dim * 36 * res.cols();

      const auto & A1 = cdata.A1_world;
      const auto & A2 = cdata.A2_world;

      const auto & A = cdata.A_world;
      for (Eigen::DenseIndex jj = 0; jj < model.nv; ++jj)
      {
        if (!(colwise_joint1_sparsity[jj] || colwise_joint2_sparsity[jj]))
          continue;
        Vector3 AxSi;

        typedef typename Data::Matrix6x::ConstColXpr ConstColXpr;
        const ConstColXpr Jcol = data.J.col(jj);

        if (colwise_joint1_sparsity[jj] && colwise_joint2_sparsity[jj])
        {
          AxSi.noalias() = A * Jcol;
        }
        else if (colwise_joint1_sparsity[jj])
          AxSi.noalias() = A1 * Jcol;
        else
          AxSi.noalias() = A2 * Jcol;

        if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
          res.noalias() -= AxSi * mat.row(jj);
        else // AddTo, SetTo
          res.noalias() += AxSi * mat.row(jj);
      }
    }

    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Derived>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef typename traits<Derived>::template JacobianTransposeMatrixProductReturnType<
        InputMatrix>::type ReturnType;
      ReturnType res(model.nv, mat.cols());
      jacobianTransposeMatrixProduct(model, data, cdata, mat.derived(), res);
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
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef typename Data::Vector3 Vector3;
      OutputMatrix & res = _res.const_cast_derived();

      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), size());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.cols(), mat.cols());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), model.nv);
      PINOCCHIO_UNUSED_VARIABLE(aot);

      if (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
        res.setZero();

      const auto & A1 = cdata.A1_world;
      const auto & A2 = cdata.A2_world;

      const auto & A = cdata.A_world;
      for (Eigen::DenseIndex jj = 0; jj < model.nv; ++jj)
      {
        if (!(colwise_joint1_sparsity[jj] || colwise_joint2_sparsity[jj]))
          continue;
        Vector3 AxSi;

        typedef typename Data::Matrix6x::ConstColXpr ConstColXpr;
        const ConstColXpr Jcol = data.J.col(jj);

        if (colwise_joint1_sparsity[jj] && colwise_joint2_sparsity[jj])
        {
          AxSi.noalias() = A * Jcol;
        }
        else if (colwise_joint1_sparsity[jj])
          AxSi.noalias() = A1 * Jcol;
        else
          AxSi.noalias() = A2 * Jcol;

        if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
          res.row(jj).noalias() -= AxSi.transpose() * mat;
        else
          res.row(jj).noalias() += AxSi.transpose() * mat;
      }
    }

    using RootBase::jacobian;

    ///  \brief Evaluate the Jacobian associated to the constraint at the given state stored in data
    /// and cdata.  The results Jacobian is evaluated in the jacobian input/output matrix.
    template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobian(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const
    {
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      JacobianMatrix & jacobian_matrix = _jacobian_matrix.const_cast_derived();

      //      const PointConstraintModelBase & cmodel = *this;

      const SE3 & oMc1 = cdata.oMc1;
      const SE3 & oMc2 = cdata.oMc2;
      const SE3 & c1Mc2 = cdata.c1Mc2;
      const auto & position_error = cdata.constraint_position_error;

      for (Eigen::DenseIndex jj = 0; jj < model.nv; ++jj)
      {
        if (colwise_joint1_sparsity[jj] || colwise_joint2_sparsity[jj])
        {
          typedef typename Data::Matrix6x::ConstColXpr ConstColXpr;
          const ConstColXpr Jcol = data.J.col(jj);
          const MotionRef<const ConstColXpr> Jcol_motion(Jcol);

          jacobian_matrix.col(jj).setZero();
          if (colwise_joint1_sparsity[jj])
          {
            const Motion Jcol_local(oMc1.actInv(Jcol_motion)); // TODO: simplify computations
            jacobian_matrix.col(jj).noalias() -= Jcol_local.linear();
            jacobian_matrix.col(jj).noalias() += -Jcol_local.angular().cross(position_error);
          }

          if (colwise_joint2_sparsity[jj])
          {
            const Motion Jcol_local(oMc2.actInv(Jcol_motion)); // TODO: simplify computations
            jacobian_matrix.col(jj) += c1Mc2.rotation() * Jcol_local.linear();
          }
        }
      }
    }

    ///
    /// \copydoc Base::mapConstraintForceToJointForces(const ModelTpl<Scalar, Options,
    /// JointCollectionTpl> &, const DataTpl<Scalar, Options, JointCollectionTpl> &, const
    /// ConstraintData &, const Eigen::MatrixBase<ForceLike> &, std::vector<ForceTpl<Scalar,
    /// Options>, ForceAllocator> &, ReferenceFrameTag<rf>)
    ///
    template<
      template<typename, int> class JointCollectionTpl,
      typename ForceLike,
      typename ForceAllocator,
      ReferenceFrame rf>
    void mapConstraintForceToJointForces(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
      ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_forces.size(), size_t(model.njoints));
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), size());
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);

      // Todo: optimize code
      const auto & A1 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A1_world : cdata.A1_local;
      const auto & A2 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A2_world : cdata.A2_local;

      if (joint1_id > 0)
        joint_forces[this->joint1_id].toVector().noalias() += A1.transpose() * constraint_forces;
      if (joint2_id > 0)
        joint_forces[this->joint2_id].toVector().noalias() += A2.transpose() * constraint_forces;
    }

    ///
    /// \copydoc Base::mapJointMotionsToConstraintMotion(const ModelTpl<Scalar, Options,
    /// JointCollectionTpl> &, const DataTpl<Scalar, Options, JointCollectionTpl> &, const
    /// ConstraintData &, const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> &, const
    /// Eigen::MatrixBase<VectorLike> &, ReferenceFrameTag<rf>)
    ///
    template<
      template<typename, int> class JointCollectionTpl,
      typename MotionAllocator,
      typename VectorLike,
      ReferenceFrame rf>
    void mapJointMotionsToConstraintMotion(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_accelerations,
      const Eigen::MatrixBase<VectorLike> & constraint_motion,
      ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_accelerations.size(), size_t(model.njoints));
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motion.rows(), activeSize());
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);

      const auto & A1 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A1_world : cdata.A1_local;
      const auto & A2 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A2_world : cdata.A2_local;

      if (joint1_id > 0 && joint2_id > 0)
        constraint_motion.const_cast_derived().noalias() =
          A1 * joint_accelerations[this->joint1_id].toVector()
          + A2 * joint_accelerations[this->joint2_id].toVector();
      else if (joint1_id > 0)
        constraint_motion.const_cast_derived().noalias() =
          A1 * joint_accelerations[this->joint1_id].toVector();
      else if (joint2_id > 0)
        constraint_motion.const_cast_derived().noalias() =
          A2 * joint_accelerations[this->joint2_id].toVector();
      else
        constraint_motion.const_cast_derived().setZero();
    }

    static int size()
    {
      return 3;
    }

    static int activeSize()
    {
      return size();
    }

    /// \returns An expression of *this with the Scalar type casted to NewScalar.
    template<typename NewScalar, typename OtherDerived>
    void cast(PointConstraintModelBase<OtherDerived> & res) const
    {
      Base::cast(res);
      BaseCommonParameters::template cast<NewScalar>(res);

      res.joint1_id = joint1_id;
      res.joint2_id = joint2_id;
      res.joint1_placement = joint1_placement.template cast<NewScalar>();
      res.joint2_placement = joint2_placement.template cast<NewScalar>();
      res.desired_constraint_offset = desired_constraint_offset.template cast<NewScalar>();
      res.desired_constraint_velocity = desired_constraint_velocity.template cast<NewScalar>();
      res.desired_constraint_acceleration =
        desired_constraint_acceleration.template cast<NewScalar>();
      res.colwise_joint1_sparsity = colwise_joint1_sparsity;
      res.colwise_joint2_sparsity = colwise_joint2_sparsity;
      res.joint1_span_indexes = joint1_span_indexes;
      res.joint2_span_indexes = joint2_span_indexes;
      res.colwise_sparsity = colwise_sparsity;
      res.colwise_span_indexes = colwise_span_indexes;
      res.nv = nv;
      res.depth_joint1 = depth_joint1;
      res.depth_joint2 = depth_joint2;
      res.loop_span_indexes = loop_span_indexes;
    }

  protected:
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void init(const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model)
    {
      nv = model.nv;
      depth_joint1 = static_cast<size_t>(model.supports[joint1_id].size());
      depth_joint2 = static_cast<size_t>(model.supports[joint2_id].size());

      typedef ModelTpl<Scalar, OtherOptions, JointCollectionTpl> Model;
      typedef typename Model::JointModel JointModel;
      static const bool default_sparsity_value = false;
      colwise_joint1_sparsity.fill(default_sparsity_value);
      colwise_joint2_sparsity.fill(default_sparsity_value);

      joint1_span_indexes.reserve(size_t(model.njoints));
      joint2_span_indexes.reserve(size_t(model.njoints));

      JointIndex current1_id = 0;
      if (joint1_id > 0)
        current1_id = joint1_id;

      JointIndex current2_id = 0;
      if (joint2_id > 0)
        current2_id = joint2_id;

      while (current1_id != current2_id)
      {
        if (current1_id > current2_id)
        {
          const JointModel & joint1 = model.joints[current1_id];
          joint1_span_indexes.push_back((Eigen::DenseIndex)current1_id);
          Eigen::DenseIndex current1_col_id = joint1.idx_v();
          for (int k = 0; k < joint1.nv(); ++k, ++current1_col_id)
          {
            colwise_joint1_sparsity[current1_col_id] = true;
          }
          current1_id = model.parents[current1_id];
        }
        else
        {
          const JointModel & joint2 = model.joints[current2_id];
          joint2_span_indexes.push_back((Eigen::DenseIndex)current2_id);
          Eigen::DenseIndex current2_col_id = joint2.idx_v();
          for (int k = 0; k < joint2.nv(); ++k, ++current2_col_id)
          {
            colwise_joint2_sparsity[current2_col_id] = true;
          }
          current2_id = model.parents[current2_id];
        }
      }
      assert(current1_id == current2_id && "current1_id should be equal to current2_id");

      {
        JointIndex current_id = current1_id;
        while (current_id > 0)
        {
          const JointModel & joint = model.joints[current_id];
          joint1_span_indexes.push_back((Eigen::DenseIndex)current_id);
          joint2_span_indexes.push_back((Eigen::DenseIndex)current_id);
          Eigen::DenseIndex current_row_id = joint.idx_v();
          for (int k = 0; k < joint.nv(); ++k, ++current_row_id)
          {
            colwise_joint1_sparsity[current_row_id] = true;
            colwise_joint2_sparsity[current_row_id] = true;
          }
          current_id = model.parents[current_id];
        }
      }
      std::reverse(joint1_span_indexes.begin(), joint1_span_indexes.end());
      std::reverse(joint2_span_indexes.begin(), joint2_span_indexes.end());
      colwise_span_indexes.reserve((size_t)model.nv);
      colwise_sparsity.resize(model.nv);
      colwise_sparsity.setZero();
      loop_span_indexes.reserve((size_t)model.nv);
      for (Eigen::DenseIndex col_id = 0; col_id < model.nv; ++col_id)
      {
        if (colwise_joint1_sparsity[col_id] || colwise_joint2_sparsity[col_id])
        {
          colwise_span_indexes.push_back(col_id);
          colwise_sparsity[col_id] = true;
        }

        if (colwise_joint1_sparsity[col_id] != colwise_joint2_sparsity[col_id])
        {
          loop_span_indexes.push_back(col_id);
        }
      }

      // Set compliance and baumgarte parameters
      m_compliance = ComplianceVectorType::Zero(size());
      // CHOICE: right now we use the scalar Baumgarte
      // m_baumgarte_vector_parameters = BaumgarteCorrectorVectorParameters(size());
      m_baumgarte_parameters = BaumgarteCorrectorParameters();
    }
  }; // PointConstraintModelBase<Derived>

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_point_constraint_model_base_hpp__
