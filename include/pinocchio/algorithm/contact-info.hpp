//
// Copyright (c) 2019-2024 INRIA CNRS
//

#ifndef __pinocchio_algorithm_contact_info_hpp__
#define __pinocchio_algorithm_contact_info_hpp__

#include <algorithm>

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/spatial/skew.hpp"
#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"

namespace pinocchio
{
  ///  \brief Type of contact
  enum ContactType
  {
    CONTACT_3D = 0,   /// \brief Point contact model
    CONTACT_6D,       ///  \brief Frame contact model
    CONTACT_UNDEFINED ///  \brief The default contact is undefined
  };

  template<ContactType contact_type>
  struct contact_dim
  {
    enum
    {
      value = 0
    };
  };

  template<>
  struct contact_dim<CONTACT_3D>
  {
    enum
    {
      value = 3
    };
  };

  template<>
  struct contact_dim<CONTACT_6D>
  {
    enum
    {
      value = 6
    };
  };

  // TODO Remove when API is stabilized
  PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
  PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, RigidConstraintModelTpl<Scalar, Options>>
  {
    typedef RigidConstraintModelTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<RigidConstraintModelTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef RigidConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef boost::blank ConstraintSet;
  };

  template<typename _Scalar, int _Options>
  struct traits<RigidConstraintDataTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef RigidConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef boost::blank ConstraintSet;
  };

  ///
  ///  \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct PINOCCHIO_UNSUPPORTED_MESSAGE("The API will change towards more flexibility")
    RigidConstraintModelTpl : ConstraintModelBase<RigidConstraintModelTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef ConstraintModelBase<RigidConstraintModelTpl<_Scalar, _Options>> Base;

    template<typename NewScalar, int NewOptions>
    friend struct RigidConstraintModelTpl;

    using Base::base;

    typedef RigidConstraintModelTpl ContactModel;
    typedef RigidConstraintDataTpl<Scalar, Options> ContactData;
    typedef RigidConstraintDataTpl<Scalar, Options> ConstraintData;

    typedef SE3Tpl<Scalar, Options> SE3;
    typedef MotionTpl<Scalar, Options> Motion;
    typedef ForceTpl<Scalar, Options> Force;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;
    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;
    typedef Eigen::Matrix<Scalar, 3, 6, Options> Matrix36;
    typedef Eigen::Matrix<Scalar, 6, 6, Options> Matrix6;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef Eigen::Matrix<Scalar, 6, 1, Options> Vector6;

    ///  \brief Type of the contact.
    ContactType type;

    /// \brief Index of the first joint in the model tree
    JointIndex joint1_id;

    /// \brief Index of the second joint in the model tree
    JointIndex joint2_id;

    /// \brief Relative placement with respect to the frame of joint1.
    SE3 joint1_placement;

    /// \brief Relative placement with respect to the frame of joint2.
    SE3 joint2_placement;

    /// \brief Reference frame where the constraint is expressed (LOCAL_WORLD_ALIGNED or LOCAL)
    ReferenceFrame reference_frame;

    /// \brief Desired contact placement
    SE3 desired_contact_placement;

    /// \brief Desired contact spatial velocity
    Motion desired_contact_velocity;

    /// \brief Desired contact spatial acceleration
    Motion desired_contact_acceleration;

    ///  \brief Corrector parameters
    BaumgarteCorrectorParameters corrector;

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

    /// \brief Compliance associated with the contact model
    Vector3 compliance;

  protected:
    ///
    ///  \brief Default constructor
    ///
    RigidConstraintModelTpl()
    : nv(-1)
    , depth_joint1(0)
    , depth_joint2(0)
    {
    }

  public:
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
    RigidConstraintModelTpl(
      const ContactType type,
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement,
      const JointIndex joint2_id,
      const SE3 & joint2_placement,
      const ReferenceFrame & reference_frame = LOCAL)
    : Base(model)
    , type(type)
    , joint1_id(joint1_id)
    , joint2_id(joint2_id)
    , joint1_placement(joint1_placement)
    , joint2_placement(joint2_placement)
    , reference_frame(reference_frame)
    , desired_contact_placement(SE3::Identity())
    , desired_contact_velocity(Motion::Zero())
    , desired_contact_acceleration(Motion::Zero())
    , corrector(size())
    , colwise_joint1_sparsity(model.nv)
    , colwise_joint2_sparsity(model.nv)
    , joint1_span_indexes((size_t)model.njoints)
    , joint2_span_indexes((size_t)model.njoints)
    , loop_span_indexes((size_t)model.nv)
    , compliance(Vector3::Zero())
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
    RigidConstraintModelTpl(
      const ContactType type,
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement,
      const ReferenceFrame & reference_frame = LOCAL)
    : Base(model)
    , type(type)
    , joint1_id(joint1_id)
    , joint2_id(0)
    , joint1_placement(joint1_placement)
    , joint2_placement(SE3::Identity())
    , reference_frame(reference_frame)
    , desired_contact_placement(SE3::Identity())
    , desired_contact_velocity(Motion::Zero())
    , desired_contact_acceleration(Motion::Zero())
    , corrector(size())
    , colwise_joint1_sparsity(model.nv)
    , colwise_joint2_sparsity(model.nv)
    , joint1_span_indexes((size_t)model.njoints)
    , joint2_span_indexes((size_t)model.njoints)
    , loop_span_indexes((size_t)model.nv)
    , compliance(Vector3::Zero())
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
    RigidConstraintModelTpl(
      const ContactType type,
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const JointIndex joint2_id,
      const ReferenceFrame & reference_frame = LOCAL)
    : Base(model)
    , type(type)
    , joint1_id(joint1_id)
    , joint2_id(joint2_id)
    , joint1_placement(SE3::Identity())
    , joint2_placement(SE3::Identity())
    , reference_frame(reference_frame)
    , desired_contact_placement(SE3::Identity())
    , desired_contact_velocity(Motion::Zero())
    , desired_contact_acceleration(Motion::Zero())
    , corrector(size())
    , colwise_joint1_sparsity(model.nv)
    , colwise_joint2_sparsity(model.nv)
    , joint1_span_indexes((size_t)model.njoints)
    , joint2_span_indexes((size_t)model.njoints)
    , loop_span_indexes((size_t)model.nv)
    , compliance(Vector3::Zero())
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
    RigidConstraintModelTpl(
      const ContactType type,
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const ReferenceFrame & reference_frame = LOCAL)
    : Base(model)
    , type(type)
    , joint1_id(joint1_id)
    , joint2_id(0) // set to be the Universe
    , joint1_placement(SE3::Identity())
    , joint2_placement(SE3::Identity())
    , reference_frame(reference_frame)
    , desired_contact_placement(SE3::Identity())
    , desired_contact_velocity(Motion::Zero())
    , desired_contact_acceleration(Motion::Zero())
    , corrector(size())
    , colwise_joint1_sparsity(model.nv)
    , colwise_joint2_sparsity(model.nv)
    , joint1_span_indexes((size_t)model.njoints)
    , joint2_span_indexes((size_t)model.njoints)
    , loop_span_indexes((size_t)model.nv)
    , compliance(Vector3::Zero())
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
    const BooleanVector & getRowSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return colwise_sparsity;
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    const EigenIndexVector & getRowActiveIndexes(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return colwise_span_indexes;
    }

    ///
    ///  \brief Comparison operator
    ///
    /// \param[in] other Other RigidConstraintModelTpl to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributs
    /// must be the same).
    ///
    template<int OtherOptions>
    bool operator==(const RigidConstraintModelTpl<Scalar, OtherOptions> & other) const
    {
      return base() == other.base() && type == other.type && joint1_id == other.joint1_id
             && joint2_id == other.joint2_id && joint1_placement == other.joint1_placement
             && joint2_placement == other.joint2_placement
             && reference_frame == other.reference_frame && corrector == other.corrector
             && colwise_joint1_sparsity == other.colwise_joint1_sparsity
             && colwise_joint2_sparsity == other.colwise_joint2_sparsity
             && joint1_span_indexes == other.joint1_span_indexes
             && joint2_span_indexes == other.joint2_span_indexes && nv == other.nv
             && depth_joint1 == other.depth_joint1 && depth_joint2 == other.depth_joint2
             && colwise_sparsity == other.colwise_sparsity
             && colwise_span_indexes == other.colwise_span_indexes
             && loop_span_indexes == other.loop_span_indexes && compliance == other.compliance;
    }

    ///
    ///  \brief Oposite of the comparison operator.
    ///
    /// \param[in] other Other RigidConstraintModelTpl to compare with.
    ///
    /// \returns false if the two *this is not equal to other (at least type, joint1_id or placement
    /// attributs is different).
    ///
    template<int OtherOptions>
    bool operator!=(const RigidConstraintModelTpl<Scalar, OtherOptions> & other) const
    {
      return !(*this == other);
    }

    /// \brief Evaluate the constraint values at the current state given by data and store the
    /// results in cdata.
    template<template<typename, int> class JointCollectionTpl>
    void calc(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      RigidConstraintDataTpl<Scalar, Options> & cdata) const
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
    }

    /// \brief Returns the constraint projector associated with joint 1.
    /// This matrix transforms a spatial velocity expressed at the origin to the first component of
    /// the constraint associated with joint 1.
    template<ReferenceFrame rf>
    Matrix36
    getA1(const RigidConstraintDataTpl<Scalar, Options> & cdata, ReferenceFrameTag<rf>) const
    {
      Matrix36 res;
      typedef typename SE3::Vector3 Vector3;

      if (std::is_same<ReferenceFrameTag<rf>, WorldFrame>::value)
      {
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = oM1.rotation().transpose() * v_tmp;

        const SE3 & oM1 = cdata.oMc1;
        Vector3 v_tmp;
        res.template leftCols<3>() = oM1.rotation().transpose();
        INTERNAL_LOOP(0, oM1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, oM1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, oM1.translation(), res.template rightCols<3>());

#undef INTERNAL_LOOP
      }
      else if (std::is_same<ReferenceFrameTag<rf>, LocalFrame>::value)
      {
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = M1.rotation().transpose() * v_tmp;

        const SE3 & M1 = this->joint1_placement;
        Vector3 v_tmp;
        res.template leftCols<3>() = M1.rotation().transpose();
        INTERNAL_LOOP(0, M1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, M1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, M1.translation(), res.template rightCols<3>());

#undef INTERNAL_LOOP
      }

      return res;
    }

    /// \brief Returns the constraint projector associated with joint 2.
    /// This matrix transforms a spatial velocity expressed at the origin to the first component of
    /// the constraint associated with joint 2.
    template<ReferenceFrame rf>
    Matrix36
    getA2(const RigidConstraintDataTpl<Scalar, Options> & cdata, ReferenceFrameTag<rf>) const
    {
      Matrix36 res;
      typedef typename SE3::Vector3 Vector3;

      if (std::is_same<ReferenceFrameTag<rf>, WorldFrame>::value)
      {
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = oM1.rotation().transpose() * v_tmp;

        const SE3 & oM1 = cdata.oMc1;
        const SE3 & oM2 = cdata.oMc2;
        res.template leftCols<3>() = -oM1.rotation().transpose();
        Vector3 v_tmp;
        INTERNAL_LOOP(0, -oM2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, -oM2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, -oM2.translation(), res.template rightCols<3>());

#undef INTERNAL_LOOP
      }
      else if (std::is_same<ReferenceFrameTag<rf>, LocalFrame>::value)
      {
        const SE3 & j2Mc2 = this->joint2_placement;
        const SE3 & c1Mc2 = cdata.c1Mc2;
        const typename SE3::Matrix3 c1Rj2 = c1Mc2.rotation() * j2Mc2.rotation().transpose();
        res.template leftCols<3>() = -c1Rj2;
        Vector3 v_tmp;
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = -c1Rj2 * v_tmp;

        INTERNAL_LOOP(0, j2Mc2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, j2Mc2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, j2Mc2.translation(), res.template rightCols<3>());

#undef INTERNAL_LOOP
      }

      return res;
    }

    ///
    /// @brief This function computes the spatial inertia associated with the constraint.
    /// This function is useful to express the constraint inertia associated with the constraint for
    /// AL settings.
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
      template<typename, int>
      class JointCollectionTpl,
      typename Vector3Like,
      typename Matrix6Like,
      typename Matrix6LikeAllocator>
    void appendConstraintDiagonalInertiaToJointInertias(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const RigidConstraintDataTpl<Scalar, Options> & cdata,
      const Eigen::MatrixBase<Vector3Like> & diagonal_constraint_inertia,
      std::vector<Matrix6Like, Matrix6LikeAllocator> & inertias) const
    {
      EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Vector3Like, Vector3);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      PINOCCHIO_CHECK_ARGUMENT_SIZE(inertias.size(), size_t(model.njoints));
      assert(
        ((joint1_id > 0 && joint2_id == 0) || (joint1_id == 0 && joint2_id > 0))
        && "The behavior is only defined for this context");

      if (this->joint1_id != 0)
      {
        const SE3 & placement = this->joint1_placement;
        inertias[this->joint1_id] +=
          computeConstraintSpatialInertia(placement, diagonal_constraint_inertia);
      }

      if (this->joint2_id != 0)
      {
        const SE3 & placement = this->joint2_placement;
        inertias[this->joint2_id] +=
          computeConstraintSpatialInertia(placement, diagonal_constraint_inertia);
      }
    }

    template<
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int>
      class JointCollectionTpl>
    void jacobian_matrix_product(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      RigidConstraintDataTpl<Scalar, Options> & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res) const
    {
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef typename Data::Vector3 Vector3;
      OutputMatrix & res = _res.const_cast_derived();

      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), model.nv);
      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.cols(), res.cols());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), size());
      res.setZero();

      //      const Eigen::DenseIndex constraint_dim = size();
      //
      //      const Eigen::DenseIndex
      //      complexity_strategy_1 = 6 * res.cols() * 36 + constraint_dim * 36 * res.cols(),
      //      complexity_strategy_2 = 36 * constraint_dim * 6 + constraint_dim * 36 * res.cols();

      const Matrix36 A1 = getA1(cdata, WorldFrame());
      const Matrix36 A2 = getA2(cdata, WorldFrame());
      for (Eigen::DenseIndex jj = 0; jj < model.nv; ++jj)
      {
        if (!(colwise_joint1_sparsity[jj] || colwise_joint2_sparsity[jj]))
          continue;
        Matrix36 A;
        Vector3 AxSi;

        typedef typename Data::Matrix6x::ConstColXpr ConstColXpr;
        const ConstColXpr Jcol = data.J.col(jj);

        if (colwise_joint1_sparsity[jj] && colwise_joint2_sparsity[jj])
        {
          A = A1 + A2;
          AxSi.noalias() = A * Jcol;
        }
        else if (colwise_joint1_sparsity[jj])
          AxSi.noalias() = A1 * Jcol;
        else
          AxSi.noalias() = A2 * Jcol;

        res.noalias() += AxSi * mat.row(jj);
      }
    }

    ///  \brief Evaluate the Jacobian associated to the constraint at the given state stored in data
    /// and cdata.  The results Jacobian is evaluated in the jacobian input/output matrix.
    template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobian(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      RigidConstraintDataTpl<Scalar, Options> & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const
    {
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      JacobianMatrix & jacobian_matrix = _jacobian_matrix.const_cast_derived();

      const RigidConstraintModelTpl & cmodel = *this;

      const SE3 & oMc1 = cdata.oMc1;
      const SE3 & oMc2 = cdata.oMc2;
      const SE3 & c1Mc2 = cdata.c1Mc2;

      for (Eigen::DenseIndex jj = 0; jj < model.nv; ++jj)
      {

        if (colwise_joint1_sparsity[jj] || colwise_joint2_sparsity[jj])
        {
          const int sign = colwise_joint1_sparsity[jj] != colwise_joint2_sparsity[jj]
                             ? colwise_joint1_sparsity[jj] ? +1 : -1
                             : 0; // specific case for CONTACT_3D

          typedef typename Data::Matrix6x::ConstColXpr ConstColXpr;
          const ConstColXpr Jcol = data.J.col(jj);
          const MotionRef<const ConstColXpr> Jcol_motion(Jcol);

          switch (cmodel.type)
          {
          case CONTACT_3D: {
            switch (cmodel.reference_frame)
            {
            case LOCAL: {
              if (sign == 0)
              {
                const Motion Jcol_local1(oMc1.actInv(Jcol_motion)); // TODO: simplify computations
                const Motion Jcol_local2(oMc2.actInv(Jcol_motion)); // TODO: simplify computations
                const typename Motion::Vector3 Jdiff_linear =
                  Jcol_local1.linear() - c1Mc2.rotation() * Jcol_local2.linear();
                jacobian_matrix.col(jj) = Jdiff_linear;
                break;
              }
              else if (sign == 1)
              {
                const Motion Jcol_local(oMc1.actInv(Jcol_motion));
                jacobian_matrix.col(jj) = Jcol_local.linear();
                break;
              }
              else // sign == -1
              {
                Motion Jcol_local(oMc2.actInv(Jcol_motion)); // TODO: simplify computations
                Jcol_local.linear() =
                  c1Mc2.rotation() * Jcol_local.linear(); // TODO: simplify computations
                jacobian_matrix.col(jj) = -Jcol_local.linear();
                break;
              }
            }
            case LOCAL_WORLD_ALIGNED: {
              if (sign == 0)
              {
                const typename Motion::Vector3 Jdiff_linear =
                  (oMc2.translation() - oMc1.translation()).cross(Jcol_motion.angular());
                jacobian_matrix.col(jj) = Jdiff_linear;
                break;
              }
              else
              {
                typename Motion::Vector3 Jcol_local_world_aligned_linear(Jcol_motion.linear());
                if (sign == 1)
                  Jcol_local_world_aligned_linear -=
                    oMc1.translation().cross(Jcol_motion.angular());
                else
                  Jcol_local_world_aligned_linear -=
                    oMc2.translation().cross(Jcol_motion.angular());
                jacobian_matrix.col(jj) = Jcol_local_world_aligned_linear * Scalar(sign);
                break;
              }
            }
            case WORLD: {
              PINOCCHIO_THROW_PRETTY(
                std::invalid_argument, "Contact3D in world frame is not managed");
            }
            }
            break;
          }

          case CONTACT_6D: {
            switch (cmodel.reference_frame)
            {
            case LOCAL: {
              const Motion Jcol_local(oMc1.actInv(Jcol_motion));
              jacobian_matrix.col(jj) = Jcol_local.toVector() * Scalar(sign);
              break;
            }
            case LOCAL_WORLD_ALIGNED: {
              Motion Jcol_local_world_aligned(Jcol_motion);
              Jcol_local_world_aligned.linear() -=
                oMc1.translation().cross(Jcol_local_world_aligned.angular());
              jacobian_matrix.col(jj) = Jcol_local_world_aligned.toVector() * Scalar(sign);
              break;
            }
            case WORLD: {
              PINOCCHIO_THROW_PRETTY(
                std::invalid_argument, "Contact6D in world frame is not managed");
            }
            }
            break;
          }

          default:
            assert(false && "must never happened");
            break;
          }
        }
      }
    }

    /// \brief Map the constraint forces (aka constraint Lagrange multipliers) to the forces
    /// supported by the joints.
    template<
      template<typename, int>
      class JointCollectionTpl,
      typename ForceLike,
      typename ForceAllocator>
    void mapConstraintForceToJointForces(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const RigidConstraintDataTpl<Scalar, Options> & cdata,
      const Eigen::MatrixBase<ForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_forces.size(), size_t(model.njoints));
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), size());
      PINOCCHIO_UNUSED_VARIABLE(data);

      assert(this->type == CONTACT_3D);

      // Todo: optimize code
      const Matrix36 A1 = getA1(cdata, LocalFrame()), A2 = getA2(cdata, LocalFrame());
      joint_forces[this->joint1_id].toVector().noalias() += A1.transpose() * constraint_forces;
      joint_forces[this->joint2_id].toVector().noalias() += A2.transpose() * constraint_forces;
    }

    /// \brief Map the joint accelerations to constraint value
    template<
      template<typename, int>
      class JointCollectionTpl,
      typename MotionAllocator,
      typename VectorLike>
    void mapJointMotionsToConstraintMotion(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const RigidConstraintDataTpl<Scalar, Options> & cdata,
      const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_accelerations,
      const Eigen::MatrixBase<VectorLike> & constraint_value) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_accelerations.size(), size_t(model.njoints));
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_value.rows(), size());
      PINOCCHIO_UNUSED_VARIABLE(data);

      assert(this->type == CONTACT_3D);

      // Todo: optimize code

      if (this->joint1_id != 0 && this->joint2_id != 0)
      {
        const Matrix36 A1 = getA1(cdata, LocalFrame()), A2 = getA2(cdata, LocalFrame());
        constraint_value.const_cast_derived().noalias() =
          A1 * joint_accelerations[this->joint1_id].toVector()
          + A2 * joint_accelerations[this->joint2_id].toVector();
      }
      else if (this->joint1_id != 0)
      {
        const Matrix36 A1 = getA1(cdata, LocalFrame());
        constraint_value.const_cast_derived().noalias() =
          A1 * joint_accelerations[this->joint1_id].toVector();
      }
      else if (this->joint2_id != 0)
      {
        const Matrix36 A2 = getA2(cdata, LocalFrame());
        constraint_value.const_cast_derived().noalias() =
          A2 * joint_accelerations[this->joint2_id].toVector();
      }
      else
        constraint_value.const_cast_derived().setZero();
    }

    int size() const
    {
      switch (type)
      {
      case CONTACT_3D:
        return contact_dim<CONTACT_3D>::value;
      case CONTACT_6D:
        return contact_dim<CONTACT_6D>::value;
      default:
        return contact_dim<CONTACT_UNDEFINED>::value;
      }
      return -1;
    }

    /// \returns An expression of *this with the Scalar type casted to NewScalar.
    template<typename NewScalar>
    RigidConstraintModelTpl<NewScalar, Options> cast() const
    {
      typedef RigidConstraintModelTpl<NewScalar, Options> ReturnType;
      ReturnType res;
      res.base() = base();
      res.type = type;
      res.joint1_id = joint1_id;
      res.joint2_id = joint2_id;
      res.joint1_placement = joint1_placement.template cast<NewScalar>();
      res.joint2_placement = joint2_placement.template cast<NewScalar>();
      res.reference_frame = reference_frame;
      res.desired_contact_placement = desired_contact_placement.template cast<NewScalar>();
      res.desired_contact_velocity = desired_contact_velocity.template cast<NewScalar>();
      res.desired_contact_acceleration = desired_contact_acceleration.template cast<NewScalar>();
      res.corrector = corrector.template cast<NewScalar>();
      res.colwise_joint1_sparsity = colwise_joint1_sparsity;
      res.colwise_joint2_sparsity = colwise_joint2_sparsity;
      res.nv = nv;
      res.depth_joint1 = depth_joint1;
      res.depth_joint2 = depth_joint2;
      res.loop_span_indexes = loop_span_indexes;

      return res;
    }

  protected:
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void init(const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        reference_frame == LOCAL || reference_frame == LOCAL_WORLD_ALIGNED,
        "reference_frame should be LOCAL or LOCAL_WORLD_ALIGNED");
      nv = model.nv;
      depth_joint1 = static_cast<size_t>(model.supports[joint1_id].size());
      depth_joint2 = static_cast<size_t>(model.supports[joint2_id].size());

      typedef ModelTpl<Scalar, OtherOptions, JointCollectionTpl> Model;
      typedef typename Model::JointModel JointModel;
      static const bool default_sparsity_value = false;
      colwise_joint1_sparsity.fill(default_sparsity_value);
      colwise_joint2_sparsity.fill(default_sparsity_value);
      joint1_span_indexes.reserve((size_t)model.njoints);
      joint2_span_indexes.reserve((size_t)model.njoints);

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

      std::rotate(
        joint1_span_indexes.rbegin(), joint1_span_indexes.rbegin() + 1, joint1_span_indexes.rend());
      std::rotate(
        joint2_span_indexes.rbegin(), joint2_span_indexes.rbegin() + 1, joint2_span_indexes.rend());
      colwise_span_indexes.reserve((size_t)model.nv);
      colwise_sparsity.resize(model.nv);
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
    }
  };

  ///
  /// \brief Computes the sum of the sizes of the constraints contained in the input
  /// `contact_models` vector.
  template<template<typename T> class Holder, class ConstraintModel, class ConstraintModelAllocator>
  Eigen::DenseIndex getTotalConstraintSize(
    const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & constraint_models)
  {
    Eigen::DenseIndex total_size = 0;
    for (const auto & wrapper : constraint_models)
    {
      const ConstraintModel & constraint_model = wrapper;
      total_size += constraint_model.size();
    }

    return total_size;
  }

  ///
  /// \brief Computes the sum of the sizes of the constraints contained in the input
  /// `contact_models` vector.
  template<class ConstraintModel, class ConstraintModelAllocator>
  Eigen::DenseIndex getTotalConstraintSize(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & contact_models)
  {
    typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
    typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

    WrappedConstraintModelVector wrapped_constraint_models(
      contact_models.cbegin(), contact_models.cend());

    return getTotalConstraintSize(wrapped_constraint_models);
  }

  ///
  ///  \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct RigidConstraintDataTpl : ConstraintDataBase<RigidConstraintDataTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef RigidConstraintModelTpl<Scalar, Options> ContactModel;
    typedef RigidConstraintDataTpl ContactData;

    typedef SE3Tpl<Scalar, Options> SE3;
    typedef MotionTpl<Scalar, Options> Motion;
    typedef ForceTpl<Scalar, Options> Force;
    typedef Eigen::Matrix<Scalar, 6, 6, Options> Matrix6;
    typedef PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(Matrix6) VectorOfMatrix6;
    typedef Eigen::Matrix<Scalar, 6, Eigen::Dynamic, Options> Matrix6x;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> MatrixX;

    // data

    /// \brief Resulting contact forces
    Force contact_force;

    /// \brief Placement of the constraint frame 1 with respect to the WORLD frame
    SE3 oMc1;

    /// \brief Placement of the constraint frame 2 with respect to the WORLD frame
    SE3 oMc2;

    /// \brief Relative displacement between the two frames
    SE3 c1Mc2;

    /// \brief Current contact placement error
    Motion contact_placement_error;

    /// \brief Current contact spatial velocity of the constraint 1
    Motion contact1_velocity;

    /// \brief Current contact spatial velocity of the constraint 2
    Motion contact2_velocity;

    /// \brief Current contact velocity error
    Motion contact_velocity_error;

    /// \brief Current contact spatial acceleration
    Motion contact_acceleration;

    /// \brief Contact spatial acceleration desired
    Motion contact_acceleration_desired;

    /// \brief Current contact spatial error (due to the integration step).
    Motion contact_acceleration_error;

    /// \brief Current contact drift acceleration (acceleration only due to the Coriolis and
    /// centrifugal effects) for the constraint frame 1.
    Motion contact1_acceleration_drift;

    /// \brief Current contact drift acceleration (acceleration only due to the Coriolis and
    /// centrifugal effects) for the constraint frame 2.
    Motion contact2_acceleration_drift;

    /// \brief Contact deviation from the reference acceleration (a.k.a the error)
    Motion contact_acceleration_deviation;

    VectorOfMatrix6 extended_motion_propagators_joint1;
    VectorOfMatrix6 lambdas_joint1;
    VectorOfMatrix6 extended_motion_propagators_joint2;

    Matrix6x dv1_dq, da1_dq, da1_dv, da1_da;
    Matrix6x dv2_dq, da2_dq, da2_dv, da2_da;
    MatrixX dvc_dq, dac_dq, dac_dv, dac_da;

    /// \brief Default constructor
    RigidConstraintDataTpl()
    {
    }

    explicit RigidConstraintDataTpl(const ContactModel & contact_model)
    : contact_force(Force::Zero())
    , contact_placement_error(Motion::Zero())
    , contact1_velocity(Motion::Zero())
    , contact2_velocity(Motion::Zero())
    , contact_velocity_error(Motion::Zero())
    , contact_acceleration(Motion::Zero())
    , contact_acceleration_desired(Motion::Zero())
    , contact_acceleration_error(Motion::Zero())
    , contact1_acceleration_drift(Motion::Zero())
    , contact2_acceleration_drift(Motion::Zero())
    , contact_acceleration_deviation(Motion::Zero())
    , extended_motion_propagators_joint1(contact_model.depth_joint1, Matrix6::Zero())
    , lambdas_joint1(contact_model.depth_joint1, Matrix6::Zero())
    , extended_motion_propagators_joint2(contact_model.depth_joint2, Matrix6::Zero())
    , dv1_dq(6, contact_model.nv)
    , da1_dq(6, contact_model.nv)
    , da1_dv(6, contact_model.nv)
    , da1_da(6, contact_model.nv)
    , dv2_dq(6, contact_model.nv)
    , da2_dq(6, contact_model.nv)
    , da2_dv(6, contact_model.nv)
    , da2_da(6, contact_model.nv)
    , dvc_dq(contact_model.size(), contact_model.nv)
    , dac_dq(contact_model.size(), contact_model.nv)
    , dac_dv(contact_model.size(), contact_model.nv)
    , dac_da(contact_model.size(), contact_model.nv)
    {
    }

    bool operator==(const RigidConstraintDataTpl & other) const
    {
      return contact_force == other.contact_force && oMc1 == other.oMc1 && oMc2 == other.oMc2
             && c1Mc2 == other.c1Mc2 && contact_placement_error == other.contact_placement_error
             && contact1_velocity == other.contact1_velocity
             && contact2_velocity == other.contact2_velocity
             && contact_velocity_error == other.contact_velocity_error
             && contact_acceleration == other.contact_acceleration
             && contact_acceleration_desired == other.contact_acceleration_desired
             && contact_acceleration_error == other.contact_acceleration_error
             && contact1_acceleration_drift == other.contact1_acceleration_drift
             && contact2_acceleration_drift == other.contact2_acceleration_drift
             && contact_acceleration_deviation == other.contact_acceleration_deviation
             && extended_motion_propagators_joint1 == other.extended_motion_propagators_joint1
             && lambdas_joint1 == other.lambdas_joint1
             && extended_motion_propagators_joint2 == other.extended_motion_propagators_joint2
             //
             && dv1_dq == other.dv1_dq && da1_dq == other.da1_dq && da1_dv == other.da1_dv
             && da1_da == other.da1_da
             //
             && dv2_dq == other.dv2_dq && da2_dq == other.da2_dq && da2_dv == other.da2_dv
             && da2_da == other.da2_da
             //
             && dvc_dq == other.dvc_dq && dac_dq == other.dac_dq && dac_dv == other.dac_dv
             && dac_da == other.dac_da;
    }

    bool operator!=(const RigidConstraintDataTpl & other) const
    {
      return !(*this == other);
    }
  };
  PINOCCHIO_COMPILER_DIAGNOSTIC_POP

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_contact_info_hpp__
