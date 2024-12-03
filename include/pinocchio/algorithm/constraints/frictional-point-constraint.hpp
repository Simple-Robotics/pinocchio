//
// Copyright (c) 2019-2024 INRIA CNRS
//

#ifndef __pinocchio_algorithm_constraints_frictional_point_constraint_hpp__
#define __pinocchio_algorithm_constraints_frictional_point_constraint_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/coulomb-friction-cone.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-data-base.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, FrictionalPointConstraintModelTpl<Scalar, Options>>
  {
    typedef FrictionalPointConstraintModelTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<FrictionalPointConstraintModelTpl<_Scalar, _Options>>
  : traits<PointConstraintModelBase<FrictionalPointConstraintModelTpl<_Scalar, _Options>>>
  {
    typedef _Scalar Scalar;

    enum
    {
      Options = _Options
    };

    typedef FrictionalPointConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef FrictionalPointConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef CoulombFrictionConeTpl<Scalar> ConstraintSet;

    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef Vector3 VectorConstraintSize;

    typedef Vector3 ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;
  };

  template<typename _Scalar, int _Options>
  struct traits<FrictionalPointConstraintDataTpl<_Scalar, _Options>>
  : traits<FrictionalPointConstraintModelTpl<_Scalar, _Options>>
  {
  };

  ///
  ///  \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct FrictionalPointConstraintModelTpl
  : PointConstraintModelBase<FrictionalPointConstraintModelTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef PointConstraintModelBase<FrictionalPointConstraintModelTpl> Base;

    template<typename NewScalar, int NewOptions>
    friend struct FrictionalPointConstraintModelTpl;

    typedef FrictionalPointConstraintModelTpl ConstraintModel;
    typedef FrictionalPointConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef CoulombFrictionConeTpl<Scalar> ConstraintSet;

    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;
    using typename Base::BooleanVector;
    using typename Base::ComplianceVectorType;
    using typename Base::EigenIndexVector;
    using typename Base::Force;
    using typename Base::Matrix36;
    using typename Base::Matrix6;
    using typename Base::Motion;
    using typename Base::SE3;
    using typename Base::Vector3;
    using typename Base::Vector6;
    using typename Base::VectorConstraintSize;

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

  protected:
    ///
    ///  \brief Default constructor
    ///
    FrictionalPointConstraintModelTpl()
    : Base()
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
    /// expressed.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    FrictionalPointConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement,
      const JointIndex joint2_id,
      const SE3 & joint2_placement)
    : Base(model, joint1_id, joint1_placement, joint2_id, joint2_placement)
    {
    }

    ///
    ///  \brief Contructor with from a given type, joint1_id and placement.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint1_placement Placement of the constraint w.r.t the frame of joint1.
    /// expressed.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    FrictionalPointConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement)
    : Base(model, joint1_id, joint1_placement)
    {
    }

    ///
    ///  \brief Contructor with from a given type and the joint ids.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint2_id Index of the joint 2 in the model tree.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    FrictionalPointConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const JointIndex joint2_id)
    : Base(model, joint1_id, joint2_id)
    {
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
    FrictionalPointConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model, const JointIndex joint1_id)
    : Base(model, joint1_id)
    {
    }

    ///
    /// \brief Create data storage associated to the constraint
    ///
    ConstraintData createData() const
    {
      return ConstraintData(*this);
    }

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, FrictionalPointConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, FrictionalPointConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::template cast<NewScalar>(res);
      res.m_set = m_set.template cast<NewScalar>();
      return res;
    }

    ///
    ///  \brief Comparison operator
    ///
    /// \param[in] other Other FrictionalPointConstraintModelTpl to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributs
    /// must be the same).
    ///
    bool operator==(const FrictionalPointConstraintModelTpl & other) const
    {
      return base() == other.base();
    }

    ///
    ///  \brief Oposite of the comparison operator.
    ///
    /// \param[in] other Other FrictionalPointConstraintModelTpl to compare with.
    ///
    /// \returns false if the two *this is not equal to other (at least type, joint1_id or placement
    /// attributs is different).
    ///
    bool operator!=(const FrictionalPointConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    const ConstraintSet & set() const
    {
      return m_set;
    }

    ConstraintSet & set()
    {
      return m_set;
    }

    using Base::compliance;

  protected:
    ConstraintSet m_set = ConstraintSet();

  }; // struct FrictionalPointConstraintModelTpl<_Scalar,_Options>

  ///
  ///  \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct FrictionalPointConstraintDataTpl
  : PointConstraintDataBase<FrictionalPointConstraintDataTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef FrictionalPointConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef FrictionalPointConstraintDataTpl ConstraintData;
    typedef PointConstraintDataBase<FrictionalPointConstraintDataTpl> Base;

    using typename Base::Force;
    using typename Base::Matrix6;
    using typename Base::Matrix6x;
    using typename Base::MatrixX;
    using typename Base::Motion;
    using typename Base::SE3;
    using typename Base::Vector3;
    using typename Base::VectorOfMatrix6;

    /// \brief Default constructor
    FrictionalPointConstraintDataTpl()
    {
    }

    explicit FrictionalPointConstraintDataTpl(const ConstraintModel & constraint_model)
    : Base(constraint_model)
    {
    }

    bool operator==(const FrictionalPointConstraintDataTpl & other) const
    {
      return base() == other.base();
    }

    bool operator!=(const FrictionalPointConstraintDataTpl & other) const
    {
      return !(*this == other);
    }

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

  }; // struct FrictionalPointConstraintDataTpl<_Scalar,_Options>

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_frictional_point_constraint_hpp__
