//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_set_hpp__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_set_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/orthant-cone.hpp"

namespace pinocchio
{

  template<typename Scalar>
  struct JointLimitConstraintConeTpl;

  template<typename NewScalar, typename Scalar>
  struct CastType<NewScalar, JointLimitConstraintConeTpl<Scalar>>
  {
    typedef JointLimitConstraintConeTpl<NewScalar> type;
  };

  template<typename _Scalar>
  struct traits<JointLimitConstraintConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
    typedef JointLimitConstraintConeTpl<Scalar> DualCone;
  };

  template<typename _Scalar>
  struct JointLimitConstraintConeTpl : ConeBase<JointLimitConstraintConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
    typedef ConeBase<JointLimitConstraintConeTpl> Base;
    typedef typename traits<JointLimitConstraintConeTpl>::DualCone DualCone;

    typedef PositiveOrthantConeTpl<Scalar> PositiveOrthantCone;
    typedef NegativeOrthantConeTpl<Scalar> NegativeOrthantCone;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

    using Base::project;

    JointLimitConstraintConeTpl() = default;

    JointLimitConstraintConeTpl(
      const Eigen::DenseIndex negative_orthant_size, const Eigen::DenseIndex positive_orthant_size)
    : negative_orthant(negative_orthant_size)
    , positive_orthant(positive_orthant_size)
    {
    }

    ///  \brief Copy constructor
    JointLimitConstraintConeTpl(const JointLimitConstraintConeTpl & other)
    : negative_orthant(other.negative_orthant)
    , positive_orthant(other.positive_orthant)
    {
    }

    void resize(
      const Eigen::DenseIndex negative_orthant_size, const Eigen::DenseIndex positive_orthant_size)
    {
      negative_orthant.resize(negative_orthant_size);
      positive_orthant.resize(positive_orthant_size);
    }

    /// \brief Cast operator
    template<typename NewScalar>
    JointLimitConstraintConeTpl<NewScalar> cast() const
    {
      typedef JointLimitConstraintConeTpl<NewScalar> ReturnType;
      return ReturnType(negative_orthant.size(), positive_orthant.size());
    }

    /// \brief Returns the dual cone associated with this.
    ///
    /// \remarks Orthant cones are by definition self dual.
    DualCone dual() const
    {
      return JointLimitConstraintConeTpl(*this);
    }

    /// \brief Returns the dimension of the box.
    Eigen::DenseIndex dim() const
    {
      return negative_orthant.size() + positive_orthant.size();
    }

    Eigen::DenseIndex size() const
    {
      return dim();
    }

    /// \brief Comparison operator
    bool operator==(const JointLimitConstraintConeTpl & other) const
    {
      return negative_orthant == other.negative_orthant
        && positive_orthant == other.positive_orthant;
    }

    /// \brief Difference  operator
    bool operator!=(const JointLimitConstraintConeTpl & other) const
    {
      return !(*this == other);
    }

    /// \brief Check whether a vector x lies within the orthant.
    ///
    /// \param[in] x vector to check .
    ///
    template<typename VectorLike>
    bool isInside(const Eigen::MatrixBase<VectorLike> & x, const Scalar prec = Scalar(0)) const
    {
      assert(x.size() == size());
      return negative_orthant.isInsidex(x.head(negative_orthant.size()), prec)
             && positive_orthant.isInsidex(x.tail(positive_orthant.size()), prec);
    }

    /// \brief Project a vector x into orthant.
    ///
    /// \param[in] x a vector to project.
    /// \param[in] res result of the projection.
    ///
    template<typename VectorLikeIn, typename VectorLikeOut>
    void project(
      const Eigen::MatrixBase<VectorLikeIn> & x,
      const Eigen::MatrixBase<VectorLikeOut> & res_) const
    {
      auto & res = res_.const_cast_derived();
      negative_orthant.project(x.head(negative_orthant.size()), res.head(negative_orthant.size()));
      positive_orthant.project(x.tail(positive_orthant.size()), res.tail(positive_orthant.size()));
    }

    /// \brief Project the value given as input for the given row index.
    Scalar rowiseProject(const Eigen::DenseIndex row_id, const Scalar value) const
    {
      assert(row_id < size());
      if (row_id < negative_orthant.size())
      {
        return negative_orthant.rowiseProject(row_id, value);
      }
      else
      {
        return positive_orthant.rowiseProject(row_id - negative_orthant.size(), value);
      }
    }

    /// \brief Returns a const reference to the negative orthant.
    ///
    const NegativeOrthantCone & getNegativeOrthant() const
    {
      return negative_orthant;
    }

    /// \brief Returns a const reference to the positive orthant.
    ///
    const PositiveOrthantCone & getPositiveOrthant() const
    {
      return positive_orthant;
    }

  protected:
    NegativeOrthantCone negative_orthant;
    PositiveOrthantCone positive_orthant;
  };
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_set_hpp__
