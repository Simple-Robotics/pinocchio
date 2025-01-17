//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_unbounded_set_hpp__
#define __pinocchio_algorithm_constraints_unbounded_set_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/cone-base.hpp"

namespace pinocchio
{

  template<typename Scalar, int _Options>
  struct NullSetTpl;

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, UnboundedSetTpl<Scalar, Options>>
  {
    typedef UnboundedSetTpl<NewScalar, Options> type;
  };

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, NullSetTpl<Scalar, Options>>
  {
    typedef NullSetTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<UnboundedSetTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef NullSetTpl<Scalar, _Options> DualCone;
  };

  template<typename _Scalar, int _Options>
  struct traits<NullSetTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;

    enum
    {
      Options = _Options
    };
    typedef UnboundedSetTpl<Scalar, _Options> DualCone;
  };

  ///  \brief Unbounded set covering the whole space
  template<typename _Scalar, int _Options>
  struct UnboundedSetTpl : ConeBase<UnboundedSetTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef SetBase<UnboundedSetTpl> Base;

    /// \brief Constructor from a given size
    ///
    explicit UnboundedSetTpl(const Eigen::DenseIndex size)
    : m_size(size)
    {
    }

    /// \brief Copy constructor.
    UnboundedSetTpl(const UnboundedSetTpl & other) = default;

    /// \brief Copy operator
    UnboundedSetTpl & operator=(const UnboundedSetTpl & other) = default;

    /// \brief Cast operator
    template<typename NewScalar>
    UnboundedSetTpl<NewScalar, Options> cast() const
    {
      typedef UnboundedSetTpl<NewScalar, Options> ReturnType;
      return ReturnType(this->size());
    }

    /// \brief Comparison operator
    bool operator==(const UnboundedSetTpl & other) const
    {
      return m_size == other.m_size;
    }

    /// \brief Difference  operator
    bool operator!=(const UnboundedSetTpl & other) const
    {
      return !(*this == other);
    }

    /// \brief Check whether a vector x lies within the box.
    ///
    /// \param[in] f vector to check (assimilated to a  force vector).
    ///
    template<typename VectorLike>
    bool isInside(const Eigen::MatrixBase<VectorLike> & x, const Scalar prec = Scalar(0)) const
    {
      assert(prec >= 0 && "prec should be positive");
      PINOCCHIO_UNUSED_VARIABLE(x);
      PINOCCHIO_UNUSED_VARIABLE(prec);
      return true;
    }

    using Base::project;

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
      res_.const_cast_derived() = x;
    }

    /// \brief Returns the dimension of the ambiant space.
    Eigen::DenseIndex dim() const
    {
      return m_size;
    }

    Eigen::DenseIndex size() const
    {
      return m_size;
    }

  protected:
    Eigen::DenseIndex m_size;
  }; // UnboundedSetTpl

  ///  \brief Null set containing 0.
  template<typename _Scalar, int _Options>
  struct NullSetTpl : ConeBase<NullSetTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef SetBase<NullSetTpl> Base;

    /// \brief Constructor from a given size
    ///
    explicit NullSetTpl(const Eigen::DenseIndex size)
    : m_size(size)
    {
    }

    /// \brief Copy constructor.
    NullSetTpl(const NullSetTpl & other) = default;

    /// \brief Copy operator
    NullSetTpl & operator=(const NullSetTpl & other) = default;

    /// \brief Cast operator
    template<typename NewScalar>
    NullSetTpl<NewScalar, Options> cast() const
    {
      typedef NullSetTpl<NewScalar, Options> ReturnType;
      return ReturnType(this->size());
    }

    /// \brief Comparison operator
    bool operator==(const NullSetTpl & other) const
    {
      return m_size == other.m_size;
    }

    /// \brief Difference  operator
    bool operator!=(const NullSetTpl & other) const
    {
      return !(*this == other);
    }

    /// \brief Check whether a vector x is zero.
    ///
    /// \param[in] f vector to check (assimilated to a  force vector).
    ///
    template<typename VectorLike>
    bool isInside(const Eigen::MatrixBase<VectorLike> & x, const Scalar prec = Scalar(0)) const
    {
      assert(prec >= 0 && "prec should be positive");
      return x.isZero(prec);
    }

    using Base::project;

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
      PINOCCHIO_UNUSED_VARIABLE(x);
      auto & res = res_.const_cast_derived();
      res.setZero();
    }

    /// \brief Returns the dimension of the ambiant space.
    Eigen::DenseIndex dim() const
    {
      return m_size;
    }

    Eigen::DenseIndex size() const
    {
      return m_size;
    }

  protected:
    Eigen::DenseIndex m_size;
  }; // NullSetTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_unbounded_set_hpp__
