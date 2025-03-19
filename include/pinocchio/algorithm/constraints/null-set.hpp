//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_null_set_hpp__
#define __pinocchio_algorithm_constraints_null_set_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/cone-base.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, NullSetTpl<Scalar, Options>>
  {
    typedef NullSetTpl<NewScalar, Options> type;
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

  ///  \brief Null set containing (0 singleton).
  template<typename _Scalar, int _Options>
  struct NullSetTpl : ConeBase<NullSetTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef ConeBase<NullSetTpl> Base;
    typedef typename traits<NullSetTpl>::DualCone DualCone;

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
      return base() == other.base() && m_size == other.m_size;
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

    /// \brief Project a vector x into set.
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

    /// \brief Resize by calling the resize method of Eigen.
    void resize(Eigen::DenseIndex new_size)
    {
      m_size = new_size;
    }

    /// \brief Resize by calling the conservativeResize method of Eigen.
    void conservativeResize(Eigen::DenseIndex new_size)
    {
      this->resize(new_size);
    }

    DualCone dual() const
    {
      return DualCone();
    }

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

  protected:
    Eigen::DenseIndex m_size;
  }; // NullSetTpl

} // namespace pinocchio

#endif // __pinocchio_algorithm_constraints_null_set_hpp__
