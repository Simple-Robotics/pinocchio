//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_box_set_hpp__
#define __pinocchio_algorithm_constraints_box_set_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/set-base.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar>
  struct CastType<NewScalar, PositiveOrthantConeTpl<Scalar>>
  {
    typedef PositiveOrthantConeTpl<NewScalar> type;
  };

  template<typename _Scalar>
  struct traits<PositiveOrthantConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
  };

  template<typename NewScalar, typename Scalar>
  struct CastType<NewScalar, NegativeOrthantConeTpl<Scalar>>
  {
    typedef NegativeOrthantConeTpl<NewScalar> type;
  };

  template<typename _Scalar>
  struct traits<NegativeOrthantConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
  };

  ///  \brief Box set defined by a lower and an upper bounds [lb;ub].
  template<typename Derived>
  struct OrthantConeBase : SetBase<Derived>
  {
    typedef typename traits<Derived>::Scalar Scalar;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef SetBase<Derived> Base;

    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }
    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    using Base::derived;
    using Base::project;

    /// \brief Default constructor
    ///
    OrthantConeBase()
    : m_size(0)
    {
    }

    /// \brief Constructor from a given size
    ///
    explicit OrthantConeBase(const Eigen::DenseIndex size)
    : m_size(size)
    {
    }

    /// \brief Copy constructor.
    OrthantConeBase(const OrthantConeBase & other) = default;

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, Derived>::type cast() const
    {
      typedef typename CastType<NewScalar, Derived>::type ReturnType;
      return ReturnType(size());
    }

    /// \brief Copy operator
    OrthantConeBase & operator=(const OrthantConeBase & other) = default;

    /// \brief Comparison operator
    bool operator==(const OrthantConeBase & other) const
    {
      return m_size == other.m_size;
    }

    /// \brief Difference  operator
    bool operator!=(const OrthantConeBase & other) const
    {
      return !(*this == other);
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

    /// \brief Check whether a vector x lies within the orthant.
    ///
    /// \param[in] x vector to check .
    ///
    template<typename VectorLike>
    bool isInside(const Eigen::MatrixBase<VectorLike> & x, const Scalar prec = Scalar(0)) const
    {
      assert(prec >= 0 && "prec should be positive");
      return (x - project(x)).norm() <= prec;
    }

    /// \brief Returns the dual cone associated with this.
    ///
    /// \remarks Orthant cone are by definition self dual.
    const Derived & dual() const
    {
      return derived();
    }

    /// \brief Returns the dimension of the box.
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
  }; // OrthantConeBase

  template<typename _Scalar>
  struct PositiveOrthantConeTpl : OrthantConeBase<PositiveOrthantConeTpl<_Scalar>>
  {
    typedef OrthantConeBase<PositiveOrthantConeTpl> Base;

    typedef _Scalar Scalar;

    /// \brief Default constructor
    ///
    PositiveOrthantConeTpl() = default;

    /// \brief Constructor from a given size
    ///
    explicit PositiveOrthantConeTpl(const Eigen::DenseIndex size)
    : Base(size)
    {
    }

    using Base::project;
    using Base::operator==;
    using Base::operator!=;

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
      res_.const_cast_derived() = x.array().max(Scalar(0)).matrix();
    }

  }; // struct PositiveOrthantTpl

  ///  \brief Negative orthant
  template<typename _Scalar>
  struct NegativeOrthantConeTpl : OrthantConeBase<NegativeOrthantConeTpl<_Scalar>>
  {
    typedef OrthantConeBase<NegativeOrthantConeTpl> Base;

    typedef _Scalar Scalar;

    /// \brief Default constructor
    ///
    NegativeOrthantConeTpl() = default;

    /// \brief Constructor from a given size
    ///
    explicit NegativeOrthantConeTpl(const Eigen::DenseIndex size)
    : Base(size)
    {
    }

    using Base::project;
    using Base::operator==;
    using Base::operator!=;

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
      res_.const_cast_derived() = x.array().max(Scalar(0)).matrix();
    }

  }; // struct PositiveOrthantTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_box_set_hpp__
