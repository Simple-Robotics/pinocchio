//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_box_set_hpp__
#define __pinocchio_algorithm_constraints_box_set_hpp__

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/set-base.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, BoxSetTpl<Scalar, Options>>
  {
    typedef BoxSetTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<BoxSetTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
  };

  ///  \brief Box set defined by a lower and an upper bounds [lb;ub].
  template<typename _Scalar, int _Options>
  struct BoxSetTpl : SetBase<BoxSetTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef SetBase<BoxSetTpl> Base;

    /// \brief Default constructor
    ///
    BoxSetTpl()
    {
    }

    /// \brief Constructor from a given size
    ///
    explicit BoxSetTpl(const Eigen::DenseIndex size)
    : m_lb(Vector::Constant(size, -std::numeric_limits<Scalar>::infinity()))
    , m_ub(Vector::Constant(size, +std::numeric_limits<Scalar>::infinity()))
    {
    }

    template<typename V1, typename V2>
    BoxSetTpl(const Eigen::MatrixBase<V1> & lb, const Eigen::MatrixBase<V2> & ub)
    : m_lb(lb)
    , m_ub(ub)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        (m_lb.array() <= m_ub.array()).all(), "Some components of lb are greater than ub");
    }

    /// \brief Copy constructor.
    BoxSetTpl(const BoxSetTpl & other) = default;

    /// \brief Cast operator
    template<typename NewScalar>
    BoxSetTpl<NewScalar, Options> cast() const
    {
      typedef BoxSetTpl<NewScalar, Options> ReturnType;
      return ReturnType(
        this->m_lb.template cast<NewScalar>(), this->m_ub.template cast<NewScalar>());
    }

    /// \brief Copy operator
    BoxSetTpl & operator=(const BoxSetTpl & other) = default;

    /// \brief Comparison operator
    bool operator==(const BoxSetTpl & other) const
    {
      return m_lb == other.m_lb && m_ub == other.m_ub;
    }

    /// \brief Difference  operator
    bool operator!=(const BoxSetTpl & other) const
    {
      return !(*this == other);
    }

    /// \brief Resize by calling the resize method of Eigen.
    void resize(Eigen::DenseIndex new_size)
    {
      m_lb.resize(new_size);
      m_ub.resize(new_size);
    }

    /// \brief Resize by calling the conservativeResize method of Eigen.
    void conservativeResize(Eigen::DenseIndex new_size)
    {
      m_lb.conservativeResize(new_size);
      m_ub.conservativeResize(new_size);
    }

    /// \brief Check whether a vector x lies within the box.
    ///
    /// \param[in] f vector to check (assimilated to a  force vector).
    ///
    template<typename VectorLike>
    bool isInside(const Eigen::MatrixBase<VectorLike> & x, const Scalar prec = Scalar(0)) const
    {
      assert(prec >= 0 && "prec should be positive");
      return (x - project(x)).norm() <= prec;
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
      res_.const_cast_derived() = x.array().max(m_lb.array()).min(m_ub.array());
    }

    /// \brief Returns the dimension of the box.
    Eigen::DenseIndex dim() const
    {
      return m_lb.size();
    }

    Eigen::DenseIndex size() const
    {
      return m_lb.size();
    }

    const Vector & lb() const
    {
      return m_lb;
    }
    Vector & lb()
    {
      return m_lb;
    }

    const Vector & ub() const
    {
      return m_ub;
    }
    Vector & ub()
    {
      return m_ub;
    }

    /// \brief Check whether lb <= ub for all components
    bool isValid() const
    {
      (m_lb.array() <= m_ub.array).all();
    }

    /// \brief Project the value given as input for the given row index.
    Scalar rowiseProject(const Eigen::DenseIndex row_id, const Scalar value) const
    {
      assert(row_id < size());
      return math::max(m_lb[row_id], math::min(m_ub[row_id], value));
    }

  protected:
    Vector m_lb, m_ub;
  }; // BoxSetTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_box_set_hpp__
