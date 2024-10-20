//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_box_set_hpp__
#define __pinocchio_algorithm_constraints_box_set_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/set-base.hpp"

namespace pinocchio
{

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

    /// \brief Project a vector x into the box.
    ///
    /// \param[in] x a vector to project.
    ///
    template<typename VectorLike>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike)
      project(const Eigen::MatrixBase<VectorLike> & x) const
    {
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike) VectorPlain;
      return VectorPlain(x.array().max(m_lb.array()).min(m_ub.array()));
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

  protected:
    Vector m_lb, m_ub;
  }; // BoxSetTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_box_set_hpp__
