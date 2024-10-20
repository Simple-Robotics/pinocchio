//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_unbounded_set_hpp__
#define __pinocchio_algorithm_constraints_unbounded_set_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/set-base.hpp"

namespace pinocchio
{

  template<typename _Scalar, int _Options>
  struct traits<UnboundedSetTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
  };

  ///  \brief Unbounded set covering the whole space
  template<typename _Scalar, int _Options>
  struct UnboundedSetTpl : SetBase<UnboundedSetTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;

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

    /// \brief Project a vector x into the box.
    ///
    /// \param[in] x a vector to project.
    ///
    template<typename VectorLike>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike)
      project(const Eigen::MatrixBase<VectorLike> & x) const
    {
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike) VectorPlain;
      return VectorPlain(x); // make a copy
    }

    /// \brief Returns the dimension of the ambiant space.
    Eigen::DenseIndex dim() const
    {
      return m_size;
    }

  protected:
    Eigen::DenseIndex m_size;
  }; // UnboundedSetTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_unbounded_set_hpp__
