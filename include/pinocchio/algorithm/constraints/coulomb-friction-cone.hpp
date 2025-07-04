//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_coulomb_friction_cone_hpp__
#define __pinocchio_algorithm_constraints_coulomb_friction_cone_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/cone-base.hpp"
#include "pinocchio/math/fwd.hpp"
#include "pinocchio/math/comparison-operators.hpp"

namespace pinocchio
{

  template<typename Scalar>
  struct DualCoulombFrictionConeTpl;

  template<typename NewScalar, typename Scalar>
  struct CastType<NewScalar, CoulombFrictionConeTpl<Scalar>>
  {
    typedef CoulombFrictionConeTpl<NewScalar> type;
  };

  template<typename NewScalar, typename Scalar>
  struct CastType<NewScalar, DualCoulombFrictionConeTpl<Scalar>>
  {
    typedef DualCoulombFrictionConeTpl<NewScalar> type;
  };

  template<typename _Scalar>
  struct traits<CoulombFrictionConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
    typedef DualCoulombFrictionConeTpl<Scalar> DualCone;
  };

  template<typename _Scalar>
  struct traits<DualCoulombFrictionConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
    typedef CoulombFrictionConeTpl<Scalar> DualCone;
  };

  ///  \brief 3d Coulomb friction cone.
  template<typename _Scalar>
  struct CoulombFrictionConeTpl : ConeBase<CoulombFrictionConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
    typedef typename traits<CoulombFrictionConeTpl>::DualCone DualCone;
    typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
    typedef ConeBase<CoulombFrictionConeTpl> Base;

    ///
    /// \brief Default constructor
    ///
    CoulombFrictionConeTpl()
    : mu(+std::numeric_limits<Scalar>::infinity())
    {
    }

    ///
    /// \brief Constructor from a friction coefficient mu
    ///
    /// \param[in] mu Friction coefficient.
    ///
    explicit CoulombFrictionConeTpl(const Scalar mu)
    : mu(mu)
    {
      assert(mu >= 0 && "mu must be positive");
    }

    /// \brief Copy constructor.
    CoulombFrictionConeTpl(const CoulombFrictionConeTpl & other) = default;

    /// \brief Copy operator
    CoulombFrictionConeTpl & operator=(const CoulombFrictionConeTpl & other) = default;

    /// \brief Cast operator
    template<typename NewScalar>
    CoulombFrictionConeTpl<NewScalar> cast() const
    {
      typedef CoulombFrictionConeTpl<NewScalar> ReturnType;
      return ReturnType(NewScalar(this->mu));
    }

    /// \brief Comparison operator
    bool operator==(const CoulombFrictionConeTpl & other) const
    {
      return base() == other.base() && mu == other.mu;
    }

    /// \brief Difference  operator
    bool operator!=(const CoulombFrictionConeTpl & other) const
    {
      return !(*this == other);
    }

    /// \brief Check whether a vector x lies within the cone.
    ///
    /// \param[in] f vector to check (assimilated to a  force vector).
    ///
    template<typename Vector3Like>
    bool isInside(const Eigen::MatrixBase<Vector3Like> & f, const Scalar prec = Scalar(0)) const
    {
      assert(mu >= 0 && "mu must be positive");
      assert(prec >= 0 && "prec should be positive");
      //      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Vector3Like, 3);
      assert(f.size() == 3 && "The input vector is of wrong size.");
      const Vector3 f_normalized = f.normalized();
      return f_normalized.template head<2>().norm() <= mu * f_normalized[2] + prec;
    }

    using Base::project;

    /// \brief Project a vector x onto the cone.
    ///
    /// \param[in] x a 3d vector to project.
    ///
    template<typename Vector3LikeIn, typename Vector3LikeOut>
    void project(
      const Eigen::MatrixBase<Vector3LikeIn> & x,
      const Eigen::MatrixBase<Vector3LikeOut> & res_) const
    {
      assert(mu >= 0 && "mu must be positive");
      //      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Vector3Like, 3);
      assert(x.size() == 3 && "The input vector is of wrong size.");
      typedef Eigen::Matrix<Scalar, 2, 1> Vector2Plain;

      const Scalar & z = x[2];
      const Scalar mu_z = mu * z;

      auto & res = res_.const_cast_derived();

      const Vector2Plain t = x.template head<2>();
      const Scalar t_norm = t.norm();

      if (mu * t_norm <= -z)
      {
        res.setZero();
        return;
      }
      else if (t_norm <= mu_z)
      {
        res = x;
        return;
      }
      else
      {
        res.template head<2>() = (mu / t_norm) * t;
        res[2] = 1;
        res.normalize();
        const Scalar scale = x.dot(res);
        res *= scale;
        return;
      }
    }

    /// \brief Project a vector x onto the cone with a matric specified by the diagonal matrix R.
    ///
    /// \param[in] x a 3d vector to project.
    /// \param[in] R a 3d vector representing the diagonal of the weight matrix. The tangential
    /// components (the first two) of R should be equal, assuming an isotropic scaling.
    ///
    template<typename Vector3Like1, typename Vector3Like2>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(Vector3Like1) weightedProject(
      const Eigen::MatrixBase<Vector3Like1> & x, const Eigen::MatrixBase<Vector3Like2> & R) const
    {
      assert(mu >= 0 && "mu must be positive");
      //      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Vector3Like, 3);
      assert(x.size() == 3 && "The input vector is of wrong size.");
      assert(R(2) > 0 && "R(2) must be strictly positive");
      assert(R(0) == R(1) && "R(0) must be equal to R(1)");

      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Vector3Like1) Vector3Plain;

      const Scalar weighted_mu = mu * math::sqrt(R(0) / R(2));
      const CoulombFrictionConeTpl weighted_cone(weighted_mu);
      const Vector3Plain R_sqrt = R.cwiseSqrt();
      const Vector3Plain R_sqrt_times_x = R_sqrt.array() * x.array();
      Vector3Plain res = weighted_cone.project(R_sqrt_times_x).array() / R_sqrt.array();
      return res;
    }

    /// \brief Compute the complementary shift associted to the Coulomb friction cone for
    /// complementarity satisfaction in complementary problems.
    ///
    /// \param[in] v a dual vector.
    ///
    template<typename Vector3Like>
    typename Eigen::Matrix<Scalar, 3, 1>
    computeNormalCorrection(const Eigen::MatrixBase<Vector3Like> & v) const
    {
      //      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Vector3Like, 3);
      assert(v.size() == 3);
      typedef Eigen::Matrix<Scalar, 3, 1> Vector3Plain;

      Vector3Plain res;
      res.template head<2>().setZero();
      res[2] = mu * v.template head<2>().norm();

      return res;
    }

    /// \brief Compute the radial projection associted to the Coulomb friction cone.
    ///
    /// \param[in] f a force vector.
    ///
    template<typename Vector3Like>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(Vector3Like)
      computeRadialProjection(const Eigen::MatrixBase<Vector3Like> & f) const
    {
      //      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Vector3Like, 3);
      assert(f.size() == 3 && "The input vector is of wrong size.");
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Vector3Like) Vector3Plain;

      Vector3Plain res;
      const auto & ft = f.template head<2>();
      const Scalar ft_norm = ft.norm();

      res[2] = math::max(Scalar(0), f[2]);
      const Scalar mu_fz = mu * res[2];
      if (ft_norm > mu_fz)
      {
        res.template head<2>() = Scalar(mu_fz / ft_norm) * ft;
      }
      else
        res.template head<2>() = ft;

      return res;
    }

    template<typename Vector3Like1, typename Vector3Like2>
    Scalar computeContactComplementarity(
      const Eigen::MatrixBase<Vector3Like1> & v, const Eigen::MatrixBase<Vector3Like2> & f) const
    {
      typedef Eigen::Matrix<Scalar, 3, 1> Vector3Plain;
      return math::fabs(f.dot(Vector3Plain(v + computeNormalCorrection(v))));
    }

    template<typename Vector3Like1, typename Vector3Like2>
    Scalar computeConicComplementarity(
      const Eigen::MatrixBase<Vector3Like1> & v, const Eigen::MatrixBase<Vector3Like2> & f) const
    {
      return math::fabs(f.dot(v));
    }

    /// \brief Returns the dual cone associated to this.
    DualCone dual() const
    {
      return DualCone(mu);
    }

    /// \brief Returns the dimension of the cone.
    static int dim()
    {
      return 3;
    }

    int size() const
    {
      return dim();
    }

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    /// \var Friction coefficient
    Scalar mu;
  }; // CoulombFrictionConeTpl

  ///  \brief Dual of the 3d Coulomb friction cone.
  template<typename _Scalar>
  struct DualCoulombFrictionConeTpl : ConeBase<DualCoulombFrictionConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
    typedef typename traits<DualCoulombFrictionConeTpl>::DualCone DualCone;
    typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
    typedef ConeBase<DualCoulombFrictionConeTpl> Base;

    ///
    /// \brief Default constructor
    ///
    DualCoulombFrictionConeTpl()
    : mu(+std::numeric_limits<Scalar>::infinity())
    {
    }

    ///
    /// \brief Constructor from a friction coefficient mu
    ///
    /// \param[in] mu Friction coefficient.
    ///
    explicit DualCoulombFrictionConeTpl(const Scalar mu)
    : mu(mu)
    {
      assert(mu >= 0 && "mu must be positive");
    }

    /// \brief Copy constructor.
    DualCoulombFrictionConeTpl(const DualCoulombFrictionConeTpl & other) = default;

    /// \brief Copy operator
    DualCoulombFrictionConeTpl & operator=(const DualCoulombFrictionConeTpl & other) = default;

    /// \brief Comparison operator
    bool operator==(const DualCoulombFrictionConeTpl & other) const
    {
      return base() == other.base() && mu == other.mu;
    }

    /// \brief Difference  operator
    bool operator!=(const DualCoulombFrictionConeTpl & other) const
    {
      return !(*this == other);
    }

    /// \brief Check whether a vector v lies within the cone.
    ///
    /// \param[in] v vector to check (assimilated to a linear velocity).
    ///
    template<typename Vector3Like>
    bool isInside(const Eigen::MatrixBase<Vector3Like> & v, const Scalar prec = Scalar(0)) const
    {
      assert(mu >= 0 && "mu must be positive");
      //      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Vector3Like, 3);
      assert(v.size() == 3 && "The input vector is of wrong size.");
      const Vector3 v_normalized = v.normalized();
      return mu * v_normalized.template head<2>().norm() <= v_normalized[2] + prec;
    }

    using Base::project;
    /// \brief Project a vector x onto the cone
    template<typename Vector3LikeIn, typename Vector3LikeOut>
    void project(
      const Eigen::MatrixBase<Vector3LikeIn> & x,
      const Eigen::MatrixBase<Vector3LikeOut> & res_) const
    {
      assert(mu >= 0 && "mu must be positive");
      //      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Vector3Like, 3);
      assert(x.size() == 3 && "The input vector is of wrong size.");
      const Scalar & z = x[2];

      auto & res = res_.const_cast_derived();

      const Eigen::Matrix<Scalar, 2, 1> t = x.template head<2>();
      const Scalar t_norm = t.norm();

      if (t_norm <= -mu * z)
      {
        res.setZero();
        return;
      }
      else if (mu * t_norm <= z)
      {
        res = x;
        return;
      }
      else
      {
        res.template head<2>() = t;
        res[2] = mu * t_norm;
        res.normalize();
        const Scalar scale = x.dot(res);
        res *= scale;
        return;
      }
    }

    /// \brief Returns the dimension of the cone.
    static int dim()
    {
      return 3;
    }

    int size() const
    {
      return dim();
    }

    /// \brief Returns the dual cone associated to this.    ///
    DualCone dual() const
    {
      return DualCone(mu);
    }

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    /// \var Friction coefficient
    Scalar mu;

  }; // DualCoulombFrictionConeTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_coulomb_friction_cone_hpp__
