//
// Copyright (c) 2015-2020 CNRS INRIA
// Copyright (c) 2015-2016 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_multibody_joint_spherical_ZYX_hpp__
#define __pinocchio_multibody_joint_spherical_ZYX_hpp__

#include "pinocchio/macros.hpp"
#include "pinocchio/multibody/joint/joint-base.hpp"
#include "pinocchio/multibody/joint/joint-spherical.hpp"
#include "pinocchio/multibody/joint-motion-subspace.hpp"
#include "pinocchio/math/sincos.hpp"
#include "pinocchio/math/matrix.hpp"
#include "pinocchio/spatial/inertia.hpp"
#include "pinocchio/spatial/skew.hpp"

namespace pinocchio
{
  template<typename Scalar, int Options>
  struct JointMotionSubspaceSphericalZYXTpl;

  template<typename _Scalar, int _Options>
  struct traits<JointMotionSubspaceSphericalZYXTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    enum
    {
      LINEAR = 0,
      ANGULAR = 3
    };

    typedef MotionSphericalTpl<Scalar, Options> JointMotion;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> JointForce;
    typedef Eigen::Matrix<Scalar, 6, 3, Options> DenseBase;
    typedef Eigen::Matrix<Scalar, 3, 3, Options> ReducedSquaredMatrix;

    typedef DenseBase MatrixReturnType;
    typedef const DenseBase ConstMatrixReturnType;

    typedef ReducedSquaredMatrix StDiagonalMatrixSOperationReturnType;
  }; // struct traits struct ConstraintRotationalSubspace

  template<typename _Scalar, int _Options>
  struct JointMotionSubspaceSphericalZYXTpl
  : public JointMotionSubspaceBase<JointMotionSubspaceSphericalZYXTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    PINOCCHIO_CONSTRAINT_TYPEDEF_TPL(JointMotionSubspaceSphericalZYXTpl)

    enum
    {
      NV = 3
    };
    typedef Eigen::Matrix<Scalar, 3, 3, Options> Matrix3;

    JointMotionSubspaceSphericalZYXTpl()
    {
    }

    template<typename Matrix3Like>
    JointMotionSubspaceSphericalZYXTpl(const Eigen::MatrixBase<Matrix3Like> & subspace)
    : m_S(subspace)
    {
      EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix3Like, 3, 3);
    }

    template<typename Vector3Like>
    JointMotion __mult__(const Eigen::MatrixBase<Vector3Like> & v) const
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Vector3Like, 3);
      return JointMotion(m_S * v);
    }

    Matrix3 & operator()()
    {
      return m_S;
    }
    const Matrix3 & operator()() const
    {
      return m_S;
    }

    int nv_impl() const
    {
      return NV;
    }

    struct ConstraintTranspose
    : JointMotionSubspaceTransposeBase<JointMotionSubspaceSphericalZYXTpl>
    {
      const JointMotionSubspaceSphericalZYXTpl & ref;
      ConstraintTranspose(const JointMotionSubspaceSphericalZYXTpl & ref)
      : ref(ref)
      {
      }

      template<typename Derived>
      const typename MatrixMatrixProduct<
        Eigen::Transpose<const Matrix3>,
        Eigen::Block<const typename ForceDense<Derived>::Vector6, 3, 1>>::type
      operator*(const ForceDense<Derived> & phi) const
      {
        return ref.m_S.transpose() * phi.angular();
      }

      /* [CRBA]  MatrixBase operator* (Constraint::Transpose S, ForceSet::Block) */
      template<typename D>
      const typename MatrixMatrixProduct<
        typename Eigen::Transpose<const Matrix3>,
        typename Eigen::MatrixBase<const D>::template NRowsBlockXpr<3>::Type>::type
      operator*(const Eigen::MatrixBase<D> & F) const
      {
        EIGEN_STATIC_ASSERT(
          D::RowsAtCompileTime == 6, THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE)
        return ref.m_S.transpose() * F.template middleRows<3>(ANGULAR);
      }
    }; // struct ConstraintTranspose

    ConstraintTranspose transpose() const
    {
      return ConstraintTranspose(*this);
    }

    DenseBase matrix_impl() const
    {
      DenseBase S;
      S.template middleRows<3>(LINEAR).setZero();
      S.template middleRows<3>(ANGULAR) = m_S;
      return S;
    }

    //      const typename Eigen::ProductReturnType<
    //      const ConstraintDense,
    //      const Matrix3
    //      >::Type
    template<typename S1, int O1>
    Eigen::Matrix<Scalar, 6, 3, Options> se3Action(const SE3Tpl<S1, O1> & m) const
    {
      //        Eigen::Matrix <Scalar,6,3,Options> X_subspace;
      //        X_subspace.template block <3,3> (Motion::LINEAR, 0) = skew (m.translation ()) *
      //        m.rotation (); X_subspace.template block <3,3> (Motion::ANGULAR, 0) = m.rotation ();
      //
      //        return (X_subspace * m_S).eval();

      Eigen::Matrix<Scalar, 6, 3, Options> result;

      // ANGULAR
      result.template middleRows<3>(ANGULAR).noalias() = m.rotation() * m_S;

      // LINEAR
      cross(
        m.translation(), result.template middleRows<3>(Motion::ANGULAR),
        result.template middleRows<3>(LINEAR));

      return result;
    }

    template<typename S1, int O1>
    Eigen::Matrix<Scalar, 6, 3, Options> se3ActionInverse(const SE3Tpl<S1, O1> & m) const
    {
      Eigen::Matrix<Scalar, 6, 3, Options> result;

      // LINEAR
      cross(m.translation(), m_S, result.template middleRows<3>(ANGULAR));
      result.template middleRows<3>(LINEAR).noalias() =
        -m.rotation().transpose() * result.template middleRows<3>(ANGULAR);

      // ANGULAR
      result.template middleRows<3>(ANGULAR).noalias() = m.rotation().transpose() * m_S;

      return result;
    }

    template<typename MotionDerived>
    DenseBase motionAction(const MotionDense<MotionDerived> & m) const
    {
      const typename MotionDerived::ConstLinearType v = m.linear();
      const typename MotionDerived::ConstAngularType w = m.angular();

      DenseBase res;
      cross(v, m_S, res.template middleRows<3>(LINEAR));
      cross(w, m_S, res.template middleRows<3>(ANGULAR));

      return res;
    }

    Matrix3 & angularSubspace()
    {
      return m_S;
    }
    const Matrix3 & angularSubspace() const
    {
      return m_S;
    }

    bool isEqual(const JointMotionSubspaceSphericalZYXTpl & other) const
    {
      return internal::comparison_eq(m_S, other.m_S);
    }

  protected:
    Matrix3 m_S;

  }; // struct JointMotionSubspaceSphericalZYXTpl

  namespace details
  {
    template<typename Scalar, int Options>
    struct StDiagonalMatrixSOperation<JointMotionSubspaceSphericalZYXTpl<Scalar, Options>>
    {
      typedef JointMotionSubspaceSphericalZYXTpl<Scalar, Options> Constraint;
      typedef typename traits<Constraint>::StDiagonalMatrixSOperationReturnType ReturnType;

      static ReturnType run(const JointMotionSubspaceBase<Constraint> & constraint)
      {
        return constraint.matrix().transpose() * constraint.matrix();
      }
    };
  } // namespace details

  /* [CRBA] ForceSet operator* (Inertia Y,Constraint S) */
  template<typename S1, int O1, typename S2, int O2>
  Eigen::Matrix<S1, 6, 3, O1>
  operator*(const InertiaTpl<S1, O1> & Y, const JointMotionSubspaceSphericalZYXTpl<S2, O2> & S)
  {
    typedef typename InertiaTpl<S1, O1>::Symmetric3 Symmetric3;
    typedef JointMotionSubspaceSphericalZYXTpl<S2, O2> Constraint;
    Eigen::Matrix<S1, 6, 3, O1> M;
    alphaSkew(-Y.mass(), Y.lever(), M.template middleRows<3>(Constraint::LINEAR));
    M.template middleRows<3>(Constraint::ANGULAR) =
      (Y.inertia() - typename Symmetric3::AlphaSkewSquare(Y.mass(), Y.lever())).matrix();

    return (M * S.angularSubspace()).eval();
  }

  /* [ABA] Y*S operator (Inertia Y,Constraint S) */
  //  inline Eigen::Matrix<context::Scalar,6,3>
  template<typename Matrix6Like, typename S2, int O2>
  const typename MatrixMatrixProduct<
    typename Eigen::internal::remove_const<
      typename SizeDepType<3>::ColsReturn<Matrix6Like>::ConstType>::type,
    typename JointMotionSubspaceSphericalZYXTpl<S2, O2>::Matrix3>::type
  operator*(
    const Eigen::MatrixBase<Matrix6Like> & Y, const JointMotionSubspaceSphericalZYXTpl<S2, O2> & S)
  {
    EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix6Like, 6, 6);
    return Y.derived().template middleCols<3>(Inertia::ANGULAR) * S.angularSubspace();
  }

  template<typename S1, int O1>
  struct SE3GroupAction<JointMotionSubspaceSphericalZYXTpl<S1, O1>>
  {
    //      typedef const typename Eigen::ProductReturnType<
    //      Eigen::Matrix <context::Scalar,6,3,0>,
    //      Eigen::Matrix <context::Scalar,3,3,0>
    //      >::Type Type;
    typedef Eigen::Matrix<S1, 6, 3, O1> ReturnType;
  };

  template<typename S1, int O1, typename MotionDerived>
  struct MotionAlgebraAction<JointMotionSubspaceSphericalZYXTpl<S1, O1>, MotionDerived>
  {
    typedef Eigen::Matrix<S1, 6, 3, O1> ReturnType;
  };

  template<typename Scalar, int Options>
  struct JointSphericalZYXTpl;

  template<typename _Scalar, int _Options>
  struct traits<JointSphericalZYXTpl<_Scalar, _Options>>
  {
    enum
    {
      NQ = 3,
      NV = 3,
      NVExtended = 3
    };
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef JointDataSphericalZYXTpl<Scalar, Options> JointDataDerived;
    typedef JointModelSphericalZYXTpl<Scalar, Options> JointModelDerived;
    typedef JointMotionSubspaceSphericalZYXTpl<Scalar, Options> Constraint_t;
    typedef SE3Tpl<Scalar, Options> Transformation_t;
    typedef MotionSphericalTpl<Scalar, Options> Motion_t;
    typedef MotionSphericalTpl<Scalar, Options> Bias_t;

    // [ABA]
    typedef Eigen::Matrix<Scalar, 6, NV, Options> U_t;
    typedef Eigen::Matrix<Scalar, NV, NV, Options> D_t;
    typedef Eigen::Matrix<Scalar, 6, NV, Options> UD_t;

    typedef Eigen::Matrix<Scalar, NQ, 1, Options> ConfigVector_t;
    typedef Eigen::Matrix<Scalar, NV, 1, Options> TangentVector_t;

    typedef std::false_type is_mimicable_t;

    PINOCCHIO_JOINT_DATA_BASE_ACCESSOR_DEFAULT_RETURN_TYPE
  };

  template<typename _Scalar, int _Options>
  struct traits<JointDataSphericalZYXTpl<_Scalar, _Options>>
  {
    typedef JointSphericalZYXTpl<_Scalar, _Options> JointDerived;
    typedef _Scalar Scalar;
  };

  template<typename _Scalar, int _Options>
  struct traits<JointModelSphericalZYXTpl<_Scalar, _Options>>
  {
    typedef JointSphericalZYXTpl<_Scalar, _Options> JointDerived;
    typedef _Scalar Scalar;
  };

  template<typename _Scalar, int _Options>
  struct JointDataSphericalZYXTpl
  : public JointDataBase<JointDataSphericalZYXTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef JointSphericalZYXTpl<_Scalar, _Options> JointDerived;
    PINOCCHIO_JOINT_DATA_TYPEDEF_TEMPLATE(JointDerived);
    PINOCCHIO_JOINT_DATA_BASE_DEFAULT_ACCESSOR

    ConfigVector_t joint_q;
    TangentVector_t joint_v;

    Constraint_t S;
    Transformation_t M;
    Motion_t v;
    Bias_t c;

    // [ABA] specific data
    U_t U;
    D_t Dinv;
    UD_t UDinv;
    D_t StU;

    JointDataSphericalZYXTpl()
    : joint_q(ConfigVector_t::Zero())
    , joint_v(TangentVector_t::Zero())
    , S(Constraint_t::Matrix3::Zero())
    , M(Transformation_t::Identity())
    , v(Motion_t::Vector3::Zero())
    , c(Bias_t::Vector3::Zero())
    , U(U_t::Zero())
    , Dinv(D_t::Zero())
    , UDinv(UD_t::Zero())
    , StU(D_t::Zero())
    {
    }

    static std::string classname()
    {
      return std::string("JointDataSphericalZYX");
    }
    std::string shortname() const
    {
      return classname();
    }

  }; // strcut JointDataSphericalZYXTpl

  PINOCCHIO_JOINT_CAST_TYPE_SPECIALIZATION(JointModelSphericalZYXTpl);
  template<typename _Scalar, int _Options>
  struct JointModelSphericalZYXTpl
  : public JointModelBase<JointModelSphericalZYXTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef JointSphericalZYXTpl<_Scalar, _Options> JointDerived;
    PINOCCHIO_JOINT_TYPEDEF_TEMPLATE(JointDerived);

    typedef JointModelBase<JointModelSphericalZYXTpl> Base;
    using Base::id;
    using Base::idx_q;
    using Base::idx_v;
    using Base::idx_vExtended;
    using Base::setIndexes;

    JointDataDerived createData() const
    {
      return JointDataDerived();
    }

    const std::vector<bool> hasConfigurationLimit() const
    {
      return {true, true, true};
    }

    const std::vector<bool> hasConfigurationLimitInTangent() const
    {
      return {true, true, true};
    }

    template<typename ConfigVector>
    void calc(JointDataDerived & data, const typename Eigen::MatrixBase<ConfigVector> & qs) const
    {
      data.joint_q = qs.template segment<NQ>(idx_q());

      Scalar c0, s0;
      SINCOS(data.joint_q(0), &s0, &c0);
      Scalar c1, s1;
      SINCOS(data.joint_q(1), &s1, &c1);
      Scalar c2, s2;
      SINCOS(data.joint_q(2), &s2, &c2);

      data.M.rotation() << c0 * c1, c0 * s1 * s2 - s0 * c2, c0 * s1 * c2 + s0 * s2, s0 * c1,
        s0 * s1 * s2 + c0 * c2, s0 * s1 * c2 - c0 * s2, -s1, c1 * s2, c1 * c2;

      data.S.angularSubspace() << -s1, Scalar(0), Scalar(1), c1 * s2, c2, Scalar(0), c1 * c2, -s2,
        Scalar(0);
    }

    template<typename TangentVector>
    void
    calc(JointDataDerived & data, const Blank, const typename Eigen::MatrixBase<TangentVector> & vs)
      const
    {
      data.joint_v = vs.template segment<NV>(idx_v());
      data.v().noalias() = data.S.angularSubspace() * data.joint_v;

      // TODO(jcarpent): to be done
      // #define q_dot data.joint_v
      //       data.c()(0) = -c1 * q_dot(0) * q_dot(1);
      //       data.c()(1) = -s1 * s2 * q_dot(0) * q_dot(1) + c1 * c2 * q_dot(0) * q_dot(2) - s2 *
      //       q_dot(1) * q_dot(2); data.c()(2) = -s1 * c2 * q_dot(0) * q_dot(1) - c1 * s2 *
      //       q_dot(0) * q_dot(2) - c2 * q_dot(1) * q_dot(2);
      // #undef q_dot
    }

    template<typename ConfigVector, typename TangentVector>
    void calc(
      JointDataDerived & data,
      const typename Eigen::MatrixBase<ConfigVector> & qs,
      const typename Eigen::MatrixBase<TangentVector> & vs) const
    {
      data.joint_q = qs.template segment<NQ>(idx_q());

      Scalar c0, s0;
      SINCOS(data.joint_q(0), &s0, &c0);
      Scalar c1, s1;
      SINCOS(data.joint_q(1), &s1, &c1);
      Scalar c2, s2;
      SINCOS(data.joint_q(2), &s2, &c2);

      data.M.rotation() << c0 * c1, c0 * s1 * s2 - s0 * c2, c0 * s1 * c2 + s0 * s2, s0 * c1,
        s0 * s1 * s2 + c0 * c2, s0 * s1 * c2 - c0 * s2, -s1, c1 * s2, c1 * c2;

      data.S.angularSubspace() << -s1, Scalar(0), Scalar(1), c1 * s2, c2, Scalar(0), c1 * c2, -s2,
        Scalar(0);

      data.joint_v = vs.template segment<NV>(idx_v());
      data.v().noalias() = data.S.angularSubspace() * data.joint_v;

#define q_dot data.joint_v
      data.c()(0) = -c1 * q_dot(0) * q_dot(1);
      data.c()(1) =
        -s1 * s2 * q_dot(0) * q_dot(1) + c1 * c2 * q_dot(0) * q_dot(2) - s2 * q_dot(1) * q_dot(2);
      data.c()(2) =
        -s1 * c2 * q_dot(0) * q_dot(1) - c1 * s2 * q_dot(0) * q_dot(2) - c2 * q_dot(1) * q_dot(2);
#undef q_dot
    }

    template<typename VectorLike, typename Matrix6Like>
    void calc_aba(
      JointDataDerived & data,
      const Eigen::MatrixBase<VectorLike> & armature,
      const Eigen::MatrixBase<Matrix6Like> & I,
      const bool update_I) const
    {
      data.U.noalias() = I.template middleCols<3>(Motion::ANGULAR) * data.S.angularSubspace();
      data.StU.noalias() =
        data.S.angularSubspace().transpose() * data.U.template middleRows<3>(Motion::ANGULAR);
      data.StU.diagonal() += armature;

      matrix_inversion(data.StU, data.Dinv);

      data.UDinv.noalias() = data.U * data.Dinv;

      if (update_I)
        I.const_cast_derived().noalias() -= data.UDinv * data.U.transpose();
    }

    static std::string classname()
    {
      return std::string("JointModelSphericalZYX");
    }
    std::string shortname() const
    {
      return classname();
    }

    /// \returns An expression of *this with the Scalar type casted to NewScalar.
    template<typename NewScalar>
    JointModelSphericalZYXTpl<NewScalar, Options> cast() const
    {
      typedef JointModelSphericalZYXTpl<NewScalar, Options> ReturnType;
      ReturnType res;
      res.setIndexes(id(), idx_q(), idx_v(), idx_vExtended());
      return res;
    }

  }; // struct JointModelSphericalZYXTpl

} // namespace pinocchio

#include <boost/type_traits.hpp>

namespace boost
{
  template<typename Scalar, int Options>
  struct has_nothrow_constructor<::pinocchio::JointModelSphericalZYXTpl<Scalar, Options>>
  : public integral_constant<bool, true>
  {
  };

  template<typename Scalar, int Options>
  struct has_nothrow_copy<::pinocchio::JointModelSphericalZYXTpl<Scalar, Options>>
  : public integral_constant<bool, true>
  {
  };

  template<typename Scalar, int Options>
  struct has_nothrow_constructor<::pinocchio::JointDataSphericalZYXTpl<Scalar, Options>>
  : public integral_constant<bool, true>
  {
  };

  template<typename Scalar, int Options>
  struct has_nothrow_copy<::pinocchio::JointDataSphericalZYXTpl<Scalar, Options>>
  : public integral_constant<bool, true>
  {
  };
} // namespace boost

#endif // ifndef __pinocchio_multibody_joint_spherical_ZYX_hpp__
