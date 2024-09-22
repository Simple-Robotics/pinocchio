//
// Copyright (c) 2015-2018 CNRS
// Copyright (c) 2018-2024 INRIA
// Copyright (c) 2015-2016 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_spatial_force_tpl_hpp__
#define __pinocchio_spatial_force_tpl_hpp__

namespace pinocchio
{
  template<typename _Scalar, int _Options>
  struct traits<ForceTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1, _Options> Vector3;
    typedef Eigen::Matrix<Scalar, 6, 1, _Options> Vector6;
    typedef Eigen::Matrix<Scalar, 6, 6, _Options> Matrix6;
    typedef typename PINOCCHIO_EIGEN_REF_CONST_TYPE(Vector6) ToVectorConstReturnType;
    typedef typename PINOCCHIO_EIGEN_REF_TYPE(Vector6) ToVectorReturnType;
    typedef typename Vector6::template FixedSegmentReturnType<3>::Type LinearType;
    typedef typename Vector6::template FixedSegmentReturnType<3>::Type AngularType;
    typedef typename Vector6::template ConstFixedSegmentReturnType<3>::Type ConstLinearType;
    typedef typename Vector6::template ConstFixedSegmentReturnType<3>::Type ConstAngularType;
    typedef ForceTpl<Scalar, _Options> ForcePlain;
    enum
    {
      LINEAR = 0,
      ANGULAR = 3,
      Options = _Options
    };

    typedef ForceRef<Vector6> ForceRefType;
  }; // traits ForceTpl

  ///
  /// \brief This class represents a SpatialForce composed of a linear component representing a 3d
  /// force vector and an angular component corresponding to a 3d torque vector. Internally, the
  /// data are stored in a compact 6d vector.
  ///
  /// \tparam _Scalar Scalar type of the ForceTpl
  /// \tparam _Options Eigen Options related to the internal 6d vector storing the data. In most
  /// cases, you don't have to worry about this templated quantities.
  ///
  template<typename _Scalar, int _Options>
  class ForceTpl : public ForceDense<ForceTpl<_Scalar, _Options>>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef ForceDense<ForceTpl> Base;
    FORCE_TYPEDEF_TPL(ForceTpl);
    enum
    {
      Options = _Options
    };

    using Base::operator=;
    using Base::operator!=;
    using Base::angular;
    using Base::linear;

    /// \brief Default constructor
    ForceTpl()
    {
    }

    ///
    /// \brief Constructor from a force vector and a torque vector, both considered as 3d vectors.
    ///
    /// \tparam ForceVectorLike type of the force vector.
    /// \tparam TorqueVectorLike type of the torque vector.
    ///
    /// \param[in] force force vector associated with the linear component of the spatial force.
    /// \param[in] torque torque vector associated with the angular component of the spatial force.
    ///
    template<typename ForceVectorLike, typename TorqueVectorLike>
    ForceTpl(
      const Eigen::MatrixBase<ForceVectorLike> & force,
      const Eigen::MatrixBase<TorqueVectorLike> & torque)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(ForceVectorLike, 3);
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(TorqueVectorLike, 3);
      linear() = force;
      angular() = torque;
    }

    ///
    /// \brief Constructor from a 6d vector representing a spatial force.
    ///
    /// \tparam Vector6Like type of the vector
    ///
    /// \param[in] f 6d vector whose values represent the data associated with a spatial vector.
    ///
    template<typename Vector6Like>
    explicit ForceTpl(const Eigen::MatrixBase<Vector6Like> & f)
    : m_data(f)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Vector6Like, 6);
    }

    ///
    /// \brief Copy constructor
    ///
    /// \param[in] other to copy into *this.
    ///
    ForceTpl(const ForceTpl & other)
    : m_data(other.toVector())
    {
    }

    ///
    /// \brief Copy constructor with a cast
    ///
    ///
    /// \param[in] other to copy into *this after a cast.
    ///
    template<typename S2, int O2>
    explicit ForceTpl(const ForceTpl<S2, O2> & other)
    {
      *this = other.template cast<Scalar>();
    }

    ///
    /// \brief Copy constructor from a the base class.
    ///
    /// \param[in] other to copy into *this.
    ///
    template<typename ForceDerived>
    explicit ForceTpl(const ForceBase<ForceDerived> & other)
    {
      *this = other.derived();
    }

    ///
    /// \brief Copy constructor from a spatial dense force.
    ///
    /// \param[in] other to copy into *this.
    ///
    template<typename M2>
    explicit ForceTpl(const ForceDense<M2> & other)
    {
      linear() = other.linear();
      angular() = other.angular();
    }

    ///
    /// \brief Copy operator
    ///
    /// \param[in] other to copy into *this.
    ///
    ForceTpl & operator=(const ForceTpl & other) // Copy assignment operator
    {
      m_data = other.toVector();
      return *this;
    }

    ///
    /// \brief Returns a spatial force set to zero.
    ///
    static ForceTpl Zero()
    {
      return ForceTpl(Vector6::Zero());
    }

    ///
    /// \brief Returns a spatial force set to random.
    ///
    /// \remark This static constructor calls the random function of Eigen.
    ///
    static ForceTpl Random()
    {
      return ForceTpl(Vector6::Random());
    }

    ToVectorConstReturnType toVector_impl() const
    {
      return m_data;
    }
    ToVectorReturnType toVector_impl()
    {
      return m_data;
    }

    // Getters
    ConstAngularType angular_impl() const
    {
      return m_data.template segment<3>(ANGULAR);
    }
    ConstLinearType linear_impl() const
    {
      return m_data.template segment<3>(LINEAR);
    }
    AngularType angular_impl()
    {
      return m_data.template segment<3>(ANGULAR);
    }
    LinearType linear_impl()
    {
      return m_data.template segment<3>(LINEAR);
    }

    template<typename V3>
    void angular_impl(const Eigen::MatrixBase<V3> & w)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(V3, 3);
      angular_impl() = w;
    }
    template<typename V3>
    void linear_impl(const Eigen::MatrixBase<V3> & v)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(V3, 3);
      linear_impl() = v;
    }

    /// \returns *this as a ForceRef vector
    ForceRef<Vector6> ref()
    {
      return ForceRef<Vector6>(m_data);
    }

    /// \returns a cast version of *this.
    template<typename NewScalar>
    ForceTpl<NewScalar, Options> cast() const
    {
      typedef ForceTpl<NewScalar, Options> ReturnType;
      ReturnType res(linear().template cast<NewScalar>(), angular().template cast<NewScalar>());
      return res;
    }

  protected:
    Vector6 m_data;

  }; // class ForceTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_spatial_force_tpl_hpp__
