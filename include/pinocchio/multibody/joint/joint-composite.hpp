//
// Copyright (c) 2016-2021 CNRS INRIA
//

#ifndef __pinocchio_multibody_joint_composite_hpp__
#define __pinocchio_multibody_joint_composite_hpp__

#include "pinocchio/multibody/joint/fwd.hpp"
#include "pinocchio/multibody/joint/joint-collection.hpp"
#include "pinocchio/multibody/joint/joint-basic-visitors.hpp"
#include "pinocchio/container/aligned-vector.hpp"
#include "pinocchio/spatial/act-on-set.hpp"

#include "pinocchio/serialization/fwd.hpp"

namespace pinocchio
{

  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  struct JointCompositeTpl;

  template<typename _Scalar, int _Options, template<typename S, int O> class JointCollectionTpl>
  struct traits<JointCompositeTpl<_Scalar, _Options, JointCollectionTpl>>
  {
    typedef _Scalar Scalar;

    enum
    {
      Options = _Options,
      NQ = Eigen::Dynamic,
      NV = Eigen::Dynamic,
      NVExtended = Eigen::Dynamic
    };

    typedef JointCollectionTpl<Scalar, Options> JointCollection;
    typedef JointDataCompositeTpl<Scalar, Options, JointCollectionTpl> JointDataDerived;
    typedef JointModelCompositeTpl<Scalar, Options, JointCollectionTpl> JointModelDerived;
    typedef JointMotionSubspaceTpl<Eigen::Dynamic, Scalar, Options> Constraint_t;
    typedef SE3Tpl<Scalar, Options> Transformation_t;
    typedef MotionTpl<Scalar, Options> Motion_t;
    typedef MotionTpl<Scalar, Options> Bias_t;

    // [ABA]
    typedef Eigen::Matrix<Scalar, 6, Eigen::Dynamic, Options> U_t;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> D_t;
    typedef Eigen::Matrix<Scalar, 6, Eigen::Dynamic, Options> UD_t;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> ConfigVector_t;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> TangentVector_t;

    typedef std::false_type is_mimicable_t;

    PINOCCHIO_JOINT_DATA_BASE_ACCESSOR_DEFAULT_RETURN_TYPE
  };

  template<typename _Scalar, int _Options, template<typename S, int O> class JointCollectionTpl>
  struct traits<JointModelCompositeTpl<_Scalar, _Options, JointCollectionTpl>>
  {
    typedef JointCompositeTpl<_Scalar, _Options, JointCollectionTpl> JointDerived;
    typedef _Scalar Scalar;
  };

  template<typename _Scalar, int _Options, template<typename S, int O> class JointCollectionTpl>
  struct traits<JointDataCompositeTpl<_Scalar, _Options, JointCollectionTpl>>
  {
    typedef JointCompositeTpl<_Scalar, _Options, JointCollectionTpl> JointDerived;
    typedef _Scalar Scalar;
  };

  template<typename _Scalar, int _Options, template<typename S, int O> class JointCollectionTpl>
  struct JointDataCompositeTpl
  : public JointDataBase<JointDataCompositeTpl<_Scalar, _Options, JointCollectionTpl>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef JointDataBase<JointDataCompositeTpl> Base;
    typedef JointCompositeTpl<_Scalar, _Options, JointCollectionTpl> JointDerived;
    PINOCCHIO_JOINT_DATA_TYPEDEF_TEMPLATE(JointDerived);
    PINOCCHIO_JOINT_DATA_BASE_DEFAULT_ACCESSOR

    typedef JointCollectionTpl<Scalar, Options> JointCollection;
    typedef JointDataTpl<Scalar, Options, JointCollectionTpl> JointDataVariant;

    typedef PINOCCHIO_ALIGNED_STD_VECTOR(JointDataVariant) JointDataVector;

    // JointDataComposite()  {} // can become necessary if we want a vector of JointDataComposite ?

    JointDataCompositeTpl()
    : joints()
    , iMlast(0)
    , pjMi(0)
    , joint_q(ConfigVector_t::Zero(0))
    , joint_v(TangentVector_t::Zero(0))
    , S(0)
    , M(Transformation_t::Identity())
    , v(Motion_t::Zero())
    , c(Motion_t::Zero())
    , U(6, 0)
    , Dinv(0, 0)
    , UDinv(6, 0)
    , StU(0, 0)
    {
    }

    JointDataCompositeTpl(const JointDataVector & joint_data, const int nq, const int nv)
    : joints(joint_data)
    , iMlast(joint_data.size())
    , pjMi(joint_data.size())
    , joint_q(ConfigVector_t::Zero(nq))
    , joint_v(TangentVector_t::Zero(nv))
    , S(Constraint_t::Zero(nv))
    , M(Transformation_t::Identity())
    , v(Motion_t::Zero())
    , c(Motion_t::Zero())
    , U(U_t::Zero(6, nv))
    , Dinv(D_t::Zero(nv, nv))
    , UDinv(UD_t::Zero(6, nv))
    , StU(D_t::Zero(nv, nv))
    {
    }

    /// \brief Vector of joints
    JointDataVector joints;

    /// \brief Transforms from previous joint to last joint
    PINOCCHIO_ALIGNED_STD_VECTOR(Transformation_t) iMlast;

    /// \brief Transforms from previous joint to joint i
    PINOCCHIO_ALIGNED_STD_VECTOR(Transformation_t) pjMi;

    ConfigVector_t joint_q;
    TangentVector_t joint_v;

    Constraint_t S;
    Transformation_t M;
    Motion_t v;
    Bias_t c;

    // // [ABA] specific data
    U_t U;
    D_t Dinv;
    UD_t UDinv;
    D_t StU;

    static std::string classname()
    {
      return std::string("JointDataComposite");
    }
    std::string shortname() const
    {
      return classname();
    }

    void disp(std::ostream & os) const
    {
      os << "JointDataComposite containing following models:\n";
      for (typename JointDataVector::const_iterator it = joints.begin(); it != joints.end(); ++it)
        os << "  " << it->shortname() << std::endl;
    }
  };

  template<
    typename NewScalar,
    typename Scalar,
    int Options,
    template<typename S, int O> class JointCollectionTpl>
  struct CastType<NewScalar, JointModelCompositeTpl<Scalar, Options, JointCollectionTpl>>
  {
    typedef JointModelCompositeTpl<NewScalar, Options, JointCollectionTpl> type;
  };

  template<typename _Scalar, int _Options, template<typename S, int O> class JointCollectionTpl>
  struct JointModelCompositeTpl
  : public JointModelBase<JointModelCompositeTpl<_Scalar, _Options, JointCollectionTpl>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef JointModelBase<JointModelCompositeTpl> Base;
    typedef JointCompositeTpl<_Scalar, _Options, JointCollectionTpl> JointDerived;
    PINOCCHIO_JOINT_TYPEDEF_TEMPLATE(JointDerived);

    typedef JointCollectionTpl<Scalar, Options> JointCollection;
    typedef JointModelTpl<Scalar, Options, JointCollectionTpl> JointModelVariant;

    typedef SE3Tpl<Scalar, Options> SE3;
    typedef MotionTpl<Scalar, Options> Motion;
    typedef InertiaTpl<Scalar, Options> Inertia;

    typedef PINOCCHIO_ALIGNED_STD_VECTOR(JointModelVariant) JointModelVector;

    using Base::id;
    using Base::idx_q;
    using Base::idx_v;
    using Base::idx_vExtended;
    using Base::nq;
    using Base::nv;
    using Base::nvExtended;
    using Base::setIndexes;

    /// \brief Default contructor
    JointModelCompositeTpl()
    : joints()
    , jointPlacements()
    , m_nq(0)
    , m_nv(0)
    , m_nvExtended(0)
    , njoints(0)
    {
    }

    /// \brief Default contructor with a defined size
    JointModelCompositeTpl(const size_t size)
    : joints()
    , jointPlacements()
    , m_nq(0)
    , m_nv(0)
    , m_nvExtended(0)
    , njoints(0)
    {
      joints.reserve(size);
      jointPlacements.reserve(size);
      m_idx_q.reserve(size);
      m_idx_v.reserve(size);
      m_idx_vExtended.reserve(size);
      m_nqs.reserve(size);
      m_nvs.reserve(size);
      m_nvExtendeds.reserve(size);
    }

    ///
    /// \brief Constructor with one joint.
    ///
    /// \param jmodel Model of the first joint.
    /// \param placement Placement of the first joint w.r.t. the joint origin.
    ///
    template<typename JointModel>
    JointModelCompositeTpl(
      const JointModelBase<JointModel> & jmodel, const SE3 & placement = SE3::Identity())
    : joints(1, (JointModelVariant)jmodel.derived())
    , jointPlacements(1, placement)
    , m_nq(jmodel.nq())
    , m_nv(jmodel.nv())
    , m_nvExtended(jmodel.nvExtended())
    , m_idx_q(1, 0)
    , m_nqs(1, jmodel.nq())
    , m_idx_v(1, 0)
    , m_nvs(1, jmodel.nv())
    , m_idx_vExtended(1, 0)
    , m_nvExtendeds(1, jmodel.nvExtended())
    , njoints(1)
    {
    }

    ///
    /// \brief Copy constructor.
    ///
    /// \param other JointModel to copy.
    ///
    JointModelCompositeTpl(const JointModelCompositeTpl & other)
    : Base(other)
    , joints(other.joints)
    , jointPlacements(other.jointPlacements)
    , m_nq(other.m_nq)
    , m_nv(other.m_nv)
    , m_nvExtended(other.m_nvExtended)
    , m_idx_q(other.m_idx_q)
    , m_nqs(other.m_nqs)
    , m_idx_v(other.m_idx_v)
    , m_nvs(other.m_nvs)
    , m_idx_vExtended(other.m_idx_vExtended)
    , m_nvExtendeds(other.m_nvExtendeds)
    , njoints(other.njoints)
    {
    }

    ///
    /// \brief Add a joint to the vector of joints.
    ///
    /// \param jmodel Model of the joint to add.
    /// \param placement Placement of the joint relatively to its predecessor.
    ///
    /// \return A reference to *this
    ///
    template<typename JointModel>
    JointModelDerived &
    addJoint(const JointModelBase<JointModel> & jmodel, const SE3 & placement = SE3::Identity())
    {
      joints.push_back((JointModelVariant)jmodel.derived());
      jointPlacements.push_back(placement);

      m_nq += jmodel.nq();
      m_nv += jmodel.nv();
      m_nvExtended += jmodel.nvExtended();

      updateJointIndexes();
      njoints++;

      return *this;
    }

    JointDataDerived createData() const
    {
      typename JointDataDerived::JointDataVector jdata(joints.size());
      for (int i = 0; i < (int)joints.size(); ++i)
        jdata[(size_t)i] =
          ::pinocchio::createData<Scalar, Options, JointCollectionTpl>(joints[(size_t)i]);
      return JointDataDerived(jdata, nq(), nv());
    }

    const std::vector<bool> hasConfigurationLimit() const
    {
      std::vector<bool> vec;
      for (size_t i = 0; i < joints.size(); ++i)
      {
        const std::vector<bool> & joint_cf_limit = joints[i].hasConfigurationLimit();
        vec.insert(vec.end(), joint_cf_limit.begin(), joint_cf_limit.end());
      }
      return vec;
    }

    const std::vector<bool> hasConfigurationLimitInTangent() const
    {
      std::vector<bool> vec;
      for (size_t i = 0; i < joints.size(); ++i)
      {
        const std::vector<bool> & joint_cf_limit = joints[i].hasConfigurationLimitInTangent();
        vec.insert(vec.end(), joint_cf_limit.begin(), joint_cf_limit.end());
      }
      return vec;
    }

    template<typename, int, template<typename S, int O> class, typename>
    friend struct JointCompositeCalcZeroOrderStep;

    template<typename ConfigVectorType>
    void calc(JointDataDerived & data, const Eigen::MatrixBase<ConfigVectorType> & qs) const;

    template<typename, int, template<typename S, int O> class, typename, typename>
    friend struct JointCompositeCalcFirstOrderStep;

    template<typename ConfigVectorType, typename TangentVectorType>
    void calc(
      JointDataDerived & data,
      const Eigen::MatrixBase<ConfigVectorType> & qs,
      const Eigen::MatrixBase<TangentVectorType> & vs) const;

    template<typename TangentVectorType>
    void calc(
      JointDataDerived & data,
      const Blank blank,
      const Eigen::MatrixBase<TangentVectorType> & vs) const;

    template<typename VectorLike, typename Matrix6Like>
    void calc_aba(
      JointDataDerived & data,
      const Eigen::MatrixBase<VectorLike> & armature,
      const Eigen::MatrixBase<Matrix6Like> & I,
      const bool update_I) const
    {
      data.U.noalias() = I * data.S.matrix();
      data.StU.noalias() = data.S.matrix().transpose() * data.U;
      data.StU.diagonal() += armature;

      matrix_inversion(data.StU, data.Dinv);
      data.UDinv.noalias() = data.U * data.Dinv;

      if (update_I)
        I.const_cast_derived().noalias() -= data.UDinv * data.U.transpose();
    }

    // Declaration of overload : must be define after Lie group and joint visitors
    template<typename LieGroupMap>
    typename LieGroupMap::template operation<JointModelCompositeTpl>::type lieGroup_impl() const;

    int nv_impl() const
    {
      return m_nv;
    }
    int nq_impl() const
    {
      return m_nq;
    }
    int nvExtended_impl() const
    {
      return m_nvExtended;
    }

    /**
     * @brief      Update the indexes of subjoints in the stack
     */
    void setIndexes_impl(JointIndex id, int q, int v, int vExtended)
    {
      Base::setIndexes_impl(id, q, v, vExtended);
      updateJointIndexes();
    }

    static std::string classname()
    {
      return std::string("JointModelComposite");
    }
    std::string shortname() const
    {
      return classname();
    }

    JointModelCompositeTpl & operator=(const JointModelCompositeTpl & other)
    {
      Base::operator=(other);
      m_nq = other.m_nq;
      m_nv = other.m_nv;
      m_nvExtended = other.m_nvExtended;
      m_idx_q = other.m_idx_q;
      m_idx_v = other.m_idx_v;
      m_idx_vExtended = other.m_idx_vExtended;
      m_nqs = other.m_nqs;
      m_nvs = other.m_nvs;
      m_nvExtendeds = other.m_nvExtendeds;
      joints = other.joints;
      jointPlacements = other.jointPlacements;
      njoints = other.njoints;

      return *this;
    }

    using Base::isEqual;
    bool isEqual(const JointModelCompositeTpl & other) const
    {
      return Base::isEqual(other) && internal::comparison_eq(nq(), other.nq())
             && internal::comparison_eq(nv(), other.nv())
             && internal::comparison_eq(nvExtended(), other.nvExtended())
             && internal::comparison_eq(m_idx_q, other.m_idx_q)
             && internal::comparison_eq(m_idx_v, other.m_idx_v)
             && internal::comparison_eq(m_idx_vExtended, other.m_idx_vExtended)
             && internal::comparison_eq(m_nqs, other.m_nqs)
             && internal::comparison_eq(m_nvs, other.m_nvs)
             && internal::comparison_eq(m_nvExtendeds, other.m_nvExtendeds)
             && internal::comparison_eq(joints, other.joints)
             && internal::comparison_eq(jointPlacements, other.jointPlacements)
             && internal::comparison_eq(njoints, other.njoints);
    }

    /// \returns An expression of *this with the Scalar type casted to NewScalar.
    template<typename NewScalar>
    JointModelCompositeTpl<NewScalar, Options, JointCollectionTpl> cast() const
    {
      typedef JointModelCompositeTpl<NewScalar, Options, JointCollectionTpl> ReturnType;
      ReturnType res((size_t)njoints);
      res.setIndexes(id(), idx_q(), idx_v(), idx_vExtended());
      res.m_nq = m_nq;
      res.m_nv = m_nv;
      res.m_nvExtended = m_nvExtended;
      res.m_idx_q = m_idx_q;
      res.m_idx_v = m_idx_v;
      res.m_idx_vExtended = m_idx_vExtended;
      res.m_nqs = m_nqs;
      res.m_nvs = m_nvs;
      res.m_nvExtendeds = m_nvExtendeds;
      res.njoints = njoints;

      res.joints.resize(joints.size());
      res.jointPlacements.resize(jointPlacements.size());
      for (size_t k = 0; k < jointPlacements.size(); ++k)
      {
        res.joints[k] = joints[k].template cast<NewScalar>();
        res.jointPlacements[k] = jointPlacements[k].template cast<NewScalar>();
      }

      return res;
    }

    /// \brief Vector of joints contained in the joint composite.
    JointModelVector joints;
    /// \brief Vector of joint placements. Those placements correspond to the origin of the joint
    /// relatively to their parent.
    PINOCCHIO_ALIGNED_STD_VECTOR(SE3) jointPlacements;

    template<typename D>
    typename SizeDepType<NQ>::template SegmentReturn<D>::ConstType
    JointMappedConfigSelector(const Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_q, nq());
    }
    template<typename D>
    typename SizeDepType<NQ>::template SegmentReturn<D>::Type
    JointMappedConfigSelector(Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_q, nq());
    }

    template<typename D>
    typename SizeDepType<NQ>::template SegmentReturn<D>::ConstType
    jointConfigSelector(const Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_q, nq());
    }
    template<typename D>
    typename SizeDepType<NQ>::template SegmentReturn<D>::Type
    jointConfigSelector(Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_q, nq());
    }

    template<typename D>
    typename SizeDepType<NV>::template SegmentReturn<D>::ConstType
    JointMappedVelocitySelector(const Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_v, nv());
    }
    template<typename D>
    typename SizeDepType<NV>::template SegmentReturn<D>::Type
    JointMappedVelocitySelector(Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_v, nv());
    }

    template<typename D>
    typename SizeDepType<NV>::template SegmentReturn<D>::ConstType
    jointVelocitySelector(const Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_v, nv());
    }
    template<typename D>
    typename SizeDepType<NV>::template SegmentReturn<D>::Type
    jointVelocitySelector(Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_v, nv());
    }

    template<typename D>
    typename SizeDepType<NV>::template ColsReturn<D>::ConstType
    jointCols(const Eigen::MatrixBase<D> & A) const
    {
      return A.middleCols(Base::i_v, nv());
    }
    template<typename D>
    typename SizeDepType<NV>::template ColsReturn<D>::ConstType
    jointExtendedModelCols(const Eigen::MatrixBase<D> & A) const
    {
      return A.middleCols(Base::i_vExtended, nvExtended());
    }
    template<typename D>
    typename SizeDepType<NV>::template ColsReturn<D>::Type jointCols(Eigen::MatrixBase<D> & A) const
    {
      return A.middleCols(Base::i_v, nv());
    }
    template<typename D>
    typename SizeDepType<NV>::template ColsReturn<D>::Type
    jointExtendedModelCols(Eigen::MatrixBase<D> & A) const
    {
      return A.middleCols(Base::i_vExtended, nvExtended());
    }

    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template SegmentReturn<D>::ConstType
    jointConfigSelector_impl(const Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_q, nq());
    }
    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template SegmentReturn<D>::Type
    jointConfigSelector_impl(Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_q, nq());
    }

    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template SegmentReturn<D>::ConstType
    JointMappedConfigSelector_impl(const Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_q, nq());
    }
    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template SegmentReturn<D>::Type
    JointMappedConfigSelector_impl(Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_q, nq());
    }

    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template SegmentReturn<D>::ConstType
    JointMappedVelocitySelector_impl(const Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_v, nv());
    }
    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template SegmentReturn<D>::Type
    JointMappedVelocitySelector_impl(Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_v, nv());
    }

    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template SegmentReturn<D>::ConstType
    jointVelocitySelector_impl(const Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_v, nv());
    }
    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template SegmentReturn<D>::Type
    jointVelocitySelector_impl(Eigen::MatrixBase<D> & a) const
    {
      return a.segment(Base::i_v, nv());
    }

    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template ColsReturn<D>::ConstType
    jointCols_impl(const Eigen::MatrixBase<D> & A) const
    {
      return A.middleCols(Base::i_v, nv());
    }
    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template ColsReturn<D>::ConstType
    jointExtendedModelCols_impl(const Eigen::MatrixBase<D> & A) const
    {
      return A.middleCols(Base::i_vExtended, nvExtended());
    }
    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template ColsReturn<D>::Type
    jointCols_impl(Eigen::MatrixBase<D> & A) const
    {
      return A.middleCols(Base::i_v, nv());
    }
    template<typename D>
    typename SizeDepType<Eigen::Dynamic>::template ColsReturn<D>::Type
    jointExtendedModelCols_impl(Eigen::MatrixBase<D> & A) const
    {
      return A.middleCols(Base::i_vExtended, nvExtended());
    }

    void disp(std::ostream & os) const
    {
      typedef typename JointModelCompositeTpl<Scalar, Options, JointCollectionTpl>::JointModelVector
        JointModelVector;

      os << "JointModelComposite containing following models:\n";
      for (typename JointModelVector::const_iterator it = joints.begin(); it != joints.end(); ++it)
        os << "  " << it->shortname() << std::endl;
    }

  protected:
    friend struct Serialize<JointModelCompositeTpl>;

    template<typename, int, template<typename, int> class>
    friend struct JointModelCompositeTpl;

    /// \brief Update the indexes of the joints contained in the composition according
    /// to the position of the joint composite.
    void updateJointIndexes()
    {
      int idx_q = this->idx_q();
      int idx_v = this->idx_v();
      int idx_vExtended = this->idx_vExtended();

      m_idx_q.resize(joints.size());
      m_idx_v.resize(joints.size());
      m_idx_vExtended.resize(joints.size());
      m_nqs.resize(joints.size());
      m_nvs.resize(joints.size());
      m_nvExtendeds.resize(joints.size());

      for (size_t i = 0; i < joints.size(); ++i)
      {
        JointModelVariant & joint = joints[i];

        m_idx_q[i] = idx_q;
        m_idx_v[i] = idx_v;
        m_idx_vExtended[i] = idx_vExtended;
        ::pinocchio::setIndexes(joint, i, idx_q, idx_v, idx_vExtended);
        m_nqs[i] = ::pinocchio::nq(joint);
        m_nvs[i] = ::pinocchio::nv(joint);
        m_nvExtendeds[i] = ::pinocchio::nvExtended(joint);
        idx_q += m_nqs[i];
        idx_v += m_nvs[i];
        idx_vExtended += m_nvExtendeds[i];
      }
    }

    /// \brief Dimensions of the config and tangent space of the composite joint.
    int m_nq, m_nv, m_nvExtended;

    /// Keep information of both the dimension and the position of the joints in the composition.

    /// \brief Index in the config vector
    std::vector<int> m_idx_q;
    /// \brief Dimension of the segment in the config vector
    std::vector<int> m_nqs;
    /// \brief Index in the tangent vector
    std::vector<int> m_idx_v;
    /// \brief Dimension of the segment in the tangent vector
    std::vector<int> m_nvs;
    /// \brief Index in the jacobian matrix
    std::vector<int> m_idx_vExtended;
    /// \brief Dimension of the segment in the jacobian matrix
    std::vector<int> m_nvExtendeds;

  public:
    /// \brief Number of joints contained in the JointModelComposite
    int njoints;
  };

} // namespace pinocchio

#include <boost/type_traits.hpp>

namespace boost
{
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  struct has_nothrow_constructor<
    ::pinocchio::JointModelCompositeTpl<Scalar, Options, JointCollectionTpl>>
  : public integral_constant<bool, true>
  {
  };

  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  struct has_nothrow_copy<::pinocchio::JointModelCompositeTpl<Scalar, Options, JointCollectionTpl>>
  : public integral_constant<bool, true>
  {
  };

  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  struct has_nothrow_constructor<
    ::pinocchio::JointDataCompositeTpl<Scalar, Options, JointCollectionTpl>>
  : public integral_constant<bool, true>
  {
  };

  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  struct has_nothrow_copy<::pinocchio::JointDataCompositeTpl<Scalar, Options, JointCollectionTpl>>
  : public integral_constant<bool, true>
  {
  };
} // namespace boost

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
#include "pinocchio/multibody/joint/joint-composite.hxx"

#endif // ifndef __pinocchio_multibody_joint_composite_hpp__
