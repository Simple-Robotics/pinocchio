//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_serialization_constraints_set_hpp__
#define __pinocchio_serialization_constraints_set_hpp__

#include "pinocchio/algorithm/constraints/sets.hpp"
#include "pinocchio/serialization/eigen.hpp"

#include <boost/serialization/variant.hpp>

namespace boost
{
  namespace serialization
  {

    template<typename Archive, typename Derived>
    void serialize(Archive & ar, ::pinocchio::SetBase<Derived> & set, const unsigned int version)
    {
      // Nothing to do
      PINOCCHIO_UNUSED_VARIABLE(ar);
      PINOCCHIO_UNUSED_VARIABLE(set);
      PINOCCHIO_UNUSED_VARIABLE(version);
    }

    template<typename Archive, typename Derived>
    void
    serialize(Archive & ar, ::pinocchio::ConeBase<Derived> & set, const unsigned int /*version*/)
    {
      typedef ::pinocchio::ConeBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(set));
    }

    namespace internal
    {
      template<typename Scalar, int Options>
      struct BoxSetAccessor : public ::pinocchio::BoxSetTpl<Scalar, Options>
      {
        typedef ::pinocchio::BoxSetTpl<Scalar, Options> Base;
        using Base::m_lb;
        using Base::m_ub;
      };
    } // namespace internal

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar, ::pinocchio::BoxSetTpl<Scalar, Options> & set, const unsigned int /*version*/)
    {
      typedef ::pinocchio::BoxSetTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(set));

      typedef internal::BoxSetAccessor<Scalar, Options> Accessor;
      auto & set_ = reinterpret_cast<Accessor &>(set);
      ar & make_nvp("m_lb", set_.m_lb);
      ar & make_nvp("m_ub", set_.m_ub);
    }

    namespace internal
    {
      template<typename Scalar, int Options>
      struct UnboundedSetAccessor : public ::pinocchio::UnboundedSetTpl<Scalar, Options>
      {
        typedef ::pinocchio::UnboundedSetTpl<Scalar, Options> Base;
        using Base::m_size;
      };
    } // namespace internal

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::UnboundedSetTpl<Scalar, Options> & set,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::UnboundedSetTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(set));

      typedef internal::UnboundedSetAccessor<Scalar, Options> Accessor;
      auto & set_ = reinterpret_cast<Accessor &>(set);
      ar & make_nvp("m_size", set_.m_size);
    }

    namespace internal
    {
      template<typename Scalar, int Options>
      struct NullSetAccessor : public ::pinocchio::NullSetTpl<Scalar, Options>
      {
        typedef ::pinocchio::NullSetTpl<Scalar, Options> Base;
        using Base::m_size;
      };
    } // namespace internal

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar, ::pinocchio::NullSetTpl<Scalar, Options> & set, const unsigned int /*version*/)
    {
      typedef ::pinocchio::NullSetTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(set));

      typedef internal::NullSetAccessor<Scalar, Options> Accessor;
      auto & set_ = reinterpret_cast<Accessor &>(set);
      ar & make_nvp("m_size", set_.m_size);
    }

    template<typename Archive, typename Scalar>
    void serialize(
      Archive & ar,
      ::pinocchio::CoulombFrictionConeTpl<Scalar> & set,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::CoulombFrictionConeTpl<Scalar> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(set));

      ar & make_nvp("mu", set.mu);
    }

    template<typename Archive, typename Scalar>
    void serialize(
      Archive & ar,
      ::pinocchio::DualCoulombFrictionConeTpl<Scalar> & set,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::DualCoulombFrictionConeTpl<Scalar> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(set));

      ar & make_nvp("mu", set.mu);
    }

    namespace internal
    {
      template<typename Derived>
      struct OrthantConeBaseAccessor : public ::pinocchio::OrthantConeBase<Derived>
      {
        typedef ::pinocchio::OrthantConeBase<Derived> Base;
        using Base::m_size;
      };
    } // namespace internal

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar, ::pinocchio::OrthantConeBase<Derived> & set, const unsigned int /*version*/)
    {
      typedef ::pinocchio::OrthantConeBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(set));

      typedef internal::OrthantConeBaseAccessor<Derived> Accessor;
      auto & set_ = reinterpret_cast<Accessor &>(set);
      ar & make_nvp("m_size", set_.m_size);
    }

    template<typename Archive, typename Scalar>
    void serialize(
      Archive & ar,
      ::pinocchio::PositiveOrthantConeTpl<Scalar> & set,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::PositiveOrthantConeTpl<Scalar> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(set));
    }

    template<typename Archive, typename Scalar>
    void serialize(
      Archive & ar,
      ::pinocchio::NegativeOrthantConeTpl<Scalar> & set,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::NegativeOrthantConeTpl<Scalar> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(set));
    }

    namespace internal
    {
      template<typename Scalar>
      struct JointLimitConstraintConeAccessor
      : public ::pinocchio::JointLimitConstraintConeTpl<Scalar>
      {
        typedef ::pinocchio::JointLimitConstraintConeTpl<Scalar> Base;
        using Base::negative_orthant;
        using Base::positive_orthant;
      };
    } // namespace internal

    template<typename Archive, typename Scalar>
    void serialize(
      Archive & ar,
      ::pinocchio::JointLimitConstraintConeTpl<Scalar> & set,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::JointLimitConstraintConeTpl<Scalar> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(set));

      typedef internal::JointLimitConstraintConeAccessor<Scalar> Accessor;
      auto & set_ = reinterpret_cast<Accessor &>(set);
      ar & make_nvp("positive_orthant", set_.positive_orthant);
      ar & make_nvp("negative_orthant", set_.negative_orthant);
    }

  } // namespace serialization
} // namespace boost

#endif // __pinocchio_serialization_constraints_set_hpp__
