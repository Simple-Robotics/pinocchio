//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_serialization_constraints_data_hpp__
#define __pinocchio_serialization_constraints_data_hpp__

#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/serialization/eigen.hpp"
#include "pinocchio/serialization/se3.hpp"
#include "pinocchio/serialization/boost-blank.hpp"

#include <boost/serialization/variant.hpp>

namespace boost
{
  namespace serialization
  {

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar, ::pinocchio::ConstraintDataBase<Derived> & cdata, const unsigned int version)
    {
      PINOCCHIO_UNUSED_VARIABLE(ar);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      PINOCCHIO_UNUSED_VARIABLE(version);
    }

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::JointLimitConstraintDataTpl<Scalar, Options> & cdata,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::JointLimitConstraintDataTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cdata));

      ar & make_nvp("constraint_residual", cdata.constraint_residual);
    }

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::FrictionalJointConstraintDataTpl<Scalar, Options> & cdata,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::FrictionalJointConstraintDataTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cdata));
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::PointConstraintDataBase<Derived> & cdata,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::PointConstraintDataBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cdata));

      ar & make_nvp("constraint_force", cdata.constraint_force);
      ar & make_nvp("oMc1", cdata.oMc1);
      ar & make_nvp("oMc2", cdata.oMc2);
      ar & make_nvp("c1Mc2", cdata.c1Mc2);
      ar & make_nvp("constraint_position_error", cdata.constraint_position_error);
      ar & make_nvp("constraint_velocity_error", cdata.constraint_velocity_error);
      ar & make_nvp("constraint_acceleration_error", cdata.constraint_acceleration_error);
      ar & make_nvp("constraint_acceleration_biais_term", cdata.constraint_acceleration_biais_term);

      ar & make_nvp("A1_world", cdata.A1_world);
      ar & make_nvp("A2_world", cdata.A2_world);
      ar & make_nvp("A_world", cdata.A_world);
      ar & make_nvp("A1_local", cdata.A1_local);
      ar & make_nvp("A2_local", cdata.A2_local);
      ar & make_nvp("A_local", cdata.A_local);
    }

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::BilateralPointConstraintDataTpl<Scalar, Options> & cdata,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::BilateralPointConstraintDataTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cdata));
    }

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::FrictionalPointConstraintDataTpl<Scalar, Options> & cdata,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::FrictionalPointConstraintDataTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cdata));
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::FrameConstraintDataBase<Derived> & cdata,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::FrameConstraintDataBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cdata));

      ar & make_nvp("constraint_force", cdata.constraint_force);
      ar & make_nvp("oMc1", cdata.oMc1);
      ar & make_nvp("oMc2", cdata.oMc2);
      ar & make_nvp("c1Mc2", cdata.c1Mc2);
      ar & make_nvp("constraint_position_error", cdata.constraint_position_error);
      ar & make_nvp("constraint_velocity_error", cdata.constraint_velocity_error);
      ar & make_nvp("constraint_acceleration_error", cdata.constraint_acceleration_error);
      ar & make_nvp("constraint_acceleration_biais_term", cdata.constraint_acceleration_biais_term);

      ar & make_nvp("A1_world", cdata.A1_world);
      ar & make_nvp("A2_world", cdata.A2_world);
      ar & make_nvp("A_world", cdata.A_world);
      ar & make_nvp("A1_local", cdata.A1_local);
      ar & make_nvp("A2_local", cdata.A2_local);
      ar & make_nvp("A_local", cdata.A_local);
    }

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::WeldConstraintDataTpl<Scalar, Options> & cdata,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::WeldConstraintDataTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cdata));
    }

    template<
      typename Archive,
      typename Scalar,
      int Options,
      template<typename, int> class ConstraintCollectionTpl>
    void serialize(
      Archive & ar,
      pinocchio::ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cdata));

      typedef typename ConstraintCollectionTpl<Scalar, Options>::ConstraintDataVariant
        ConstraintDataVariant;
      ar & make_nvp("base_variant", base_object<ConstraintDataVariant>(cdata));
    }

  } // namespace serialization
} // namespace boost

#endif // __pinocchio_serialization_constraints_data_hpp__
