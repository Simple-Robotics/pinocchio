//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_serialization_constraints_model_hpp__
#define __pinocchio_serialization_constraints_model_hpp__

#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/serialization/eigen.hpp"
#include "pinocchio/serialization/se3.hpp"
#include "pinocchio/serialization/eigen-storage.hpp"
#include "pinocchio/serialization/constraints-set.hpp"
#include "pinocchio/serialization/boost-blank.hpp"

#include <boost/serialization/variant.hpp>

namespace boost
{
  namespace serialization
  {

    template<typename Archive, typename Scalar>
    void serialize(
      Archive & ar,
      ::pinocchio::BaumgarteCorrectorVectorParametersTpl<Scalar> & baumgarte_vector_parameters,
      const unsigned int /*version*/)
    {
      ar & make_nvp("Kp", baumgarte_vector_parameters.Kp);
      ar & make_nvp("Kd", baumgarte_vector_parameters.Kd);
    }

    template<typename Archive, typename Scalar>
    void serialize(
      Archive & ar,
      ::pinocchio::BaumgarteCorrectorParametersTpl<Scalar> & baumgarte_parameters,
      const unsigned int /*version*/)
    {
      ar & make_nvp("Kp", baumgarte_parameters.Kp);
      ar & make_nvp("Kd", baumgarte_parameters.Kd);
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar, ::pinocchio::ConstraintModelBase<Derived> & cmodel, const unsigned int version)
    {
      PINOCCHIO_UNUSED_VARIABLE(version);
      ar & make_nvp("name", cmodel.name);
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::KinematicsConstraintModelBase<Derived> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::KinematicsConstraintModelBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));

      ar & make_nvp("joint1_id", cmodel.joint1_id);
      ar & make_nvp("joint2_id", cmodel.joint2_id);
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::JointWiseConstraintModelBase<Derived> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::JointWiseConstraintModelBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
    }

    namespace internal
    {
      template<typename Derived>
      struct ConstraintModelCommonParametersAccessor
      : public ::pinocchio::ConstraintModelCommonParameters<Derived>
      {
        typedef ::pinocchio::ConstraintModelCommonParameters<Derived> Base;
        using Base::m_compliance;
        // CHOICE: right now we use the scalar Baumgarte
        // using Base::m_baumgarte_vector_parameters;
        using Base::m_baumgarte_parameters;
      };
    } // namespace internal

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::ConstraintModelCommonParameters<Derived> & cmodel,
      const unsigned int version)
    {
      PINOCCHIO_UNUSED_VARIABLE(version);
      typedef internal::ConstraintModelCommonParametersAccessor<Derived> Accessor;
      auto & cmodel_ = reinterpret_cast<Accessor &>(cmodel);
      ar & make_nvp("m_compliance", cmodel_.m_compliance);
      // CHOICE: right now we use the scalar Baumgarte
      // ar & make_nvp("m_baumgarte_vector_parameters", cmodel_.m_baumgarte_vector_parameters);
      ar & make_nvp("m_baumgarte_parameters", cmodel_.m_baumgarte_parameters);
    }

    namespace internal
    {
      template<typename Scalar, int Options>
      struct JointLimitConstraintModelAccessor
      : public ::pinocchio::JointLimitConstraintModelTpl<Scalar, Options>
      {
        typedef ::pinocchio::JointLimitConstraintModelTpl<Scalar, Options> Base;
        using Base::activable_idx_qs;
        using Base::activable_idx_qs_reduce;
        using Base::activable_idx_rows;
        using Base::activable_idx_vs;
        using Base::activable_joints;
        using Base::activable_nvs;
        using Base::active_compliance;
        using Base::active_compliance_storage;
        using Base::active_idx_qs_reduce;
        using Base::active_idx_rows;
        using Base::active_idx_vs;
        using Base::active_nvs;
        using Base::active_set_indexes;
        using Base::bound_position_limit;
        using Base::bound_position_margin;
        using Base::lower_activable_size;
        using Base::lower_active_size;
        using Base::m_set;
        using Base::nq_reduce;
        using Base::nv_max_atom;
        using Base::row_indexes;
        using Base::row_sparsity_pattern;
      };
    } // namespace internal

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::JointLimitConstraintModelTpl<Scalar, Options> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::JointLimitConstraintModelTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
      typedef typename Self::BaseCommonParameters BaseCommonParameters;
      ar & make_nvp(
        "base_common_parameters", boost::serialization::base_object<BaseCommonParameters>(cmodel));

      typedef internal::JointLimitConstraintModelAccessor<Scalar, Options> Accessor;
      auto & cmodel_ = reinterpret_cast<Accessor &>(cmodel);
      ar & make_nvp("m_set", cmodel_.m_set);
      ar & make_nvp("activable_joints", cmodel_.activable_joints);
      ar & make_nvp("nq_reduce", cmodel_.nq_reduce);
      ar & make_nvp("nv_max_atom", cmodel_.nv_max_atom);
      ar & make_nvp("lower_activable_size", cmodel_.lower_activable_size);
      ar & make_nvp("lower_active_size", cmodel_.lower_active_size);
      ar & make_nvp("row_sparsity_pattern", cmodel_.row_sparsity_pattern);
      ar & make_nvp("row_indexes", cmodel_.row_indexes);
      ar & make_nvp("bound_position_limit", cmodel_.bound_position_limit);
      ar & make_nvp("bound_position_margin", cmodel_.bound_position_margin);
      ar & make_nvp("activable_idx_qs", cmodel_.activable_idx_qs);
      ar & make_nvp("active_set_indexes", cmodel_.active_set_indexes);
      ar & make_nvp("activable_idx_rows", cmodel_.activable_idx_rows);
      ar & make_nvp("activable_idx_qs_reduce", cmodel_.activable_idx_qs_reduce);
      ar & make_nvp("activable_nvs", cmodel_.activable_nvs);
      ar & make_nvp("activable_idx_vs", cmodel_.activable_idx_vs);
      ar & make_nvp("active_idx_rows", cmodel_.active_idx_rows);
      ar & make_nvp("active_idx_qs_reduce", cmodel_.active_idx_qs_reduce);
      ar & make_nvp("active_nvs", cmodel_.active_nvs);
      ar & make_nvp("active_idx_vs", cmodel_.active_idx_vs);
      ar & make_nvp("m_set", cmodel_.m_set);
      ar & make_nvp("active_compliance_storage", cmodel_.active_compliance_storage);

      if (Archive::is_loading::value)
      {
        cmodel_.active_compliance = cmodel_.active_compliance_storage.map();
      }
    }

    namespace internal
    {
      template<typename Scalar, int Options>
      struct FrictionalJointConstraintModelAccessor
      : public ::pinocchio::FrictionalJointConstraintModelTpl<Scalar, Options>
      {
        typedef ::pinocchio::FrictionalJointConstraintModelTpl<Scalar, Options> Base;
        using Base::active_dofs;
        using Base::active_joints;
        using Base::m_set;
        using Base::row_active_indexes;
        using Base::row_sparsity_pattern;
      };
    } // namespace internal

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::FrictionalJointConstraintModelTpl<Scalar, Options> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::FrictionalJointConstraintModelTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
      typedef typename Self::BaseCommonParameters BaseCommonParameters;
      ar & make_nvp(
        "base_common_parameters", boost::serialization::base_object<BaseCommonParameters>(cmodel));

      typedef internal::FrictionalJointConstraintModelAccessor<Scalar, Options> Accessor;
      auto & cmodel_ = reinterpret_cast<Accessor &>(cmodel);
      ar & make_nvp("active_joints", cmodel_.active_joints);
      ar & make_nvp("active_dofs", cmodel_.active_dofs);
      ar & make_nvp("row_sparsity_pattern", cmodel_.row_sparsity_pattern);
      ar & make_nvp("row_active_indexes", cmodel_.row_active_indexes);
      ar & make_nvp("m_set", cmodel_.m_set);
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::PointConstraintModelBase<Derived> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::PointConstraintModelBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
      typedef typename Self::BaseCommonParameters BaseCommonParameters;
      ar & make_nvp(
        "base_common_parameters", boost::serialization::base_object<BaseCommonParameters>(cmodel));

      ar & make_nvp("joint1_placement", cmodel.joint1_placement);
      ar & make_nvp("joint2_placement", cmodel.joint2_placement);
      ar & make_nvp("desired_constraint_offset", cmodel.desired_constraint_offset);
      ar & make_nvp("desired_constraint_velocity", cmodel.desired_constraint_velocity);
      ar & make_nvp("desired_constraint_acceleration", cmodel.desired_constraint_acceleration);
      ar & make_nvp("colwise_joint1_sparsity", cmodel.colwise_joint1_sparsity);
      ar & make_nvp("colwise_joint2_sparsity", cmodel.colwise_joint2_sparsity);
      ar & make_nvp("joint1_span_indexes", cmodel.joint1_span_indexes);
      ar & make_nvp("joint2_span_indexes", cmodel.joint2_span_indexes);
      ar & make_nvp("loop_span_indexes", cmodel.loop_span_indexes);
      ar & make_nvp("colwise_sparsity", cmodel.colwise_sparsity);
      ar & make_nvp("colwise_span_indexes", cmodel.colwise_span_indexes);
      ar & make_nvp("nv", cmodel.nv);
      ar & make_nvp("depth_joint1", cmodel.depth_joint1);
      ar & make_nvp("depth_joint2", cmodel.depth_joint2);
    }

    namespace internal
    {
      template<typename Scalar, int Options>
      struct BilateralPointConstraintModelAccessor
      : public ::pinocchio::BilateralPointConstraintModelTpl<Scalar, Options>
      {
        typedef ::pinocchio::BilateralPointConstraintModelTpl<Scalar, Options> Base;
        using Base::m_set;
      };
    } // namespace internal

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::BilateralPointConstraintModelTpl<Scalar, Options> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::BilateralPointConstraintModelTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));

      typedef internal::BilateralPointConstraintModelAccessor<Scalar, Options> Accessor;
      auto & cmodel_ = reinterpret_cast<Accessor &>(cmodel);
      ar & make_nvp("m_set", cmodel_.m_set);
    }

    namespace internal
    {
      template<typename Scalar, int Options>
      struct FrictionalPointConstraintModelAccessor
      : public ::pinocchio::FrictionalPointConstraintModelTpl<Scalar, Options>
      {
        typedef ::pinocchio::FrictionalPointConstraintModelTpl<Scalar, Options> Base;
        using Base::m_set;
      };
    } // namespace internal

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::FrictionalPointConstraintModelTpl<Scalar, Options> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::FrictionalPointConstraintModelTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));

      typedef internal::FrictionalPointConstraintModelAccessor<Scalar, Options> Accessor;
      auto & cmodel_ = reinterpret_cast<Accessor &>(cmodel);
      ar & make_nvp("m_set", cmodel_.m_set);
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::FrameConstraintModelBase<Derived> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::FrameConstraintModelBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
      typedef typename Self::BaseCommonParameters BaseCommonParameters;
      ar & make_nvp(
        "base_common_parameters", boost::serialization::base_object<BaseCommonParameters>(cmodel));

      // TODO: point/frame constraint models data structure are identical, factor them
      ar & make_nvp("joint1_placement", cmodel.joint1_placement);
      ar & make_nvp("joint2_placement", cmodel.joint2_placement);
      ar & make_nvp("desired_constraint_offset", cmodel.desired_constraint_offset);
      ar & make_nvp("desired_constraint_velocity", cmodel.desired_constraint_velocity);
      ar & make_nvp("desired_constraint_acceleration", cmodel.desired_constraint_acceleration);
      ar & make_nvp("colwise_joint1_sparsity", cmodel.colwise_joint1_sparsity);
      ar & make_nvp("colwise_joint2_sparsity", cmodel.colwise_joint2_sparsity);
      ar & make_nvp("joint1_span_indexes", cmodel.joint1_span_indexes);
      ar & make_nvp("joint2_span_indexes", cmodel.joint2_span_indexes);
      ar & make_nvp("loop_span_indexes", cmodel.loop_span_indexes);
      ar & make_nvp("colwise_sparsity", cmodel.colwise_sparsity);
      ar & make_nvp("colwise_span_indexes", cmodel.colwise_span_indexes);
      ar & make_nvp("nv", cmodel.nv);
      ar & make_nvp("depth_joint1", cmodel.depth_joint1);
      ar & make_nvp("depth_joint2", cmodel.depth_joint2);
    }

    namespace internal
    {
      template<typename Scalar, int Options>
      struct WeldConstraintModelAccessor
      : public ::pinocchio::WeldConstraintModelTpl<Scalar, Options>
      {
        typedef ::pinocchio::WeldConstraintModelTpl<Scalar, Options> Base;
        using Base::m_set;
      };
    } // namespace internal

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::WeldConstraintModelTpl<Scalar, Options> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::WeldConstraintModelTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));

      typedef internal::WeldConstraintModelAccessor<Scalar, Options> Accessor;
      auto & cmodel_ = reinterpret_cast<Accessor &>(cmodel);
      ar & make_nvp("m_set", cmodel_.m_set);
    }

    template<
      typename Archive,
      typename Scalar,
      int Options,
      template<typename, int> class ConstraintCollectionTpl>
    void serialize(
      Archive & ar,
      pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));

      typedef typename ConstraintCollectionTpl<Scalar, Options>::ConstraintModelVariant
        ConstraintModelVariant;
      ar & make_nvp(
        "base_variant", boost::serialization::base_object<ConstraintModelVariant>(cmodel));
    }

  } // namespace serialization
} // namespace boost

#endif // __pinocchio_serialization_constraints_model_hpp__
