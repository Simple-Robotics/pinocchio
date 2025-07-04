//
// Copyright (c) 2015-2021 CNRS INRIA
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_multibody_model_hxx__
#define __pinocchio_multibody_model_hxx__

#include "pinocchio/utils/string-generator.hpp"
#include "pinocchio/multibody/liegroup/liegroup-algo.hpp"
#include "pinocchio/algorithm/model.hpp"

/// @cond DEV

namespace pinocchio
{
  namespace details
  {
    struct FilterFrame
    {
      const std::string & name;
      const FrameType & typeMask;

      FilterFrame(const std::string & name, const FrameType & typeMask)
      : name(name)
      , typeMask(typeMask)
      {
      }

      template<typename Scalar, int Options>
      bool operator()(const FrameTpl<Scalar, Options> & frame) const
      {
        return (typeMask & frame.type) && (name == frame.name);
      }
    };
  } // namespace details

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  const typename ModelTpl<Scalar, Options, JointCollectionTpl>::Vector3
    ModelTpl<Scalar, Options, JointCollectionTpl>::gravity981((Scalar)0, (Scalar)0, (Scalar)-9.81);

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline std::ostream &
  operator<<(std::ostream & os, const ModelTpl<Scalar, Options, JointCollectionTpl> & model)
  {
    typedef typename ModelTpl<Scalar, Options, JointCollectionTpl>::Index Index;

    os << "Nb joints = " << model.njoints << " (nq=" << model.nq << ",nv=" << model.nv << ")"
       << std::endl;
    for (Index i = 0; i < (Index)(model.njoints); ++i)
    {
      os << "  Joint " << i << " " << model.names[i] << ": parent=" << model.parents[i]
         << std::endl;
    }

    return os;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & min_effort,
    const VectorXs & max_effort,
    const VectorXs & min_velocity,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config,
    const VectorXs & min_joint_friction,
    const VectorXs & max_joint_friction,
    const VectorXs & joint_damping)
  {
    const VectorXs config_limit_margin =
      VectorXs::Constant(joint_model.nq(), static_cast<Scalar>(0));
    return addJoint(
      parent, joint_model, joint_placement, joint_name, min_effort, max_effort, min_velocity,
      max_velocity, min_config, max_config, config_limit_margin, min_joint_friction,
      max_joint_friction, joint_damping);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & min_effort,
    const VectorXs & max_effort,
    const VectorXs & min_velocity,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config,
    const VectorXs & config_limit_margin,
    const VectorXs & min_joint_friction,
    const VectorXs & max_joint_friction,
    const VectorXs & joint_damping)
  {
    assert(
      (njoints == (int)joints.size()) && (njoints == (int)inertias.size())
      && (njoints == (int)parents.size()) && (njoints == (int)jointPlacements.size()));
    assert((joint_model.nq() >= 0) && (joint_model.nv() >= 0) && (joint_model.nvExtended() >= 0));
    assert(joint_model.nq() >= joint_model.nv());

    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      min_effort.size(), joint_model.nv(), "The joint minimal effort vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      min_joint_friction.size(), joint_model.nv(),
      "The joint minimal dry friction vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      min_velocity.size(), joint_model.nv(),
      "The joint minimal velocity vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      max_effort.size(), joint_model.nv(), "The joint maximum effort vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      max_joint_friction.size(), joint_model.nv(),
      "The joint maximum dry friction vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      max_velocity.size(), joint_model.nv(),
      "The joint maximum velocity vector is not of right size");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      (min_effort.array() <= max_effort.array()).all(),
      "Some components of min_effort are greater than max_effort");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      (min_joint_friction.array() <= max_joint_friction.array()).all(),
      "Some components of min_dry_friction are greater than max_dry_friction");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      (min_velocity.array() <= max_velocity.array()).all(),
      "Some components of min_velocity are greater than max_velocity");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      min_config.size(), joint_model.nq(),
      "The joint lower configuration bound is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      max_config.size(), joint_model.nq(),
      "The joint upper configuration bound is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      config_limit_margin.size(), joint_model.nq(),
      "The joint config limit margin is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      joint_damping.size(), joint_model.nv(), "The joint damping vector is not of right size");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      parent < (JointIndex)njoints, "The index of the parent joint is not valid.");

    JointIndex joint_id = (JointIndex)(njoints++);

    joints.push_back(JointModel(joint_model.derived()));
    JointModel & jmodel = joints.back();
    jmodel.setIndexes(joint_id, nq, nv, nvExtended);

    const int joint_nq = jmodel.nq();
    const int joint_idx_q = jmodel.idx_q();
    const int joint_nv = jmodel.nv();
    const int joint_idx_v = jmodel.idx_v();
    const int joint_nvExtended = jmodel.nvExtended();
    const int joint_idx_vExtended = jmodel.idx_vExtended();

    assert(joint_idx_q >= 0);
    assert(joint_idx_v >= 0);
    assert(joint_idx_vExtended >= 0);

    inertias.push_back(Inertia::Zero());
    parents.push_back(parent);
    children.push_back(IndexVector());
    children[parent].push_back(joint_id);
    jointPlacements.push_back(joint_placement);
    names.push_back(joint_name);

    nq += joint_nq;
    nqs.push_back(joint_nq);
    idx_qs.push_back(joint_idx_q);
    nv += joint_nv;
    nvs.push_back(joint_nv);
    idx_vs.push_back(joint_idx_v);
    nvExtended += joint_nvExtended;
    nvExtendeds.push_back(joint_nvExtended);
    idx_vExtendeds.push_back(joint_idx_vExtended);

    if (joint_nq > 0 && joint_nv > 0)
    {
      upperEffortLimit.conservativeResize(nv);
      jmodel.jointVelocitySelector(upperEffortLimit) = max_effort;
      lowerEffortLimit.conservativeResize(nv);
      jmodel.jointVelocitySelector(lowerEffortLimit) = min_effort;
      upperVelocityLimit.conservativeResize(nv);
      jmodel.jointVelocitySelector(upperVelocityLimit) = max_velocity;
      lowerVelocityLimit.conservativeResize(nv);
      jmodel.jointVelocitySelector(lowerVelocityLimit) = max_velocity;
      lowerPositionLimit.conservativeResize(nq);
      jmodel.jointConfigSelector(lowerPositionLimit) = min_config;
      upperPositionLimit.conservativeResize(nq);
      jmodel.jointConfigSelector(upperPositionLimit) = max_config;
      positionLimitMargin.conservativeResize(nq);
      jmodel.jointConfigSelector(positionLimitMargin) = config_limit_margin;

      armature.conservativeResize(nv);
      jmodel.jointVelocitySelector(armature).setZero();
      rotorInertia.conservativeResize(nv);
      jmodel.jointVelocitySelector(rotorInertia).setZero();
      rotorGearRatio.conservativeResize(nv);
      jmodel.jointVelocitySelector(rotorGearRatio).setOnes();
      upperDryFrictionLimit.conservativeResize(nv);
      jmodel.jointVelocitySelector(upperDryFrictionLimit) = max_joint_friction;
      lowerDryFrictionLimit.conservativeResize(nv);
      jmodel.jointVelocitySelector(lowerDryFrictionLimit) = min_joint_friction;
      damping.conservativeResize(nv);
      jmodel.jointVelocitySelector(damping) = joint_damping;
    }

    // Init and add joint index to its parent subtrees.
    subtrees.push_back(IndexVector(1));
    subtrees[joint_id][0] = joint_id;
    addJointIndexToParentSubtrees(joint_id);

    // Init and add joint index to the supports
    supports.push_back(supports[parent]);
    supports[joint_id].push_back(joint_id);

    mimic_joint_supports.push_back(mimic_joint_supports[parent]);
    if (
      const auto & jmodel_ =
        boost::get<JointModelMimicTpl<Scalar, Options, JointCollectionTpl>>(&jmodel))
    {
      mimicking_joints.push_back(jmodel.id());
      mimicked_joints.push_back(jmodel_->jmodel().id());
      mimic_joint_supports[joint_id].push_back(joint_id);
    }
    return joint_id;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & max_effort,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config,
    const VectorXs & config_limit_margin,
    const VectorXs & friction,
    const VectorXs & damping)
  {

    return addJoint(
      parent, joint_model, joint_placement, joint_name, -max_effort, max_effort, -max_velocity,
      max_velocity, min_config, max_config, config_limit_margin, -friction, friction, damping);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & max_effort,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config,
    const VectorXs & friction,
    const VectorXs & damping)
  {
    const VectorXs config_limit_margin =
      VectorXs::Constant(joint_model.nq(), static_cast<Scalar>(0));

    return addJoint(
      parent, joint_model, joint_placement, joint_name, -max_effort, max_effort, -max_velocity,
      max_velocity, min_config, max_config, config_limit_margin, -friction, friction, damping);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & max_effort,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config)
  {
    const VectorXs config_limit_margin =
      VectorXs::Constant(joint_model.nq(), static_cast<Scalar>(0));
    const VectorXs friction = VectorXs::Constant(joint_model.nv(), static_cast<Scalar>(0));
    const VectorXs damping = VectorXs::Constant(joint_model.nv(), static_cast<Scalar>(0));

    return addJoint(
      parent, joint_model, joint_placement, joint_name, max_effort, max_velocity, min_config,
      max_config, config_limit_margin, friction, damping);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name,
    const VectorXs & max_effort,
    const VectorXs & max_velocity,
    const VectorXs & min_config,
    const VectorXs & max_config,
    const VectorXs & config_limit_margin)
  {
    const VectorXs friction = VectorXs::Constant(joint_model.nv(), static_cast<Scalar>(0));
    const VectorXs damping = VectorXs::Constant(joint_model.nv(), static_cast<Scalar>(0));

    return addJoint(
      parent, joint_model, joint_placement, joint_name, max_effort, max_velocity, min_config,
      max_config, config_limit_margin, friction, damping);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addJoint(
    const JointIndex parent,
    const JointModel & joint_model,
    const SE3 & joint_placement,
    const std::string & joint_name)
  {
    const VectorXs max_effort =
      VectorXs::Constant(joint_model.nv(), std::numeric_limits<Scalar>::max());
    const VectorXs max_velocity =
      VectorXs::Constant(joint_model.nv(), std::numeric_limits<Scalar>::max());
    const VectorXs min_config =
      VectorXs::Constant(joint_model.nq(), -std::numeric_limits<Scalar>::max());
    const VectorXs max_config =
      VectorXs::Constant(joint_model.nq(), std::numeric_limits<Scalar>::max());
    const VectorXs config_limit_margin =
      VectorXs::Constant(joint_model.nq(), static_cast<Scalar>(0));

    return addJoint(
      parent, joint_model, joint_placement, joint_name, max_effort, max_velocity, min_config,
      max_config, config_limit_margin);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  FrameIndex ModelTpl<Scalar, Options, JointCollectionTpl>::addJointFrame(
    const JointIndex joint_index, int previous_frame_index)
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      joint_index < joints.size(),
      "The joint index is larger than the number of joints in the model.");
    if (previous_frame_index < 0)
    {
      // FIXED_JOINT is required because the parent can be the universe and its
      // type is FIXED_JOINT
      previous_frame_index =
        (int)getFrameId(names[parents[joint_index]], (FrameType)(JOINT | FIXED_JOINT));
    }
    assert((size_t)previous_frame_index < frames.size() && "Frame index out of bound");

    // Add a the joint frame attached to itself to the frame vector - redundant information but
    // useful.
    return addFrame(Frame(
      names[joint_index], joint_index, (FrameIndex)previous_frame_index, SE3::Identity(), JOINT));
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  template<typename NewScalar>
  typename CastType<NewScalar, ModelTpl<Scalar, Options, JointCollectionTpl>>::type
  ModelTpl<Scalar, Options, JointCollectionTpl>::cast() const
  {
    typedef ModelTpl<NewScalar, Options, JointCollectionTpl> ReturnType;

    ReturnType res;

    res.nq = nq;
    res.nv = nv;
    res.nvExtended = nvExtended;
    res.njoints = njoints;
    res.nbodies = nbodies;
    res.nframes = nframes;
    res.parents = parents;
    res.children = children;
    res.names = names;
    res.subtrees = subtrees;
    res.supports = supports;
    res.mimic_joint_supports = mimic_joint_supports;
    res.mimicking_joints = mimicking_joints;
    res.mimicked_joints = mimicked_joints;
    res.gravity = gravity.template cast<NewScalar>();
    res.name = name;

    res.idx_qs = idx_qs;
    res.nqs = nqs;
    res.idx_vs = idx_vs;
    res.nvs = nvs;
    res.idx_vExtendeds = idx_vExtendeds;
    res.nvExtendeds = nvExtendeds;
    // Eigen Vectors
    res.armature = armature.template cast<NewScalar>();
    res.damping = damping.template cast<NewScalar>();
    res.rotorInertia = rotorInertia.template cast<NewScalar>();
    res.rotorGearRatio = rotorGearRatio.template cast<NewScalar>();
    res.upperEffortLimit = upperEffortLimit.template cast<NewScalar>();
    res.lowerEffortLimit = lowerEffortLimit.template cast<NewScalar>();
    res.upperDryFrictionLimit = upperDryFrictionLimit.template cast<NewScalar>();
    res.lowerDryFrictionLimit = lowerDryFrictionLimit.template cast<NewScalar>();
    res.lowerVelocityLimit = lowerVelocityLimit.template cast<NewScalar>();
    res.upperVelocityLimit = upperVelocityLimit.template cast<NewScalar>();
    res.lowerPositionLimit = lowerPositionLimit.template cast<NewScalar>();
    res.upperPositionLimit = upperPositionLimit.template cast<NewScalar>();
    res.positionLimitMargin = positionLimitMargin.template cast<NewScalar>();

    typename ConfigVectorMap::const_iterator it;
    for (it = referenceConfigurations.begin(); it != referenceConfigurations.end(); it++)
    {
      res.referenceConfigurations.insert(
        std::make_pair(it->first, it->second.template cast<NewScalar>()));
    }

    // reserve vectors
    res.inertias.resize(inertias.size());
    res.jointPlacements.resize(jointPlacements.size());
    res.joints.resize(joints.size());

    // copy into vectors
    for (size_t k = 0; k < joints.size(); ++k)
    {
      res.inertias[k] = inertias[k].template cast<NewScalar>();
      res.jointPlacements[k] = jointPlacements[k].template cast<NewScalar>();
      res.joints[k] = joints[k].template cast<NewScalar>();
    }

    res.frames.resize(frames.size());
    for (size_t k = 0; k < frames.size(); ++k)
    {
      res.frames[k] = frames[k].template cast<NewScalar>();
    }

    return res;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  template<template<typename, int> class OtherJointCollectionTpl>
  ModelTpl<Scalar, Options, JointCollectionTpl> &
  ModelTpl<Scalar, Options, JointCollectionTpl>::operator=(
    const ModelTpl<Scalar, Options, OtherJointCollectionTpl> & other)
  {
    this->nq = other.nq;
    this->nv = other.nv;
    this->nvExtended = other.nvExtended;
    this->njoints = other.njoints;
    this->nbodies = other.nbodies;
    this->nframes = other.nframes;
    this->inertias = other.inertias;
    this->jointPlacements = other.jointPlacements;
    this->joints.clear();
    this->joints.reserve(other.joints.size());
    for (const auto & other_joint : other.joints)
    {
      this->joints.push_back(other_joint);
    }
    this->idx_qs = other.idx_qs;
    this->nqs = other.nqs;
    this->idx_vs = other.idx_vs;
    this->nvs = other.nvs;
    this->idx_vExtendeds = other.idx_vExtendeds;
    this->nvExtendeds = other.nvExtendeds;
    this->parents = other.parents;
    this->children = other.children;
    this->names = other.names;
    this->referenceConfigurations = other.referenceConfigurations;
    this->armature = other.armature;
    this->rotorInertia = other.rotorInertia;
    this->rotorGearRatio = other.rotorGearRatio;
    this->lowerDryFrictionLimit = other.lowerDryFrictionLimit;
    this->upperDryFrictionLimit = other.upperDryFrictionLimit;
    this->damping = other.damping;
    this->lowerEffortLimit = other.lowerEffortLimit;
    this->upperEffortLimit = other.upperEffortLimit;
    this->lowerVelocityLimit = other.lowerVelocityLimit;
    this->upperVelocityLimit = other.upperVelocityLimit;
    this->lowerPositionLimit = other.lowerPositionLimit;
    this->upperPositionLimit = other.upperPositionLimit;
    this->positionLimitMargin = other.positionLimitMargin;
    this->frames = other.frames;
    this->supports = other.supports;
    this->subtrees = other.subtrees;
    this->mimic_joint_supports = other.mimic_joint_supports;
    this->mimicking_joints = other.mimicking_joints;
    this->mimicked_joints = other.mimicked_joints;
    this->gravity = other.gravity;
    this->name = other.name;
    return *this;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  bool ModelTpl<Scalar, Options, JointCollectionTpl>::operator==(const ModelTpl & other) const
  {
    bool res = other.nq == nq && other.nv == nv && other.nvExtended == nvExtended
               && other.njoints == njoints && other.nbodies == nbodies && other.nframes == nframes
               && other.parents == parents && other.children == children && other.names == names
               && other.subtrees == subtrees && other.mimic_joint_supports == mimic_joint_supports
               && other.mimicking_joints == mimicking_joints
               && other.mimicked_joints == mimicked_joints && other.gravity == gravity
               && other.name == name;

    res &= other.idx_qs == idx_qs && other.nqs == nqs && other.idx_vs == idx_vs && other.nvs == nvs
           && other.idx_vExtendeds == idx_vExtendeds && other.nvExtendeds == nvExtendeds;

    if (other.referenceConfigurations.size() != referenceConfigurations.size())
      return false;

    typename ConfigVectorMap::const_iterator it = referenceConfigurations.begin();
    typename ConfigVectorMap::const_iterator it_other = other.referenceConfigurations.begin();
    for (long k = 0; k < (long)referenceConfigurations.size(); ++k)
    {
      if (it->second.size() != it_other->second.size())
        return false;
      if (it->second != it_other->second)
        return false;
      std::advance(it, 1);
      std::advance(it_other, 1);
    }
    if (other.armature.size() != armature.size())
      return false;
    res &= other.armature == armature;
    if (!res)
      return res;

    if (other.damping.size() != damping.size())
      return false;
    res &= other.damping == damping;
    if (!res)
      return res;

    if (other.rotorInertia.size() != rotorInertia.size())
      return false;
    res &= other.rotorInertia == rotorInertia;
    if (!res)
      return res;

    if (other.rotorGearRatio.size() != rotorGearRatio.size())
      return false;
    res &= other.rotorGearRatio == rotorGearRatio;
    if (!res)
      return res;

    if (other.lowerEffortLimit.size() != lowerEffortLimit.size())
      return false;
    res &= other.lowerEffortLimit == lowerEffortLimit;
    if (!res)
      return res;

    if (other.upperEffortLimit.size() != upperEffortLimit.size())
      return false;
    res &= other.upperEffortLimit == upperEffortLimit;
    if (!res)
      return res;

    if (other.lowerDryFrictionLimit.size() != lowerDryFrictionLimit.size())
      return false;
    res &= other.lowerDryFrictionLimit == lowerDryFrictionLimit;
    if (!res)
      return res;

    if (other.upperDryFrictionLimit.size() != upperDryFrictionLimit.size())
      return false;
    res &= other.upperDryFrictionLimit == upperDryFrictionLimit;
    if (!res)
      return res;

    if (other.lowerVelocityLimit.size() != lowerVelocityLimit.size())
      return false;
    res &= other.lowerVelocityLimit == lowerVelocityLimit;
    if (!res)
      return res;

    if (other.upperVelocityLimit.size() != upperVelocityLimit.size())
      return false;
    res &= other.upperVelocityLimit == upperVelocityLimit;
    if (!res)
      return res;

    if (other.lowerPositionLimit.size() != lowerPositionLimit.size())
      return false;
    res &= other.lowerPositionLimit == lowerPositionLimit;
    if (!res)
      return res;

    if (other.upperPositionLimit.size() != upperPositionLimit.size())
      return false;
    res &= other.upperPositionLimit == upperPositionLimit;

    if (other.positionLimitMargin.size() != positionLimitMargin.size())
      return false;
    res &= other.positionLimitMargin == positionLimitMargin;

    if (!res)
      return res;

    for (size_t k = 1; k < inertias.size(); ++k)
    {
      res &= other.inertias[k] == inertias[k];
      if (!res)
        return res;
    }

    for (size_t k = 1; k < other.jointPlacements.size(); ++k)
    {
      res &= other.jointPlacements[k] == jointPlacements[k];
      if (!res)
        return res;
    }

    res &= other.joints == joints && other.frames == frames;

    return res;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  void ModelTpl<Scalar, Options, JointCollectionTpl>::appendBodyToJoint(
    const typename ModelTpl::JointIndex joint_index, const Inertia & Y, const SE3 & body_placement)
  {
    const Inertia & iYf = Y.se3Action(body_placement);
    inertias[joint_index] += iYf;
    nbodies++;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::FrameIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addBodyFrame(
    const std::string & body_name,
    const JointIndex & parentJoint,
    const SE3 & body_placement,
    int parentFrame)
  {
    if (parentFrame < 0)
    {
      // FIXED_JOINT is required because the parent can be the universe and its
      // type is FIXED_JOINT
      parentFrame = (int)getFrameId(names[parentJoint], (FrameType)(JOINT | FIXED_JOINT));
    }
    PINOCCHIO_CHECK_INPUT_ARGUMENT((size_t)parentFrame < frames.size(), "Frame index out of bound");
    return addFrame(Frame(body_name, parentJoint, (FrameIndex)parentFrame, body_placement, BODY));
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline typename ModelTpl<Scalar, Options, JointCollectionTpl>::FrameIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::getBodyId(const std::string & name) const
  {
    return getFrameId(name, BODY);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline bool
  ModelTpl<Scalar, Options, JointCollectionTpl>::existBodyName(const std::string & name) const
  {
    return existFrame(name, BODY);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::getJointId(const std::string & name) const
  {
    typedef std::vector<std::string>::iterator::difference_type it_diff_t;
    it_diff_t res = std::find(names.begin(), names.end(), name) - names.begin();
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      (res < INT_MAX), "Id superior to int range. Should never happen.");
    return ModelTpl::JointIndex(res);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline bool
  ModelTpl<Scalar, Options, JointCollectionTpl>::existJointName(const std::string & name) const
  {
    return (names.end() != std::find(names.begin(), names.end(), name));
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline typename ModelTpl<Scalar, Options, JointCollectionTpl>::FrameIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::getFrameId(
    const std::string & name, const FrameType & type) const
  {
    typename PINOCCHIO_ALIGNED_STD_VECTOR(Frame)::const_iterator it =
      std::find_if(frames.begin(), frames.end(), details::FilterFrame(name, type));
    std::ostringstream os;
    os << "Several frames match the filter - please specify the FrameType (name=\"" << name
       << "\", type=\"" << type << "\")";
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      ((it == frames.end()
        || (std::find_if(boost::next(it), frames.end(), details::FilterFrame(name, type)) == frames.end()))),
      os.str().c_str());
    return FrameIndex(it - frames.begin());
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  inline bool ModelTpl<Scalar, Options, JointCollectionTpl>::existFrame(
    const std::string & name, const FrameType & type) const
  {
    return std::find_if(frames.begin(), frames.end(), details::FilterFrame(name, type))
           != frames.end();
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  typename ModelTpl<Scalar, Options, JointCollectionTpl>::FrameIndex
  ModelTpl<Scalar, Options, JointCollectionTpl>::addFrame(
    const Frame & frame, const bool append_inertia)
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      frame.parentJoint < (JointIndex)njoints, "The index of the parent joint is not valid.");

    //    TODO: fix it
    //    PINOCCHIO_CHECK_INPUT_ARGUMENT(frame.inertia.isValid(),
    //                                   "The input inertia is not valid.")

    // Check if the frame.name exists with the same type
    if (existFrame(frame.name, frame.type))
    {
      return getFrameId(frame.name, frame.type);
    }
    // else: we must add a new frames to the current stack
    frames.push_back(frame);
    if (append_inertia)
      inertias[frame.parentJoint] += frame.placement.act(frame.inertia);
    nframes++;
    return FrameIndex(nframes - 1);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  void ModelTpl<Scalar, Options, JointCollectionTpl>::addJointIndexToParentSubtrees(
    const JointIndex joint_id)
  {
    for (JointIndex parent = parents[joint_id]; parent > 0; parent = parents[parent])
      subtrees[parent].push_back(joint_id);

    // Also add joint_id to the universe
    subtrees[0].push_back(joint_id);
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  std::vector<bool> ModelTpl<Scalar, Options, JointCollectionTpl>::hasConfigurationLimit() const
  {
    std::vector<bool> vec;
    for (Index i = 1; i < (Index)(njoints); ++i)
    {
      const std::vector<bool> & cf_limits = joints[i].hasConfigurationLimit();
      vec.insert(vec.end(), cf_limits.begin(), cf_limits.end());
    }
    return vec;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  std::vector<bool>
  ModelTpl<Scalar, Options, JointCollectionTpl>::hasConfigurationLimitInTangent() const
  {
    std::vector<bool> vec;
    for (Index i = 1; i < (Index)(njoints); ++i)
    {
      const std::vector<bool> & cf_limits = joints[i].hasConfigurationLimitInTangent();
      vec.insert(vec.end(), cf_limits.begin(), cf_limits.end());
    }
    return vec;
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  std::vector<JointIndex> ModelTpl<Scalar, Options, JointCollectionTpl>::getChildJoints() const
  {
    std::vector<JointIndex> res;
    for (JointIndex joint_id = 1; joint_id < JointIndex(njoints); ++joint_id)
    {
      if (this->children[joint_id].empty())
        res.push_back(joint_id);
    }
    return res;
  }

} // namespace pinocchio

/// @endcond

#endif // ifndef __pinocchio_multibody_model_hxx__
