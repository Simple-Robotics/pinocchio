//
// Copyright (c) 2015-2024 CNRS INRIA
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_multibody_model_hpp__
#define __pinocchio_multibody_model_hpp__

#include "pinocchio/spatial/fwd.hpp"
#include "pinocchio/spatial/se3.hpp"
#include "pinocchio/spatial/force.hpp"
#include "pinocchio/spatial/motion.hpp"
#include "pinocchio/spatial/inertia.hpp"

#include "pinocchio/common/model-entity.hpp"

#include "pinocchio/multibody/fwd.hpp"
#include "pinocchio/multibody/frame.hpp"
#include "pinocchio/multibody/joint/joint-generic.hpp"

#include "pinocchio/container/aligned-vector.hpp"

#include "pinocchio/serialization/serializable.hpp"

#include <map>
#include <iterator>

namespace pinocchio
{
  template<
    typename NewScalar,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl>
  struct CastType<NewScalar, ModelTpl<Scalar, Options, JointCollectionTpl>>
  {
    typedef ModelTpl<NewScalar, Options, JointCollectionTpl> type;
  };

  template<typename _Scalar, int _Options, template<typename, int> class JointCollectionTpl>
  struct traits<ModelTpl<_Scalar, _Options, JointCollectionTpl>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
    typedef JointCollectionTpl<Scalar, Options> JointCollection;
  };

  template<typename _Scalar, int _Options, template<typename, int> class JointCollectionTpl>
  struct ModelTpl
  : serialization::Serializable<ModelTpl<_Scalar, _Options, JointCollectionTpl>>
  , NumericalBase<ModelTpl<_Scalar, _Options, JointCollectionTpl>>
  , ModelEntity<ModelTpl<_Scalar, _Options, JointCollectionTpl>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef typename traits<ModelTpl>::Scalar Scalar;
    enum
    {
      Options = traits<ModelTpl>::Options
    };

    typedef typename traits<ModelTpl>::JointCollection JointCollection;
    typedef typename traits<ModelTpl>::Data Data;

    typedef SE3Tpl<Scalar, Options> SE3;
    typedef MotionTpl<Scalar, Options> Motion;
    typedef ForceTpl<Scalar, Options> Force;
    typedef InertiaTpl<Scalar, Options> Inertia;
    typedef FrameTpl<Scalar, Options> Frame;

    typedef pinocchio::Index Index;
    typedef pinocchio::JointIndex JointIndex;
    typedef pinocchio::GeomIndex GeomIndex;
    typedef pinocchio::FrameIndex FrameIndex;
    typedef std::vector<Index> IndexVector;

    typedef JointModelTpl<Scalar, Options, JointCollectionTpl> JointModel;
    typedef JointDataTpl<Scalar, Options, JointCollectionTpl> JointData;

    typedef PINOCCHIO_ALIGNED_STD_VECTOR(JointModel) JointModelVector;
    typedef PINOCCHIO_ALIGNED_STD_VECTOR(JointData) JointDataVector;

    typedef PINOCCHIO_ALIGNED_STD_VECTOR(Frame) FrameVector;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> MatrixXs;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;

    typedef PINOCCHIO_ALIGNED_STD_VECTOR(Inertia) InertiaVector;
    typedef PINOCCHIO_ALIGNED_STD_VECTOR(SE3) SE3Vector;

    /// \brief Dense vectorized version of a joint configuration vector.
    typedef VectorXs ConfigVectorType;

    /// \brief Map between a string (key) and a configuration vector
    typedef std::map<std::string, ConfigVectorType> ConfigVectorMap;

    /// \brief Dense vectorized version of a joint tangent vector (e.g. velocity, acceleration,
    /// etc).
    ///        It also handles the notion of co-tangent vector (e.g. torque, etc).
    typedef VectorXs TangentVectorType;

    /// \brief Dimension of the configuration vector representation.
    int nq;

    /// \brief Dimension of the velocity vector space.
    int nv;

    /// \brief Dimension of the jacobian space.
    int nvExtended;

    /// \brief Number of joints.
    int njoints;

    /// \brief Number of bodies.
    int nbodies;

    /// \brief Number of operational frames.
    int nframes;

    /// \brief Vector of spatial inertias supported by each joint.
    InertiaVector inertias;

    /// \brief Vector of joint placements: placement of a joint *i* wrt its parent joint frame.
    SE3Vector jointPlacements;

    /// \brief Vector of joint models.
    JointModelVector joints;

    /// \brief Vector of starting index of the *i*th  joint in the configuration space.
    std::vector<int> idx_qs;

    /// \brief Vector of dimension of the  joint configuration subspace.
    std::vector<int> nqs;

    /// \brief Starting index of the *i*th joint in the tangent configuration space.
    std::vector<int> idx_vs;

    /// \brief Dimension of the *i*th joint tangent subspace.
    std::vector<int> nvs;

    /// \brief Starting index of the *i*th joint in the jacobian space.
    std::vector<int> idx_vExtendeds;

    /// \brief Dimension of the *i*th joint jacobian subspace.
    std::vector<int> nvExtendeds;

    /// \brief Vector of parent joint indexes. The parent of joint *i*, denoted *li*, corresponds to
    /// li==parents[i].
    std::vector<JointIndex> parents;

    /// \brief Vector of children index. Chidren of the *i*th joint, denoted *mu(i)* corresponds to
    /// the set (i==parents[k] for k in mu(i)).
    std::vector<IndexVector> children;

    /// \brief Vector of mimicking joints in the tree (with type MimicTpl)
    std::vector<JointIndex> mimicking_joints;

    /// \brief Vector of mimicked joints in the tree (can be any joint type)
    /// The i-th element of this vector correspond to the mimicked joint of the i-th mimicking
    /// vector in mimicking_joints
    std::vector<JointIndex> mimicked_joints;

    /// \brief Name of the joints.
    std::vector<std::string> names;

    /// \brief Map of reference configurations, indexed by user given names.
    ConfigVectorMap referenceConfigurations;

    /// \brief Vector of armature values expressed at the joint level
    /// This vector may contain the contribution of rotor inertia effects for instance.
    VectorXs armature;

    /// \brief Vector of rotor inertia parameters
    TangentVectorType rotorInertia;

    /// \brief Vector of rotor gear ratio parameters
    TangentVectorType rotorGearRatio;

    /// \brief Vector of joint friction parameters
    /// Deprecated in favor of lowerDryFrictionLimit and upperDryFrictionLimit
    PINOCCHIO_DEPRECATED TangentVectorType & friction;

    /// \brief Vector of joint friction parameters
    TangentVectorType lowerDryFrictionLimit;

    /// \brief Vector of joint friction parameters
    TangentVectorType upperDryFrictionLimit;

    /// \brief Vector of joint damping parameters
    TangentVectorType damping;

    /// \brief Vector of minimal joint torques
    TangentVectorType lowerEffortLimit;

    /// \brief Vector of maximal joint torques
    TangentVectorType upperEffortLimit;

    /// \brief Vector of maximal joint torques
    /// Deprecated in favor of lowerEffortLimit and upperEffortLimit
    PINOCCHIO_DEPRECATED TangentVectorType & effortLimit;

    /// \brief Vector of minimal joint velocities
    TangentVectorType lowerVelocityLimit;

    /// \brief Vector of maximal joint velocities
    TangentVectorType upperVelocityLimit;

    /// \brief Vector of maximal joint velocities
    /// Deprecated in favor of lowerVelocityLimit and upperVelocityLimit
    PINOCCHIO_DEPRECATED TangentVectorType & velocityLimit;

    /// \brief Lower joint configuration limit
    ConfigVectorType lowerPositionLimit;

    /// \brief Upper joint configuration limit
    ConfigVectorType upperPositionLimit;

    /// \brief Joint configuration limit margin
    ConfigVectorType positionLimitMargin;

    /// \brief Vector of operational frames registered on the model.
    FrameVector frames;

    /// \brief Vector of joint supports.
    /// supports[j] corresponds to the vector of indices of the joints located on the path between
    /// joint *j*  and "universe".
    /// The first element of supports[j] is "universe", the last one is the index of joint *j*
    /// itself.
    std::vector<IndexVector> supports;

    /// \brief Vector of mimic supports joints.
    /// mimic_joint_supports[j] corresponds to the vector of mimic joints indices located on the
    /// path between joint *j*  and "universe". The first element of mimic_joint_supports[j] is
    /// "universe". If *j* is a mimic, the last element is the index of joint *j* itself.
    std::vector<IndexVector> mimic_joint_supports;

    /// \brief Vector of joint subtrees.
    /// subtree[j] corresponds to the subtree supported by the joint *j*.
    /// The first element of subtree[j] is the index of the joint *j* itself.
    std::vector<IndexVector> subtrees;

    /// \brief Spatial gravity of the model.
    Motion gravity;

    /// \brief Default 3D gravity vector (=(0,0,-9.81)).
    static const Vector3 gravity981;

    /// \brief Model name.
    std::string name;

    /// \brief Default constructor. Builds an empty model with no joints.
    ModelTpl()
    : nq(0)
    , nv(0)
    , nvExtended(0)
    , njoints(1)
    , nbodies(1)
    , nframes(0)
    , inertias(1, Inertia::Zero())
    , jointPlacements(1, SE3::Identity())
    , joints(1)
    , idx_qs(1, 0)
    , nqs(1, 0)
    , idx_vs(1, 0)
    , nvs(1, 0)
    , idx_vExtendeds(1, 0)
    , nvExtendeds(1, 0)
    , parents(1, 0)
    , children(1)
    , names(1)
    , friction(upperDryFrictionLimit)
    , effortLimit(upperEffortLimit)
    , velocityLimit(upperVelocityLimit)
    , supports(1, IndexVector(1, 0))
    , mimic_joint_supports(1, IndexVector(1, 0))
    , subtrees(1)
    , gravity(gravity981, Vector3::Zero())
    {
      names[0] = "universe"; // Should be "universe joint (trivial)"
      // FIXME Should the universe joint be a FIXED_JOINT even if it is
      // in the list of joints ? See comment in definition of
      // Model::addJointFrame and Model::addBodyFrame
      addFrame(Frame("universe", 0, 0, SE3::Identity(), FIXED_JOINT));
    }

    ///
    /// \brief Copy constructor by casting
    ///
    /// \param[in] other model to copy to *this
    ///
    template<typename S2, int O2>
    explicit ModelTpl(const ModelTpl<S2, O2> & other)
    : friction(upperDryFrictionLimit)
    , effortLimit(upperEffortLimit)
    , velocityLimit(upperVelocityLimit)
    {
      *this = other.template cast<Scalar>();
    }

    ///
    /// \brief Copy constructor from another collection
    ///
    /// \param[in] other model to copy to *this
    ///
    template<template<typename, int> class OtherJointCollectionTpl>
    ModelTpl(const ModelTpl<Scalar, Options, OtherJointCollectionTpl> & other)
    : friction(upperDryFrictionLimit)
    , effortLimit(upperEffortLimit)
    , velocityLimit(upperVelocityLimit)
    {
      *this = other;
    }

    ///
    /// \brief Copy constructor.
    ///
    /// \param[in] other model to copy to *this
    ///
    ModelTpl(const ModelTpl & other)
    : friction(upperDryFrictionLimit)
    , effortLimit(upperEffortLimit)
    , velocityLimit(upperVelocityLimit)
    {
      *this = other;
    }

    /// \returns A new copy of *this with the Scalar type casted to NewScalar.
    template<typename NewScalar>
    typename CastType<NewScalar, ModelTpl>::type cast() const;

    ///
    /// \brief Equality comparison operator.
    ///
    /// \returns true if *this is equal to other.
    ///
    bool operator==(const ModelTpl & other) const;

    ///
    /// \brief Assignment operator from another collection.
    ///
    ///
    template<template<typename, int> class OtherJointCollectionTpl>
    ModelTpl & operator=(const ModelTpl<Scalar, Options, OtherJointCollectionTpl> & other);

    ///
    /// \brief Assignment operator.
    ///
    ///
    ModelTpl & operator=(const ModelTpl & other)
    {
      (*this).template operator= <JointCollectionTpl>(other);
      return *this;
    }

    ///
    /// \returns true if *this is NOT equal to other.
    ///
    bool operator!=(const ModelTpl & other) const
    {
      return !(*this == other);
    }

    ///
    /// \brief Add a joint to the kinematic tree with infinite bounds.
    ///
    /// \remarks This method does not add a Frame of same name to the vector of frames.
    ///         Use Model::addJointFrame.
    /// \remarks The inertia supported by the joint is set to Zero.
    /// \remark Joints need to be added to the tree in a depth-first order.
    ///
    /// \tparam JointModelDerived The type of the joint model.
    ///
    /// \param[in] parent Index of the parent joint.
    /// \param[in] joint_model The joint model.
    /// \param[in] joint_placement Placement of the joint inside its parent joint.
    /// \param[in] joint_name Name of the joint. If empty, the name is random.
    ///
    /// \return The index of the new joint.
    ///
    /// \sa Model::appendBodyToJoint
    ///
    JointIndex addJoint(
      const JointIndex parent,
      const JointModel & joint_model,
      const SE3 & joint_placement,
      const std::string & joint_name);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &)
    /// Deprecated in favor of the constructor using min and max effort/velocity
    ///
    /// \param[in] max_effort Maximal joint torque.
    /// \param[in] max_velocity Maximal joint velocity.
    /// \param[in] min_config Lower joint configuration.
    /// \param[in] max_config Upper joint configuration.
    ///
    JointIndex addJoint(
      const JointIndex parent,
      const JointModel & joint_model,
      const SE3 & joint_placement,
      const std::string & joint_name,
      const VectorXs & max_effort,
      const VectorXs & max_velocity,
      const VectorXs & min_config,
      const VectorXs & max_config);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &)
    /// Deprecated in favor of the constructor using min and max effort/velocity
    ///
    /// \param[in] max_effort Maximal joint torque.
    /// \param[in] max_velocity Maximal joint velocity.
    /// \param[in] min_config Lower joint configuration.
    /// \param[in] max_config Upper joint configuration.
    /// \param[in] config_limit_margin Joint configuration limit margin.
    ///
    JointIndex addJoint(
      const JointIndex parent,
      const JointModel & joint_model,
      const SE3 & joint_placement,
      const std::string & joint_name,
      const VectorXs & max_effort,
      const VectorXs & max_velocity,
      const VectorXs & min_config,
      const VectorXs & max_config,
      const VectorXs & config_limit_margin);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &,const VectorXs &,const VectorXs &,const VectorXs &,const VectorXs &)
    /// Deprecated in favor of the constructor using min and max effort/velocity
    ///
    /// \param[in] min_effort Minimal joint torque.
    /// \param[in] min_velocity Minimal joint velocity.
    /// \param[in] min_friction Minimal joint friction parameters.
    /// \param[in] max_friction Maximal joint friction parameters.
    /// \param[in] damping Joint damping parameters.
    ///
    JointIndex addJoint(
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
      const VectorXs & min_friction,
      const VectorXs & max_friction,
      const VectorXs & damping);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &,const VectorXs &,const VectorXs &,const VectorXs &,const VectorXs &)
    /// Deprecated in favor of the constructor using min and max effort/velocity
    ///
    /// \param[in] config_limit_margin Joint configuration limit margin.
    /// \param[in] min_effort Minimal joint torque.
    /// \param[in] min_velocity Minimal joint velocity.
    /// \param[in] min_friction Minimal joint friction parameters.
    /// \param[in] max_friction Maximal joint friction parameters.
    /// \param[in] damping Joint damping parameters.
    ///
    JointIndex addJoint(
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
      const VectorXs & min_friction,
      const VectorXs & max_friction,
      const VectorXs & damping);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &,const VectorXs &,const VectorXs &,const VectorXs &,const VectorXs &)
    ///
    /// \param[in] config_limit_margin Joint configuration limit margin.
    /// \param[in] friction Joint friction parameters.
    /// \param[in] damping Joint damping parameters.
    ///
    JointIndex addJoint(
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
      const VectorXs & damping);

    ///
    /// \copydoc ModelTpl::addJoint(const JointIndex,const JointModel &,const SE3 &,const
    /// std::string &,const VectorXs &,const VectorXs &,const VectorXs &,const VectorXs &)
    ///
    /// \param[in] friction Joint friction parameters.
    /// \param[in] damping Joint damping parameters.
    ///
    JointIndex addJoint(
      const JointIndex parent,
      const JointModel & joint_model,
      const SE3 & joint_placement,
      const std::string & joint_name,
      const VectorXs & max_effort,
      const VectorXs & max_velocity,
      const VectorXs & min_config,
      const VectorXs & max_config,
      const VectorXs & friction,
      const VectorXs & damping);

    ///
    /// \brief Add a joint to the frame tree.
    ///
    /// \param[in] jointIndex Index of the joint.
    /// \param[in] frameIndex Index of the parent frame. If negative,
    ///            the parent frame is the frame of the parent joint.
    ///
    /// \return The index of the new frame
    ///
    FrameIndex addJointFrame(const JointIndex joint_index, int previous_frame_index = -1);

    ///
    /// \brief Append a body to a given joint of the kinematic tree.
    ///
    /// \param[in] joint_index Index of the supporting joint.
    /// \param[in] Y Spatial inertia of the body.
    /// \param[in] body_placement The relative placement of the body regarding to the parent joint.
    /// Set default to the Identity placement.
    ///
    /// \sa Model::addJoint
    ///
    void appendBodyToJoint(
      const JointIndex joint_index,
      const Inertia & Y,
      const SE3 & body_placement = SE3::Identity());

    ///
    /// \brief Add a body to the frame tree.
    ///
    /// \param[in] body_name Name of the body.
    /// \param[in] parentJoint Index of the parent joint.
    /// \param[in] body_placement The relative placement of the body regarding to the parent joint.
    /// Set default to the Identity placement. \param[in] parentFrame Index of the parent frame. If
    /// negative,
    ///            the parent frame is the frame of the parent joint.
    ///
    /// \return The index of the new frame
    ///
    FrameIndex addBodyFrame(
      const std::string & body_name,
      const JointIndex & parentJoint,
      const SE3 & body_placement = SE3::Identity(),
      int parentFrame = -1);

    ///
    /// \brief Return the index of a body given by its name.
    ///
    /// \warning If no body is found, return the number of elements at time T.
    /// This can lead to errors if the model is expanded after this method is called
    /// (for example to get the id of a parent body)
    ///
    /// \param[in] name Name of the body.
    ///
    /// \return Index of the body.
    ///
    FrameIndex getBodyId(const std::string & name) const;

    ///
    /// \brief Check if a body given by its name exists.
    ///
    /// \param[in] name Name of the body.
    ///
    /// \return True if the body exists in the kinematic tree.
    ///
    bool existBodyName(const std::string & name) const;

    ///
    /// \brief Return the index of a joint given by its name.
    ///
    /// \warning If no joint is found, return the number of elements at time T.
    /// This can lead to errors if the model is expanded after this method is called
    /// (for example to get the id of a parent joint)
    /// \param[in] name Name of the joint.
    ///
    /// \return Index of the joint.
    ///
    JointIndex getJointId(const std::string & name) const;

    ///
    /// \brief Check if a joint given by its name exists.
    ///
    /// \param[in] name Name of the joint.
    ///
    /// \return True if the joint exists in the kinematic tree.
    ///
    bool existJointName(const std::string & name) const;

    ///
    /// \brief Returns the index of a frame given by its name.
    ///        \sa Model::existFrame to check if the frame exists or not.
    ///
    /// \warning If no frame is found, returns the size of the vector of Model::frames.
    /// This can lead to errors if the model is expanded after this method is called
    /// (for example to get the id of a parent frame).
    ///
    /// \param[in] name Name of the frame.
    /// \param[in] type Type of the frame.
    ///
    /// \return Index of the frame.
    ///
    FrameIndex getFrameId(
      const std::string & name,
      const FrameType & type = (FrameType)(JOINT | FIXED_JOINT | BODY | OP_FRAME | SENSOR)) const;

    ///
    /// \brief Checks if a frame given by its name exists.
    ///
    /// \param[in] name Name of the frame.
    /// \param[in] type Type of the frame.
    ///
    /// \return Returns true if the frame exists.
    ///
    bool existFrame(
      const std::string & name,
      const FrameType & type = (FrameType)(JOINT | FIXED_JOINT | BODY | OP_FRAME | SENSOR)) const;

    ///
    /// \brief Adds a frame to the kinematic tree.
    ///        The inertia stored within the frame will be happened to the inertia supported by the
    ///        joint (frame.parentJoint).
    ///
    /// \param[in] frame The frame to add to the kinematic tree.
    /// \param[in] append_inertia Append the inertia contained in the Frame to the inertia supported
    /// by the joint.
    ///
    /// \return Returns the index of the frame if it has been successfully added or if it already
    /// exists in the kinematic tree.
    ///
    FrameIndex addFrame(const Frame & frame, const bool append_inertia = true);

    ///
    /// \brief Check the validity of the attributes of Model with respect to the specification of
    /// some algorithms.
    ///
    /// The method is a template so that the checkers can be defined in each algorithms.
    /// \param[in] checker a class, typically defined in the algorithm module, that
    /// validates the attributes of model.
    ///
    /// \return true if the Model is valid, false otherwise.
    ///
    template<typename D>
    bool check(const AlgorithmCheckerBase<D> & checker) const
    {
      return checker.checkModel(*this);
    }

    ///
    /// \brief Check if joints have configuration limits
    ///
    /// \return Returns list of boolean of size model.nq.
    ///
    std::vector<bool> hasConfigurationLimit() const;

    ///
    /// \brief Check if joints have configuration limits
    ///
    /// \return Returns list of boolean of size model.nq.
    ///
    std::vector<bool> hasConfigurationLimitInTangent() const;

    /// Run check(fusion::list) with DEFAULT_CHECKERS as argument.
    bool check() const;

    ///
    /// \brief Run checkData on data and current model.
    ///
    /// \param[in] data to be checked wrt *this.
    ///
    /// \return true if the data is valid, false otherwise.
    ///
    bool check(const Data & data) const;

    ///
    /// \brief Create a Data structure associated with the current model
    ///
    Data createData() const;

    /// Returns a vector of the children joints of the kinematic tree.
    /// \remark: a child joint is a node without any child joint.
    std::vector<JointIndex> getChildJoints() const;

  protected:
    ///
    /// \brief Add the joint_id to its parent subtrees.
    ///
    /// \param[in] joint_id The id of the joint to add to the subtrees
    ///
    void addJointIndexToParentSubtrees(const JointIndex joint_id);
  };

} // namespace pinocchio

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
#include "pinocchio/multibody/model.hxx"

#if PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
  #include "pinocchio/multibody/model.txx"
#endif // PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION

#endif // ifndef __pinocchio_multibody_model_hpp__
