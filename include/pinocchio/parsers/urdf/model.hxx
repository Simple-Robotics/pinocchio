//
// Copyright (c) 2015-2021 CNRS INRIA
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_multibody_parsers_urdf_model_hxx__
#define __pinocchio_multibody_parsers_urdf_model_hxx__

#include "pinocchio/math/matrix.hpp"
#include "pinocchio/parsers/config.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <sstream>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <limits>
#include <iostream>

namespace pinocchio
{
  namespace urdf
  {
    namespace details
    {
      typedef double urdf_scalar_type;

      class UrdfVisitor
      {
      public:
        enum JointType
        {
          REVOLUTE,
          CONTINUOUS,
          PRISMATIC,
          FLOATING,
          PLANAR,
          SPHERICAL
        };

        typedef UrdfVisitor Self;

        typedef urdf_scalar_type Scalar;
        typedef SE3Tpl<Scalar, 0> SE3;
        typedef InertiaTpl<Scalar, 0> Inertia;

        typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef Eigen::Ref<Vector> VectorRef;
        typedef Eigen::Ref<const Vector> VectorConstRef;

        typedef ::pinocchio::parsers::Model Model;
        typedef typename Model::JointCollection JointCollection;
        typedef typename Model::JointModel JointModel;
        typedef typename Model::Frame Frame;

        Model & model;
        std::ostream * log;

        UrdfVisitor(Model & model)
        : model(model)
        , log(nullptr)
        {
        }

        template<typename T>
        UrdfVisitor & operator<<(const T & t)
        {
          if (log != nullptr)
            *log << t;
          return *this;
        }

        void setName(const std::string & name)
        {
          model.name = name;
        }

        void addRootJoint(
          const Inertia & Y,
          const std::string & body_name,
          const boost::optional<const JointModel &> root_joint,
          const boost::optional<const std::string &> root_joint_name)
        {
          const Frame & parent_frame = model.frames[0];
          if (root_joint.has_value())
          {
            PINOCCHIO_THROW_IF(
              !root_joint_name.has_value(), std::invalid_argument,
              "if root_joint is provided, root_joint_name must be also be provided.");
            PINOCCHIO_THROW_IF(
              model.existJointName(root_joint_name.get()), std::invalid_argument,
              "root_joint already exists as a joint in the kinematic tree.");

            const Frame & frame = model.frames[0];

            JointIndex idx = model.addJoint(
              frame.parentJoint, root_joint.get(), SE3::Identity(), root_joint_name.get()
              // TODO ,max_effort,max_velocity,min_config,max_config
            );

            FrameIndex jointFrameId = model.addJointFrame(idx, 0);
            appendBodyToJoint(jointFrameId, Y, SE3::Identity(), body_name);
            return;
          }

          // If a root joint has not been provided, we simply add a frame that represents
          // the body (which inertia and name are given as input to this method)
          model.addFrame(
            Frame(body_name, parent_frame.parentJoint, 0, parent_frame.placement, BODY, Y));
          return;
        }

        void addJointAndBody(
          JointType type,
          const Vector3 & axis,
          const FrameIndex & parentFrameId,
          const SE3 & placement,
          const std::string & joint_name,
          const Inertia & Y,
          const std::string & body_name,
          const VectorConstRef & max_effort,
          const VectorConstRef & max_velocity,
          const VectorConstRef & min_config,
          const VectorConstRef & max_config,
          const VectorConstRef & friction,
          const VectorConstRef & damping)
        {
          addJointAndBody(
            type, axis, parentFrameId, placement, joint_name, Y, SE3::Identity(), body_name,
            max_effort, max_velocity, min_config, max_config, friction, damping);
        }

        void addJointAndBody(
          JointType type,
          const Vector3 & axis,
          const FrameIndex & parentFrameId,
          const SE3 & placement,
          const std::string & joint_name,
          const Inertia & Y,
          const SE3 & frame_placement,
          const std::string & body_name,
          const VectorConstRef & max_effort,
          const VectorConstRef & max_velocity,
          const VectorConstRef & min_config,
          const VectorConstRef & max_config,
          const VectorConstRef & friction,
          const VectorConstRef & damping)
        {
          addJointAndBody(
            type, axis, parentFrameId, placement, joint_name, Y, frame_placement, body_name,
            -max_effort, max_effort, -max_velocity, max_velocity, min_config, max_config, -friction,
            friction, damping);
        }
        void addJointAndBody(
          JointType type,
          const Vector3 & axis,
          const FrameIndex & parentFrameId,
          const SE3 & placement,
          const std::string & joint_name,
          const Inertia & Y,
          const SE3 & frame_placement,
          const std::string & body_name,
          const VectorConstRef & min_effort,
          const VectorConstRef & max_effort,
          const VectorConstRef & min_velocity,
          const VectorConstRef & max_velocity,
          const VectorConstRef & min_config,
          const VectorConstRef & max_config,
          const VectorConstRef & min_dry_friction,
          const VectorConstRef & max_dry_friction,
          const VectorConstRef & damping)
        {
          const Vector config_limit_margin =
            Vector::Constant(max_config.size(), static_cast<Scalar>(0));
          addJointAndBody(
            type, axis, parentFrameId, placement, joint_name, Y, frame_placement, body_name,
            min_effort, max_effort, min_velocity, max_velocity, min_config, max_config,
            config_limit_margin, min_dry_friction, max_dry_friction, damping);
        }

        void addJointAndBody(
          JointType type,
          const Vector3 & axis,
          const FrameIndex & parentFrameId,
          const SE3 & placement,
          const std::string & joint_name,
          const Inertia & Y,
          const SE3 & frame_placement,
          const std::string & body_name,
          const VectorConstRef & min_effort,
          const VectorConstRef & max_effort,
          const VectorConstRef & min_velocity,
          const VectorConstRef & max_velocity,
          const VectorConstRef & min_config,
          const VectorConstRef & max_config,
          const VectorConstRef & config_limit_margin,
          const VectorConstRef & min_dry_friction,
          const VectorConstRef & max_dry_friction,
          const VectorConstRef & damping)
        {
          JointIndex joint_id;
          const Frame & frame = model.frames[parentFrameId];
          switch (type)
          {
          case Self::FLOATING:
            joint_id = model.addJoint(
              frame.parentJoint, typename JointCollection::JointModelFreeFlyer(),
              frame.placement * placement, joint_name, min_effort, max_effort, min_velocity,
              max_velocity, min_config, max_config, config_limit_margin, min_dry_friction,
              max_dry_friction, damping);
            break;
          case Self::REVOLUTE:
            joint_id = addJoint<
              typename JointCollection::JointModelRX, typename JointCollection::JointModelRY,
              typename JointCollection::JointModelRZ,
              typename JointCollection::JointModelRevoluteUnaligned>(
              axis, frame, placement, joint_name, min_effort, max_effort, min_velocity,
              max_velocity, min_config, max_config, config_limit_margin, min_dry_friction,
              max_dry_friction, damping);
            break;
          case Self::CONTINUOUS:
            joint_id = addJoint<
              typename JointCollection::JointModelRUBX, typename JointCollection::JointModelRUBY,
              typename JointCollection::JointModelRUBZ,
              typename JointCollection::JointModelRevoluteUnboundedUnaligned>(
              axis, frame, placement, joint_name, min_effort, max_effort, min_velocity,
              max_velocity, min_config, max_config, config_limit_margin, min_dry_friction,
              max_dry_friction, damping);
            break;
          case Self::PRISMATIC:
            joint_id = addJoint<
              typename JointCollection::JointModelPX, typename JointCollection::JointModelPY,
              typename JointCollection::JointModelPZ,
              typename JointCollection::JointModelPrismaticUnaligned>(
              axis, frame, placement, joint_name, min_effort, max_effort, min_velocity,
              max_velocity, min_config, max_config, config_limit_margin, min_dry_friction,
              max_dry_friction, damping);
            break;
          case Self::PLANAR:
            joint_id = model.addJoint(
              frame.parentJoint, typename JointCollection::JointModelPlanar(),
              frame.placement * placement, joint_name, min_effort, max_effort, min_velocity,
              max_velocity, min_config, max_config, config_limit_margin, min_dry_friction,
              max_dry_friction, damping);
            break;
          case Self::SPHERICAL:
            joint_id = model.addJoint(
              frame.parentJoint, typename JointCollection::JointModelSpherical(),
              frame.placement * placement, joint_name, min_effort, max_effort, min_velocity,
              max_velocity, min_config, max_config, config_limit_margin, min_dry_friction,
              max_dry_friction, damping);
            break;
          default:
            PINOCCHIO_CHECK_INPUT_ARGUMENT(false, "The joint type is not correct.");
          };

          FrameIndex jointFrameId = model.addJointFrame(joint_id, (int)parentFrameId);
          appendBodyToJoint(jointFrameId, Y, frame_placement, body_name);
        }

        void addFixedJointAndBody(
          const FrameIndex & parent_frame_id,
          const SE3 & joint_placement,
          const std::string & joint_name,
          const Inertia & Y,
          const std::string & body_name)
        {
          const Frame & parent_frame = model.frames[parent_frame_id];
          const JointIndex parent_frame_parent = parent_frame.parentJoint;

          const SE3 placement = parent_frame.placement * joint_placement;
          FrameIndex fid = model.addFrame(Frame(
            joint_name, parent_frame.parentJoint, parent_frame_id, placement, FIXED_JOINT, Y));

          model.addBodyFrame(body_name, parent_frame_parent, placement, (int)fid);
        }

        void appendBodyToJoint(
          const FrameIndex fid,
          const Inertia & Y,
          const SE3 & placement,
          const std::string & body_name)
        {
          const Frame & frame = model.frames[fid];
          const SE3 & p = frame.placement * placement;
          assert(frame.parentJoint >= 0);
          if (!Y.isZero(Scalar(0)))
          {
            model.appendBodyToJoint(frame.parentJoint, Y, p);
          }

          model.addBodyFrame(body_name, frame.parentJoint, p, (int)fid);
          // Reference to model.frames[fid] can has changed because the vector
          // may have been reallocated.
          assert(model.frames[fid].parentJoint >= 0);
          {
            assert(
              !hasNaN(model.inertias[model.frames[fid].parentJoint].lever())
              && !hasNaN(model.inertias[model.frames[fid].parentJoint].inertia().data()));
          }
        }

        FrameIndex getBodyId(const std::string & frame_name) const
        {

          if (model.existFrame(frame_name, BODY))
          {
            FrameIndex fid = model.getFrameId(frame_name, BODY);
            assert(model.frames[fid].type == BODY);
            return fid;
          }
          else
            throw std::invalid_argument("Model does not have any body named " + frame_name);
        }

        FrameIndex getJointId(const std::string & joint_name) const
        {

          if (model.existJointName(joint_name))
          {
            JointIndex jid = model.getJointId(joint_name);
            return jid;
          }
          else
            throw std::invalid_argument("Model does not have any joint named " + joint_name);
        }

        const std::string & getJointName(const JointIndex jointId) const
        {
          return model.names[jointId];
        }

        Frame getBodyFrame(const std::string & frame_name) const
        {

          if (model.existFrame(frame_name, BODY))
          {
            FrameIndex fid = model.getFrameId(frame_name, BODY);
            assert(model.frames[fid].type == BODY);
            return model.frames[fid];
          }
          else
            throw std::invalid_argument("Model does not have any body named " + frame_name);
        }

        JointIndex getParentId(const std::string & frame_name) const
        {

          if (model.existFrame(frame_name, BODY))
          {
            FrameIndex fid = model.getFrameId(frame_name, BODY);
            assert(model.frames[fid].type == BODY);
            return model.frames[fid].parentJoint;
          }
          else
            throw std::invalid_argument("Model does not have any body named " + frame_name);
        }

        bool existFrame(const std::string & frame_name, const FrameType type) const
        {
          return model.existFrame(frame_name, type);
        }

        template<typename TypeX, typename TypeY, typename TypeZ, typename TypeUnaligned>
        JointIndex addJoint(
          const Vector3 & axis,
          const Frame & frame,
          const SE3 & placement,
          const std::string & joint_name,
          const VectorConstRef & max_effort,
          const VectorConstRef & max_velocity,
          const VectorConstRef & min_config,
          const VectorConstRef & max_config,
          const VectorConstRef & friction,
          const VectorConstRef & damping)
        {
          return addJoint<TypeX, TypeY, TypeZ, TypeUnaligned>(
            axis, frame, placement, joint_name, -max_effort, max_effort, -max_velocity,
            max_velocity, min_config, max_config, -friction, friction, damping);
        }

        template<typename TypeX, typename TypeY, typename TypeZ, typename TypeUnaligned>
        JointIndex addJoint(
          const Vector3 & axis,
          const Frame & frame,
          const SE3 & placement,
          const std::string & joint_name,
          const VectorConstRef & min_effort,
          const VectorConstRef & max_effort,
          const VectorConstRef & min_velocity,
          const VectorConstRef & max_velocity,
          const VectorConstRef & min_config,
          const VectorConstRef & max_config,
          const VectorConstRef & min_dry_friction,
          const VectorConstRef & max_dry_friction,
          const VectorConstRef & damping)
        {
          const Vector config_limit_margin =
            Vector::Constant(max_config.size(), static_cast<Scalar>(0));
          return addJoint<TypeX, TypeY, TypeZ, TypeUnaligned>(
            axis, frame, placement, joint_name, min_effort, max_effort, min_velocity, max_velocity,
            min_config, max_config, config_limit_margin, min_dry_friction, max_dry_friction,
            damping);
        }

        template<typename TypeX, typename TypeY, typename TypeZ, typename TypeUnaligned>
        JointIndex addJoint(
          const Vector3 & axis,
          const Frame & frame,
          const SE3 & placement,
          const std::string & joint_name,
          const VectorConstRef & min_effort,
          const VectorConstRef & max_effort,
          const VectorConstRef & min_velocity,
          const VectorConstRef & max_velocity,
          const VectorConstRef & min_config,
          const VectorConstRef & max_config,
          const VectorConstRef & config_limit_margin,
          const VectorConstRef & min_dry_friction,
          const VectorConstRef & max_dry_friction,
          const VectorConstRef & damping)
        {
          CartesianAxis axisType = extractCartesianAxis(axis);
          switch (axisType)
          {
          case AXIS_X:
            return model.addJoint(
              frame.parentJoint, TypeX(), frame.placement * placement, joint_name, min_effort,
              max_effort, min_velocity, max_velocity, min_config, max_config, config_limit_margin,
              min_dry_friction, max_dry_friction, damping);
            break;

          case AXIS_Y:
            return model.addJoint(
              frame.parentJoint, TypeY(), frame.placement * placement, joint_name, min_effort,
              max_effort, min_velocity, max_velocity, min_config, max_config, config_limit_margin,
              min_dry_friction, max_dry_friction, damping);
            break;

          case AXIS_Z:
            return model.addJoint(
              frame.parentJoint, TypeZ(), frame.placement * placement, joint_name, min_effort,
              max_effort, min_velocity, max_velocity, min_config, max_config, config_limit_margin,
              min_dry_friction, max_dry_friction, damping);
            break;

          case AXIS_UNALIGNED:
            return model.addJoint(
              frame.parentJoint, TypeUnaligned(axis.normalized()), frame.placement * placement,
              joint_name, min_effort, max_effort, min_velocity, max_velocity, min_config,
              max_config, config_limit_margin, min_dry_friction, max_dry_friction, damping);
            break;
          default:
            PINOCCHIO_CHECK_INPUT_ARGUMENT(false, "The axis type of the joint is of wrong type.");
            break;
          }
        }

      private:
        ///
        /// \brief The four possible cartesian types of an 3D axis.
        ///
        enum CartesianAxis
        {
          AXIS_X = 0,
          AXIS_Y = 1,
          AXIS_Z = 2,
          AXIS_UNALIGNED
        };

        ///
        /// \brief Extract the cartesian property of a particular 3D axis.
        ///
        /// \param[in] axis The input URDF axis.
        ///
        /// \return The property of the particular axis pinocchio::urdf::CartesianAxis.
        ///
        static inline CartesianAxis extractCartesianAxis(const Vector3 & axis)
        {
          if (axis.isApprox(Vector3::UnitX()))
            return AXIS_X;
          else if (axis.isApprox(Vector3::UnitY()))
            return AXIS_Y;
          else if (axis.isApprox(Vector3::UnitZ()))
            return AXIS_Z;
          else
            return AXIS_UNALIGNED;
        }
      };

      PINOCCHIO_PARSERS_DLLAPI void parseRootTree(
        const ::urdf::ModelInterface * urdfTree,
        UrdfVisitor & model,
        const boost::optional<const ::pinocchio::parsers::JointModel &> rootJoint = boost::none,
        const boost::optional<const std::string &> rootJointName = boost::none);

      PINOCCHIO_PARSERS_DLLAPI void parseRootTree(
        const std::string & filename,
        UrdfVisitor & model,
        const boost::optional<const ::pinocchio::parsers::JointModel &> rootJoint = boost::none,
        const boost::optional<const std::string &> rootJointName = boost::none);

      PINOCCHIO_PARSERS_DLLAPI void parseRootTreeFromXML(
        const std::string & xmlString,
        UrdfVisitor & model,
        const boost::optional<const ::pinocchio::parsers::JointModel &> rootJoint = boost::none,
        const boost::optional<const std::string &> rootJointName = boost::none);
    } // namespace details

    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    ModelTpl<Scalar, Options, JointCollectionTpl> & buildModel(
      const std::string & filename,
      const typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointModel & rootJoint,
      const std::string & rootJointName,
      ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const bool verbose)
    {
      typedef ::pinocchio::parsers::Model Model;
      typedef ::pinocchio::parsers::JointModel JointModel;
      if (rootJointName.empty())
        throw std::invalid_argument(
          "rootJoint was given without a name. Please fill the argument root_joint_name");

      // copy in case incoming model is not empty
      Model urdf_model = model;
      JointModel root_joint = rootJoint;

      details::UrdfVisitor visitor(urdf_model);
      if (verbose)
        visitor.log = &std::cout;

      boost::optional<const JointModel &> root_joint_opt(root_joint);
      boost::optional<const std::string &> root_joint_name_opt(rootJointName);
      details::parseRootTree(filename, visitor, root_joint_opt, root_joint_name_opt);

      // cast back to the input model's joint collection
      model = visitor.model;
      return model;
    }

    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    ModelTpl<Scalar, Options, JointCollectionTpl> & buildModel(
      const std::string & filename,
      const typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointModel & rootJoint,
      ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const bool verbose)
    {
      return buildModel(filename, rootJoint, "root_joint", model, verbose);
    }

    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    ModelTpl<Scalar, Options, JointCollectionTpl> & buildModel(
      const std::string & filename,
      ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const bool verbose)
    {
      typedef ::pinocchio::parsers::Model Model;
      Model urdf_model = model;

      details::UrdfVisitor visitor(urdf_model);
      if (verbose)
        visitor.log = &std::cout;
      details::parseRootTree(filename, visitor);

      model = visitor.model;
      return model;
    }

    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    ModelTpl<Scalar, Options, JointCollectionTpl> & buildModelFromXML(
      const std::string & xmlStream,
      const typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointModel & rootJoint,
      const std::string & rootJointName,
      ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const bool verbose)
    {
      typedef ::pinocchio::parsers::Model Model;
      typedef ::pinocchio::parsers::JointModel JointModel;
      if (rootJointName.empty())
        throw std::invalid_argument(
          "rootJoint was given without a name. Please fill the argument rootJointName");

      Model urdf_model = model;
      JointModel root_joint = rootJoint;

      details::UrdfVisitor visitor(urdf_model);
      if (verbose)
        visitor.log = &std::cout;

      boost::optional<const JointModel &> root_joint_opt(root_joint);
      boost::optional<const std::string &> root_joint_name_opt(rootJointName);
      details::parseRootTreeFromXML(xmlStream, visitor, root_joint_opt, root_joint_name_opt);

      model = visitor.model;
      return model;
    }

    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    ModelTpl<Scalar, Options, JointCollectionTpl> & buildModelFromXML(
      const std::string & xmlStream,
      const typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointModel & rootJoint,
      ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const bool verbose)
    {
      return buildModelFromXML(xmlStream, rootJoint, "root_joint", model, verbose);
    }

    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    ModelTpl<Scalar, Options, JointCollectionTpl> & buildModelFromXML(
      const std::string & xmlStream,
      ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const bool verbose)
    {
      typedef ::pinocchio::parsers::Model Model;
      Model urdf_model = model;

      details::UrdfVisitor visitor(urdf_model);
      if (verbose)
        visitor.log = &std::cout;
      details::parseRootTreeFromXML(xmlStream, visitor);

      model = visitor.model;
      return model;
    }

    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    ModelTpl<Scalar, Options, JointCollectionTpl> & buildModel(
      const std::shared_ptr<::urdf::ModelInterface> urdfTree,
      const typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointModel & rootJoint,
      const std::string & rootJointName,
      ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const bool verbose)
    {
      typedef ::pinocchio::parsers::Model Model;
      typedef ::pinocchio::parsers::JointModel JointModel;
      if (rootJointName.empty())
        throw std::invalid_argument(
          "rootJoint was given without a name. Please fill the argument rootJointName");

      PINOCCHIO_CHECK_INPUT_ARGUMENT(urdfTree != NULL);

      Model urdf_model = model;
      JointModel root_joint = rootJoint;

      details::UrdfVisitor visitor(urdf_model);
      if (verbose)
        visitor.log = &std::cout;

      boost::optional<const JointModel &> root_joint_opt(root_joint);
      boost::optional<const std::string &> root_joint_name_opt(rootJointName);
      details::parseRootTree(urdfTree.get(), visitor, root_joint_opt, root_joint_name_opt);

      model = visitor.model;
      return model;
    }

    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    ModelTpl<Scalar, Options, JointCollectionTpl> & buildModel(
      const std::shared_ptr<::urdf::ModelInterface> urdfTree,
      const typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointModel & rootJoint,
      ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const bool verbose)
    {
      return buildModel(urdfTree, rootJoint, "root_joint", model, verbose);
    }

    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    ModelTpl<Scalar, Options, JointCollectionTpl> & buildModel(
      const std::shared_ptr<::urdf::ModelInterface> urdfTree,
      ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const bool verbose)
    {
      typedef ::pinocchio::parsers::Model Model;
      PINOCCHIO_CHECK_INPUT_ARGUMENT(urdfTree != NULL);
      Model urdf_model = model;

      details::UrdfVisitor visitor(urdf_model);
      if (verbose)
        visitor.log = &std::cout;
      details::parseRootTree(urdfTree.get(), visitor);

      model = visitor.model;
      return model;
    }

  } // namespace urdf
} // namespace pinocchio

#endif // ifndef __pinocchio_multibody_parsers_urdf_model_hxx__
