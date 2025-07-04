//
// Copyright (c) 2019 CNRS INRIA
//

#ifndef __pinocchio_algorithm_model_hpp__
#define __pinocchio_algorithm_model_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/geometry.hpp"

namespace pinocchio
{
  /**
   *  \brief Append a child model into a parent model, after a specific frame given by its index.
   *
   *  \param[in] modelA the parent model.
   *  \param[in] modelB the child model.
   *  \param[in] frameInModelA index of the frame of modelA where to append modelB.
   *  \param[in] aMb pose of modelB universe joint (index 0) in frameInModelA.
   *  \param[out] model the resulting model.
   *
   *  The order of the joints in the output models are
   *    - joints of modelA up to the parent of FrameInModelA,
   *    - all the descendents of parent of FrameInModelA,
   *    - the remaining joints of modelA.
   */
  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  void appendModel(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & modelA,
    const ModelTpl<Scalar, Options, JointCollectionTpl> & modelB,
    const FrameIndex frameInModelA,
    const SE3Tpl<Scalar, Options> & aMb,
    ModelTpl<Scalar, Options, JointCollectionTpl> & model);

  /**
   *  \brief Append a child model into a parent model, after a specific frame given by its index.
   *
   *  \param[in] modelA the parent model.
   *  \param[in] modelB the child model.
   *  \param[in] frameInModelA index of the frame of modelA where to append modelB.
   *  \param[in] aMb pose of modelB universe joint (index 0) in frameInModelA.
   *
   *  \return A new model containing the fusion of modelA and modelB.
   *
   *  The order of the joints in the output models are
   *    - joints of modelA up to the parent of FrameInModelA,
   *    - all the descendents of parent of FrameInModelA,
   *    - the remaining joints of modelA.
   */
  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  ModelTpl<Scalar, Options, JointCollectionTpl> appendModel(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & modelA,
    const ModelTpl<Scalar, Options, JointCollectionTpl> & modelB,
    const FrameIndex frameInModelA,
    const SE3Tpl<Scalar, Options> & aMb)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    Model model;

    appendModel(modelA, modelB, frameInModelA, aMb, model);

    return model;
  }

  /**
   *  \brief Append a child model into a parent model, after a specific frame given by its index.
   *
   *  \param[in] modelA the parent model.
   *  \param[in] modelB the child model.
   *  \param[in] geomModelA the parent geometry model.
   *  \param[in] geomModelB the child geometry model.
   *  \param[in] frameInModelA index of the frame of modelA where to append modelB.
   *  \param[in] aMb pose of modelB universe joint (index 0) in frameInModelA.
   *  \param[out] model the resulting model.
   *  \param[out] geomModel the resulting geometry model.
   */
  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  void appendModel(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & modelA,
    const ModelTpl<Scalar, Options, JointCollectionTpl> & modelB,
    const GeometryModel & geomModelA,
    const GeometryModel & geomModelB,
    const FrameIndex frameInModelA,
    const SE3Tpl<Scalar, Options> & aMb,
    ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    GeometryModel & geomModel);

  /**
   *
   *  \brief Build a reduced model from a given input model and a list of joint to lock.
   *
   *  \param[in] model the input model to reduce.
   *  \param[in] list_of_joints_to_lock list of joints to lock in the input model.
   *  \param[in] reference_configuration reference configuration.
   *  \param[out] reduced_model the reduced model.
   *
   *  \remarks All the joints that have been set to be fixed in the new reduced_model now appear in
   * the kinematic tree as a Frame as FIXED_JOINT.
   *
   *  \todo At the moment, the joint and geometry order is kept while the frames
   *  are re-ordered in a hard to predict way. Their order could be kept.
   *
   */
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType>
  void buildReducedModel(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    std::vector<JointIndex> list_of_joints_to_lock,
    const Eigen::MatrixBase<ConfigVectorType> & reference_configuration,
    ModelTpl<Scalar, Options, JointCollectionTpl> & reduced_model);

  /**
   *
   *  \brief Build a reduced model from a given input model and a list of joint to lock.
   *
   *  \param[in] model the input model to reduce.
   *  \param[in] list_of_joints_to_lock list of joints to lock in the input model.
   *  \param[in] reference_configuration reference configuration.
   *
   *  \returns A reduce model of the input model.
   *
   *  \remarks All the joints that have been set to be fixed in the new reduced_model now appear in
   * the kinematic tree as a Frame as FIXED_JOINT.
   *
   */
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType>
  ModelTpl<Scalar, Options, JointCollectionTpl> buildReducedModel(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const std::vector<JointIndex> & list_of_joints_to_lock,
    const Eigen::MatrixBase<ConfigVectorType> & reference_configuration)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    Model reduced_model;

    buildReducedModel(model, list_of_joints_to_lock, reference_configuration, reduced_model);

    return reduced_model;
  }

  /**
   *
   *  \brief Build a reduced model and a rededuced geometry model  from a given input model, a given
   * input geometry model and a list of joint to lock.
   *
   *  \param[in] model the input model to reduce.
   *  \param[in] geom_model the input geometry model to reduce.
   *  \param[in] list_of_joints_to_lock list of joints to lock in the input model.
   *  \param[in] reference_configuration reference configuration.
   *  \param[out] reduced_model the reduced model.
   *  \param[out] reduced_geom_model the reduced geometry model.
   *
   *  \remarks All the joints that have been set to be fixed in the new reduced_model now appear in
   * the kinematic tree as a Frame as FIXED_JOINT.
   *
   */
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType>
  void buildReducedModel(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const GeometryModel & geom_model,
    const std::vector<JointIndex> & list_of_joints_to_lock,
    const Eigen::MatrixBase<ConfigVectorType> & reference_configuration,
    ModelTpl<Scalar, Options, JointCollectionTpl> & reduced_model,
    GeometryModel & reduced_geom_model);

  /**
   *
   *  \brief Build a reduced model and a rededuced geometry model  from a given input model, a given
   * input geometry model and a list of joint to lock.
   *
   *  \param[in] model the input model to reduce.
   *  \param[in] list_of_geom_models the input geometry model to reduce (example: visual_model,
   * collision_model). \param[in] list_of_joints_to_lock list of joints to lock in the input model.
   *  \param[in] reference_configuration reference configuration.
   *  \param[out] reduced_model the reduced model.
   *  \param[out] list_of_reduced_geom_model the list of reduced geometry models.
   *
   *  \remarks All the joints that have been set to be fixed in the new reduced_model now appear in
   * the kinematic tree as a Frame as FIXED_JOINT.
   *
   */
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename GeometryModelAllocator,
    typename ConfigVectorType>
  void buildReducedModel(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const std::vector<GeometryModel, GeometryModelAllocator> & list_of_geom_models,
    const std::vector<JointIndex> & list_of_joints_to_lock,
    const Eigen::MatrixBase<ConfigVectorType> & reference_configuration,
    ModelTpl<Scalar, Options, JointCollectionTpl> & reduced_model,
    std::vector<GeometryModel, GeometryModelAllocator> & list_of_reduced_geom_models);

  /**
   *
   *  \brief Transform of a joint of the model into a mimic joint. Keep the type of the joint as it
   * was previously.
   *
   *  \param[in] model the input model to take joints from.
   *  \param[in] index_mimicked index of the joint to mimic
   *  \param[in] index_mimicking index of the joint that will mimic
   *  \param[in] scaling Scaling of joint velocity and configuration
   *  \param[in] offset Offset of joint configuration
   *  \param[out] output_model Model with the joint mimic
   *
   */
  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  void transformJointIntoMimic(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & input_model,
    const JointIndex & index_mimicked,
    const JointIndex & index_mimicking,
    const Scalar & scaling,
    const Scalar & offset,
    ModelTpl<Scalar, Options, JointCollectionTpl> & output_model);

  /**
   *
   *  \brief Transform joints of a model into mimic joints
   *
   *  \param[in] model the input model to take joints from.
   *  \param[in] index_mimicked indexes of the joint to mimic
   *  \param[in] index_mimicking indexes of the joint that will mimic
   *  \param[in] scaling Scalings of joint velocity and configuration
   *  \param[in] offset Offsets of joint configuration
   *  \param[out] output_model Model with the joint mimic
   *
   */
  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  void buildMimicModel(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & input_model,
    const std::vector<JointIndex> & index_mimicked,
    const std::vector<JointIndex> & index_mimicking,
    const std::vector<Scalar> & scaling,
    const std::vector<Scalar> & offset,
    ModelTpl<Scalar, Options, JointCollectionTpl> & output_model);

  /**
   *
   *  \brief Computes the common ancestor between two joints belonging to the same kinematic tree.
   *
   *  \param[in] model the input model.
   *  \param[in] joint1_id index of the first joint.
   *  \param[in] joint2_id index of the second joint.
   *  \param[out] index_ancestor_in_support1 index of the ancestor within model.supports[joint1_id].
   *  \param[out] index_ancestor_in_support2 index of the ancestor within model.supports[joint2_id].
   *
   */
  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  JointIndex findCommonAncestor(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    JointIndex joint1_id,
    JointIndex joint2_id,
    size_t & index_ancestor_in_support1,
    size_t & index_ancestor_in_support2);

  /**
   *
   *  \brief Computes the common ancestor between two joints belonging to the same kinematic tree.
   *
   *  \param[in] model the input model.
   *  \param[in] joint1_id index of the first joint.
   *  \param[in] joint2_id index of the second joint.
   *
   */
  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  JointIndex findCommonAncestor(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    JointIndex joint1_id,
    JointIndex joint2_id)
  {
    size_t index_ancestor_in_support1, index_ancestor_in_support2;
    return findCommonAncestor(
      model, joint1_id, joint2_id, index_ancestor_in_support1, index_ancestor_in_support2);
  }

} // namespace pinocchio

#include "pinocchio/algorithm/model.hxx"

#if PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
  #include "pinocchio/algorithm/model.txx"
#endif // PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION

#endif // ifndef __pinocchio_algorithm_model_hpp__
