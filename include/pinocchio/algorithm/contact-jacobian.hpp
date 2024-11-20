//
// Copyright (c) 2021-2024 INRIA
//

#ifndef __pinocchio_algorithm_contact_jacobian_hpp__
#define __pinocchio_algorithm_contact_jacobian_hpp__

#include "pinocchio/algorithm/constraints/constraints.hpp"

namespace pinocchio
{

  ///
  /// \brief Evaluates all the constraints by calling cmodel.calc().
  ///
  /// \remarks This function assumes that data is up-to-date.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint datas.
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  void evalConstraints(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas);

  ///
  /// \brief Maps the constraint forces expressed in the constraint space to joint forces expressed
  /// in the local frame.
  ///
  /// \remarks This function assumes that the constrained data are up-to-date.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint datas.
  /// \param[in] constraint_forces Matrix or vector containing the constraint forces.
  /// \param[out] joint_forces Vector of  joint forces (dimension model.njoints).
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename ForceMatrix,
    class ForceAllocator>
  void mapConstraintForcesToJointForces(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<ForceMatrix> & constraint_forces,
    std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces);

  ///
  /// \brief Maps the joint motions expressed in the joint space local frame to the constraint
  /// motions.
  ///
  /// \remarks This function assumes that the constrained data are up-to-date.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint datas.
  /// \param[in] joint_motions Vector of  joint motions (dimension model.njoints).
  /// \param[out] constraint_motions Resulting matrix or vector containing the constraint motions.
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    class MotionAllocator,
    typename MotionMatrix>
  void mapJointMotionsToConstraintMotions(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
    const Eigen::MatrixBase<MotionMatrix> & constraint_motions);

  ///
  /// \brief Computes the kinematic Jacobian associatied to a given constraint model.
  ///
  /// \remarks This function assumes that the a computeJointJacobians has been called first or any
  /// algorithms that computes data.J and data.oMi.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_model Constraint model.
  /// \param[in] constraint_data Constraint data.
  /// \param[out] J A reference on the Jacobian matrix where the results will be stored in (dim 6 x
  /// model.nv). You must fill J with zero elements, e.g. J.fill(0.).
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConstraintModelDerived,
    typename ConstraintDataDerived,
    typename Matrix6Like>
  PINOCCHIO_UNSUPPORTED_MESSAGE("The API will change towards more flexibility")
  void getConstraintJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintModelBase<ConstraintModelDerived> & constraint_model,
    ConstraintDataBase<ConstraintDataDerived> & constraint_data,
    const Eigen::MatrixBase<Matrix6Like> & J);

  ///
  /// \brief Computes the kinematic Jacobian associatied to a given set of constraint models.
  ///
  /// \remarks This function assumes that the a computeJointJacobians has been called first or any
  /// algorithms that computes data.J and data.oMi.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint data.
  /// \param[out] J A reference on the Jacobian matrix where the results will be stored in (dim nc x
  /// model.nv). You must fill J with zero elements, e.g. J.fill(0.).
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    template<typename T> class Holder,
    typename DynamicMatrixLike,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  void getConstraintsJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & constraint_model,
    std::vector<Holder<ConstraintData>, ConstraintDataAllocator> & constraint_data,
    const Eigen::MatrixBase<DynamicMatrixLike> & J);

  ///
  /// \brief Computes the kinematic Jacobian associatied to a given set of constraint models.
  ///
  /// \remarks This function assumes that the a computeJointJacobians has been called first or any
  /// algorithms that computes data.J and data.oMi.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint data.
  /// \param[out] J A reference on the Jacobian matrix where the results will be stored in (dim nc x
  /// model.nv). You must fill J with zero elements, e.g. J.fill(0.).
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename DynamicMatrixLike,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  void getConstraintsJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintDataAllocator> & constraint_model,
    std::vector<ConstraintData, ConstraintDataAllocator> & constraint_data,
    const Eigen::MatrixBase<DynamicMatrixLike> & J);

  ///
  /// \brief Evaluate the operation res = J.T * rhs
  ///
  /// \remarks This function assumes that the a computeJointJacobians has been called first or any
  /// algorithms that computes data.J and data.oMi.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint data.
  /// \param[in] rhs Right-hand side term.
  /// \param[out] res Results.
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename RhsMatrixType,
    typename ResultMatrixType>
  void evalConstraintJacobianTransposeProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<RhsMatrixType> & rhs,
    const Eigen::MatrixBase<ResultMatrixType> & res);

} // namespace pinocchio

#include "pinocchio/algorithm/contact-jacobian.hxx"

#if PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
  #include "pinocchio/algorithm/contact-jacobian.txx"
#endif // PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION

#endif // ifndef __pinocchio_algorithm_contact_jacobian_hpp__
