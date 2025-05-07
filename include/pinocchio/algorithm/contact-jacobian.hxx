//
// Copyright (c) 2021-2025 INRIA
//

#ifndef __pinocchio_algorithm_contact_jacobian_hxx__
#define __pinocchio_algorithm_contact_jacobian_hxx__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/check.hpp"

namespace pinocchio
{

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
    std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    const size_t num_ee = constraint_models.size();

    for (size_t ee_id = 0; ee_id < num_ee; ++ee_id)
    {
      const ConstraintModel & cmodel = constraint_models[ee_id];
      ConstraintData & cdata = constraint_datas[ee_id];

      cmodel.calc(model, data, cdata);
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename ForceMatrix,
    class ForceAllocator,
    ReferenceFrame rf>
  void mapConstraintForcesToJointForces(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<ForceMatrix> & constraint_forces,
    std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
    ReferenceFrameTag<rf> reference_frame)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_forces.size(), size_t(model.njoints));

    const Eigen::DenseIndex constraint_size = getTotalConstraintActiveSize(constraint_models);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), constraint_size);

    for (auto & force : joint_forces)
      force.setZero();

    Eigen::Index row_id = 0;
    for (size_t ee_id = 0; ee_id < constraint_models.size(); ++ee_id)
    {
      const ConstraintModel & cmodel = constraint_models[ee_id];
      const auto constraint_size = cmodel.activeSize();
      const ConstraintData & cdata = constraint_datas[ee_id];

      const auto constraint_force = constraint_forces.segment(row_id, constraint_size);
      cmodel.mapConstraintForceToJointForces(
        model, data, cdata, constraint_force, joint_forces, reference_frame);

      row_id += constraint_size;
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    class MotionAllocator,
    typename MotionMatrix,
    ReferenceFrame rf>
  void mapJointMotionsToConstraintMotions(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
    const Eigen::MatrixBase<MotionMatrix> & constraint_motions_,
    ReferenceFrameTag<rf> reference_frame)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_motions.size(), size_t(model.njoints));

    MotionMatrix & constraint_motions = constraint_motions_.const_cast_derived();
    const Eigen::DenseIndex constraint_size = getTotalConstraintActiveSize(constraint_models);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motions.rows(), constraint_size);

    Eigen::Index row_id = 0;
    for (size_t ee_id = 0; ee_id < constraint_models.size(); ++ee_id)
    {
      const ConstraintModel & cmodel = constraint_models[ee_id];
      const ConstraintData & cdata = constraint_datas[ee_id];
      const auto constraint_size = cmodel.activeSize();

      auto constraint_motion = constraint_motions.segment(row_id, constraint_size);
      cmodel.mapJointMotionsToConstraintMotion(
        model, data, cdata, joint_motions, constraint_motion, reference_frame);

      row_id += constraint_size;
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConstraintModel,
    typename ConstraintData,
    typename JacobianMatrixLike>
  void getConstraintJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintModelBase<ConstraintModel> & constraint_model_,
    ConstraintDataBase<ConstraintData> & constraint_data_,
    const Eigen::MatrixBase<JacobianMatrixLike> & J_)
  {
    JacobianMatrixLike & J = J_.const_cast_derived();
    const auto & constraint_model = constraint_model_.derived();
    auto & constraint_data = constraint_data_.derived();

    assert(model.check(data) && "data is not consistent with model.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.rows(), constraint_model.activeSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.cols(), model.nv);

    constraint_model.calc(model, data, constraint_data);
    constraint_model.jacobian(model, data, constraint_data, J);
  }

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
    const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & constraint_models,
    std::vector<Holder<ConstraintData>, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<DynamicMatrixLike> & J_)
  {
    typedef ConstraintModel ContraintModel;
    typedef ConstraintData ContraintData;

    const Eigen::DenseIndex constraint_size = getTotalConstraintActiveSize(constraint_models);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.rows(), constraint_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.cols(), model.nv);

    assert(model.check(data) && "data is not consistent with model.");
    assert(model.check(MimicChecker()) && "Function does not support mimic joints");

    DynamicMatrixLike & J = J_.const_cast_derived();
    Eigen::DenseIndex row_id = 0;
    for (size_t k = 0; k < constraint_models.size(); ++k)
    {
      const ContraintModel & cmodel = constraint_models[k];
      ContraintData & cdata = constraint_datas[k];

      getConstraintJacobian(model, data, cmodel, cdata, J.middleRows(row_id, cmodel.activeSize()));

      row_id += cmodel.size();
    }
  }

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
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<DynamicMatrixLike> & J_)
  {
    typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
    typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

    WrappedConstraintModelVector wrapped_constraint_models(
      constraint_models.cbegin(), constraint_models.cend());

    typedef std::reference_wrapper<ConstraintData> WrappedConstraintDataType;
    typedef std::vector<WrappedConstraintDataType> WrappedConstraintDataVector;

    WrappedConstraintDataVector wrapped_constraint_datas(
      constraint_datas.begin(), constraint_datas.end());

    getConstraintsJacobian(model, data, wrapped_constraint_models, wrapped_constraint_datas, J_);
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    template<typename T> class Holder,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  typename DataTpl<Scalar, Options, JointCollectionTpl>::MatrixXs getConstraintsJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & constraint_models,
    std::vector<Holder<ConstraintData>, ConstraintDataAllocator> & constraint_datas)
  {
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
    typedef typename Data::MatrixXs ReturnType;

    Eigen::DenseIndex constraint_size = 0;
    for (const ConstraintModel & cm : constraint_models)
      constraint_size += cm.size();

    ReturnType res = ReturnType::Zero(constraint_size, model.nv);
    getConstraintsJacobian(model, data, constraint_models, constraint_datas, res);

    return res;
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  typename DataTpl<Scalar, Options, JointCollectionTpl>::MatrixXs getConstraintsJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas)
  {
    typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
    typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

    WrappedConstraintModelVector wrapped_constraint_models(
      constraint_models.cbegin(), constraint_models.cend());

    typedef std::reference_wrapper<ConstraintData> WrappedConstraintDataType;
    typedef std::vector<WrappedConstraintDataType> WrappedConstraintDataVector;

    WrappedConstraintDataVector wrapped_constraint_datas(
      constraint_datas.begin(), constraint_datas.end());

    return getConstraintsJacobian(model, data, wrapped_constraint_models, wrapped_constraint_datas);
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    template<typename> class Holder,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename RhsMatrixType,
    typename ResultMatrixType>
  void evalConstraintJacobianTransposeMatrixProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & constraint_models,
    const std::vector<Holder<ConstraintData>, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<RhsMatrixType> & rhs,
    const Eigen::MatrixBase<ResultMatrixType> & res_)
  {

    const Eigen::DenseIndex constraint_size = getTotalConstraintActiveSize(constraint_models);
    ResultMatrixType & res = res_.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(rhs.rows(), constraint_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res_.rows(), model.nv);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res_.cols(), rhs.cols());

    Eigen::Index row_id = 0;
    res.setZero();
    for (size_t constraint_id = 0; constraint_id < constraint_models.size(); ++constraint_id)
    {
      const ConstraintModel & cmodel = constraint_models[constraint_id];
      const ConstraintData & cdata = constraint_datas[constraint_id];
      const auto constraint_size = cmodel.activeSize();

      const auto rhs_block = rhs.middleRows(row_id, constraint_size);
      cmodel.jacobianTransposeMatrixProduct(model, data, cdata, rhs_block, res, AddTo());

      row_id += constraint_size;
    }
  }

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
  void evalConstraintJacobianTransposeMatrixProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<RhsMatrixType> & rhs,
    const Eigen::MatrixBase<ResultMatrixType> & res_)
  {
    typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
    typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

    WrappedConstraintModelVector wrapped_constraint_models(
      constraint_models.cbegin(), constraint_models.cend());

    typedef std::reference_wrapper<const ConstraintData> WrappedConstraintDataType;
    typedef std::vector<WrappedConstraintDataType> WrappedConstraintDataVector;

    WrappedConstraintDataVector wrapped_constraint_datas(
      constraint_datas.begin(), constraint_datas.end());

    evalConstraintJacobianTransposeMatrixProduct(
      model, data, wrapped_constraint_models, wrapped_constraint_datas, rhs.derived(),
      res_.const_cast_derived());
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_contact_jacobian_hxx__
