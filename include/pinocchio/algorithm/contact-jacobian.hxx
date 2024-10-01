//
// Copyright (c) 2021-2024 INRIA
//

#ifndef __pinocchio_algorithm_contact_jacobian_hxx__
#define __pinocchio_algorithm_contact_jacobian_hxx__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"

namespace pinocchio
{

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
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
    template<typename, int>
    class JointCollectionTpl,
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
    std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_forces.size(), size_t(model.njoints));

    const Eigen::DenseIndex constraint_size = getTotalConstraintSize(constraint_models);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), constraint_size);

    for (auto & force : joint_forces)
      force.setZero();

    for (size_t ee_id = 0; ee_id < constraint_models.size(); ++ee_id)
    {
      const ConstraintModel & cmodel = constraint_models[ee_id];
      const ConstraintData & cdata = constraint_datas[ee_id];

      const auto constraint_force =
        constraint_forces.template segment<3>(Eigen::DenseIndex(ee_id * 3));
      cmodel.mapConstraintForceToJointForces(model, data, cdata, constraint_force, joint_forces);
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
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
    const Eigen::MatrixBase<MotionMatrix> & constraint_motions_)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_motions.size(), size_t(model.njoints));

    MotionMatrix & constraint_motions = constraint_motions_.const_cast_derived();
    const Eigen::DenseIndex constraint_size = getTotalConstraintSize(constraint_models);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motions.rows(), constraint_size);

    for (size_t ee_id = 0; ee_id < constraint_models.size(); ++ee_id)
    {
      const ConstraintModel & cmodel = constraint_models[ee_id];
      const ConstraintData & cdata = constraint_datas[ee_id];

      auto constraint_motion = constraint_motions.template segment<3>(Eigen::DenseIndex(ee_id * 3));
      cmodel.mapJointMotionsToConstraintMotion(
        model, data, cdata, joint_motions, constraint_motion);
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
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
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.rows(), constraint_model.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.cols(), model.nv);

    constraint_model.calc(model, data, constraint_data);
    constraint_model.jacobian(model, data, constraint_data, J);
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
    template<typename T>
    class Holder,
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

    const Eigen::DenseIndex constraint_size = getTotalConstraintSize(constraint_models);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.rows(), constraint_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.cols(), model.nv);

    DynamicMatrixLike & J = J_.const_cast_derived();
    Eigen::DenseIndex row_id = 0;
    for (size_t k = 0; k < constraint_models.size(); ++k)
    {
      const ContraintModel & cmodel = constraint_models[k];
      ContraintData & cdata = constraint_datas[k];

      getConstraintJacobian(model, data, cmodel, cdata, J.middleRows(row_id, cmodel.size()));

      row_id += cmodel.size();
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
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
    template<typename, int>
    class JointCollectionTpl,
    typename ResultMatrixType>
  struct EvalConstraintJacobianTransposeProductBackwardPass
  : public fusion::JointUnaryVisitorBase<EvalConstraintJacobianTransposeProductBackwardPass<
      Scalar,
      Options,
      JointCollectionTpl,
      ResultMatrixType>>
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;

    typedef typename Data::Force Force;
    typedef typename PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(Force) ForceVector;

    typedef boost::fusion::vector<const Model &, const Data &, ForceVector &, ResultMatrixType &>
      ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      const pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      const Data & data,
      ForceVector & f,
      const Eigen::MatrixBase<ResultMatrixType> & res_)
    {
      typedef typename Model::JointIndex JointIndex;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      jmodel.jointVelocitySelector(res_.const_cast_derived()) = jdata.S().transpose() * f[i];

      if (parent > 0)
        f[parent] += data.liMi[i].act(f[i]);
    }

  }; // struct EvalConstraintJacobianTransposeProductBackwardPass

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
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
    const Eigen::MatrixBase<ResultMatrixType> & res_)
  {

    const Eigen::DenseIndex constraint_size = getTotalConstraintSize(constraint_models);
    ResultMatrixType & res = res_.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(rhs.rows(), constraint_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res_.rows(), model.nv);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res_.cols(), rhs.cols());

    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
    typedef typename Data::Force Force;
    typedef typename PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(Force) ForceVector;

    // Temporary memory variable
    // TODO(jcarpent): remove memory allocation here
    ForceVector joint_forces(size_t(model.njoints), Force::Zero());

    for (size_t ee_id = 0; ee_id < constraint_models.size(); ++ee_id)
    {
      const ConstraintModel & cmodel = constraint_models[ee_id];
      const ConstraintData & cdata = constraint_datas[ee_id];

      const auto constraint_force = rhs.template middleRows<3>(Eigen::DenseIndex(ee_id * 3));
      cmodel.mapConstraintForceToJointForces(model, data, cdata, constraint_force, joint_forces);
    }

    res.setZero();
    typedef EvalConstraintJacobianTransposeProductBackwardPass<
      Scalar, Options, JointCollectionTpl, ResultMatrixType>
      Pass;
    for (JointIndex i = (JointIndex)model.njoints - 1; i > 0; --i)
    {
      typename Pass::ArgsType arg(model, data, joint_forces, res);
      Pass::run(model.joints[i], data.joints[i], arg);
    }
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_contact_jacobian_hxx__
