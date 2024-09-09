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
    class ConstraintModelAllocator,
    class ConstraintDataAllocator>
  void evalConstraints(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<RigidConstraintModelTpl<Scalar, Options>, ConstraintModelAllocator> &
      constraint_models,
    std::vector<RigidConstraintDataTpl<Scalar, Options>, ConstraintDataAllocator> &
      constraint_datas)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    const size_t num_ee = constraint_models.size();

    for (size_t ee_id = 0; ee_id < num_ee; ++ee_id)
    {
      const RigidConstraintModel & cmodel = constraint_models[ee_id];
      RigidConstraintData & cdata = constraint_datas[ee_id];

      cmodel.calc(model, data, cdata);
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
    class ConstraintModelAllocator,
    class ConstraintDataAllocator,
    typename ForceMatrix,
    class ForceAllocator>
  void mapConstraintForcesToJointForces(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<RigidConstraintModelTpl<Scalar, Options>, ConstraintModelAllocator> &
      constraint_models,
    const std::vector<RigidConstraintDataTpl<Scalar, Options>, ConstraintDataAllocator> &
      constraint_datas,
    const Eigen::MatrixBase<ForceMatrix> & constraint_forces,
    std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_forces.size(), size_t(model.njoints));

    const Eigen::DenseIndex constraint_size =
      Eigen::DenseIndex(getTotalConstraintSize(constraint_models));
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), constraint_size);

    for (auto & force : joint_forces)
      force.setZero();

    for (size_t ee_id = 0; ee_id < constraint_models.size(); ++ee_id)
    {
      const RigidConstraintModel & cmodel = constraint_models[ee_id];
      const RigidConstraintData & cdata = constraint_datas[ee_id];

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
    class ConstraintModelAllocator,
    class ConstraintDataAllocator,
    class MotionAllocator,
    typename MotionMatrix>
  void mapJointMotionsToConstraintMotions(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<RigidConstraintModelTpl<Scalar, Options>, ConstraintModelAllocator> &
      constraint_models,
    const std::vector<RigidConstraintDataTpl<Scalar, Options>, ConstraintDataAllocator> &
      constraint_datas,
    const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
    const Eigen::MatrixBase<MotionMatrix> & constraint_motions_)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_motions.size(), size_t(model.njoints));

    MotionMatrix & constraint_motions = constraint_motions_.const_cast_derived();
    const Eigen::DenseIndex constraint_size =
      Eigen::DenseIndex(getTotalConstraintSize(constraint_models));
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motions.rows(), constraint_size);

    for (size_t ee_id = 0; ee_id < constraint_models.size(); ++ee_id)
    {
      const RigidConstraintModel & cmodel = constraint_models[ee_id];
      const RigidConstraintData & cdata = constraint_datas[ee_id];

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
    typename Matrix6or3Like>
  void getConstraintJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const RigidConstraintModelTpl<Scalar, Options> & constraint_model,
    RigidConstraintDataTpl<Scalar, Options> & constraint_data,
    const Eigen::MatrixBase<Matrix6or3Like> & J_)
  {
    assert(model.check(data) && "data is not consistent with model.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.rows(), constraint_model.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.cols(), model.nv);

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::Motion Motion;
    typedef typename Model::SE3 SE3;
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;

    Matrix6or3Like & J = J_.const_cast_derived();

    typedef RigidConstraintModelTpl<Scalar, Options> ConstraintModel;
    const typename ConstraintModel::BooleanVector & colwise_joint1_sparsity =
      constraint_model.colwise_joint1_sparsity;
    const typename ConstraintModel::BooleanVector & colwise_joint2_sparsity =
      constraint_model.colwise_joint2_sparsity;
    const typename ConstraintModel::EigenIndexVector & colwise_span_indexes =
      constraint_model.colwise_span_indexes;

    SE3 & oMc1 = constraint_data.oMc1;
    oMc1 = data.oMi[constraint_model.joint1_id] * constraint_model.joint1_placement;
    SE3 & oMc2 = constraint_data.oMc2;
    oMc2 = data.oMi[constraint_model.joint2_id] * constraint_model.joint2_placement;
    SE3 & c1Mc2 = constraint_data.c1Mc2;
    c1Mc2 = oMc1.actInv(oMc2);

    for (size_t k = 0; k < colwise_span_indexes.size(); ++k)
    {
      const Eigen::DenseIndex col_id = colwise_span_indexes[k];

      const int sign = colwise_joint1_sparsity[col_id] != colwise_joint2_sparsity[col_id]
                         ? colwise_joint1_sparsity[col_id] ? +1 : -1
                         : 0; // specific case for CONTACT_3D

      typedef typename Data::Matrix6x::ConstColXpr ColXprIn;
      const ColXprIn Jcol_in = data.J.col(col_id);
      const MotionRef<const ColXprIn> Jcol_motion_in(Jcol_in);

      typedef typename Matrix6or3Like::ColXpr ColXprOut;
      ColXprOut Jcol_out = J.col(col_id);

      switch (constraint_model.type)
      {
      case CONTACT_3D: {
        switch (constraint_model.reference_frame)
        {
        case WORLD: {
          Jcol_out.noalias() = Jcol_motion_in.linear() * Scalar(sign);
          break;
        }
        case LOCAL: {
          if (sign == 0)
          {
            const Motion Jcol_local1(oMc1.actInv(Jcol_motion_in)); // TODO: simplify computations
            const Motion Jcol_local2(oMc2.actInv(Jcol_motion_in)); // TODO: simplify computations
            Jcol_out.noalias() = Jcol_local1.linear() - c1Mc2.rotation() * Jcol_local2.linear();
          }
          else if (sign == 1)
          {
            const Motion Jcol_local(oMc1.actInv(Jcol_motion_in));
            Jcol_out.noalias() = Jcol_local.linear();
          }
          else // sign == -1
          {
            const Motion Jcol_local(oMc2.actInv(Jcol_motion_in)); // TODO: simplify computations
            Jcol_out.noalias() =
              -c1Mc2.rotation() * Jcol_local.linear(); // TODO: simplify computations
          }
          break;
        }
        case LOCAL_WORLD_ALIGNED: {
          if (sign == 0)
          {
            Jcol_out.noalias() =
              (oMc2.translation() - oMc1.translation()).cross(Jcol_motion_in.angular());
          }
          else
          {
            if (sign == 1)
              Jcol_out.noalias() =
                Jcol_motion_in.linear() - oMc1.translation().cross(Jcol_motion_in.angular());
            else
              Jcol_out.noalias() =
                -Jcol_motion_in.linear() + oMc2.translation().cross(Jcol_motion_in.angular());
          }
          break;
        }
        }
        break;
      }
      case CONTACT_6D: {
        MotionRef<ColXprOut> Jcol_motion_out(Jcol_out);
        switch (constraint_model.reference_frame)
        {
        case WORLD: {
          Jcol_motion_out = Scalar(sign) * Jcol_motion_in;
          break;
        }
        case LOCAL: {
          Jcol_motion_out = Scalar(sign) * oMc1.actInv(Jcol_motion_in);
          break;
        }
        case LOCAL_WORLD_ALIGNED: {
          Motion Jcol_local_world_aligned(Jcol_motion_in);
          Jcol_local_world_aligned.linear() -=
            oMc1.translation().cross(Jcol_local_world_aligned.angular());
          Jcol_motion_out = Scalar(sign) * Jcol_local_world_aligned;
          break;
        }
        }
        break;
      }

      default:
        break;
      }
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
    template<typename T>
    class Holder,
    typename DynamicMatrixLike,
    class ConstraintModelAllocator,
    class ConstraintDataAllocator>
  void getConstraintsJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<
      Holder<const RigidConstraintModelTpl<Scalar, Options>>,
      ConstraintModelAllocator> & constraint_models,
    std::vector<Holder<RigidConstraintDataTpl<Scalar, Options>>, ConstraintDataAllocator> &
      constraint_datas,
    const Eigen::MatrixBase<DynamicMatrixLike> & J_)
  {
    typedef RigidConstraintModelTpl<Scalar, Options> ContraintModel;
    typedef RigidConstraintDataTpl<Scalar, Options> ContraintData;

    const Eigen::DenseIndex constraint_size =
      Eigen::DenseIndex(getTotalConstraintSize(constraint_models));
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
    class ConstraintModelAllocator,
    class ConstraintDataAllocator>
  void getConstraintsJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<RigidConstraintModelTpl<Scalar, Options>, ConstraintModelAllocator> &
      constraint_models,
    std::vector<RigidConstraintDataTpl<Scalar, Options>, ConstraintDataAllocator> &
      constraint_datas,
    const Eigen::MatrixBase<DynamicMatrixLike> & J_)
  {
    typedef std::reference_wrapper<const RigidConstraintModelTpl<Scalar, Options>>
      WrappedConstraintModelType;
    typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

    WrappedConstraintModelVector wrapped_constraint_models(
      constraint_models.cbegin(), constraint_models.cend());

    typedef std::reference_wrapper<RigidConstraintDataTpl<Scalar, Options>>
      WrappedConstraintDataType;
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
    class ConstraintModelAllocator,
    class ConstraintDataAllocator,
    typename RhsMatrixType,
    typename ResultMatrixType>
  void evalConstraintJacobianTransposeProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<RigidConstraintModelTpl<Scalar, Options>, ConstraintModelAllocator> &
      constraint_models,
    const std::vector<RigidConstraintDataTpl<Scalar, Options>, ConstraintDataAllocator> &
      constraint_datas,
    const Eigen::MatrixBase<RhsMatrixType> & rhs,
    const Eigen::MatrixBase<ResultMatrixType> & res_)
  {

    const Eigen::DenseIndex constraint_size =
      Eigen::DenseIndex(getTotalConstraintSize(constraint_models));
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
      const RigidConstraintModel & cmodel = constraint_models[ee_id];
      const RigidConstraintData & cdata = constraint_datas[ee_id];

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
