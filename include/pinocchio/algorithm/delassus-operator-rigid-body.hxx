//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__

#include "pinocchio/algorithm/check.hpp"
#include "pinocchio/algorithm/contact-info.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/contact-jacobian.hpp"

namespace pinocchio
{

  template<typename DelassusOperator, typename ConfigVectorType>
  struct DelassusOperatorRigidBodyTplComputeForwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodyTplComputeForwardPass<DelassusOperator, ConfigVectorType>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::CustomData Data;

    typedef boost::fusion::vector<const Model &, Data &, const ConfigVectorType &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      Data & data,
      const Eigen::MatrixBase<ConfigVectorType> & q)
    {
      typedef typename Model::JointIndex JointIndex;

      const JointIndex i = jmodel.id();
      jmodel.calc(jdata.derived(), q.derived());

      const JointIndex parent = model.parents[i];
      data.liMi[i] = model.jointPlacements[i] * jdata.M();
      if (parent > 0)
        data.oMi[i] = data.oMi[parent] * data.liMi[i];
      else
        data.oMi[i] = data.liMi[i];
    }
  };

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodyTplComputeBackwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodyTplComputeBackwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::CustomData Data;
    typedef typename Model::Scalar Scalar;

    typedef boost::fusion::vector<const Model &, Data &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      Data & data)
    {
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Inertia Inertia;
      typedef typename JointModel::JointDataDerived JointData;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];
      typename Inertia::Matrix6 & Ia = data.Yaba[i];
      typename Inertia::Matrix6 & Ia_augmented = data.Yaba_augmented[i];

      JointData & jdata_augmented = boost::get<JointData>(data.joints_augmented[i]);

      jmodel.calc_aba(
        jdata.derived(), jmodel.jointVelocitySelector(model.armature), Ia, parent > 0);

      jmodel.calc_aba(
        jdata_augmented, jmodel.jointVelocitySelector(model.armature), Ia_augmented, parent > 0);

      if (parent > 0)
      {
        data.Yaba[parent] += impl::internal::SE3actOn<Scalar>::run(data.liMi[i], Ia);
        data.Yaba_augmented[parent] +=
          impl::internal::SE3actOn<Scalar>::run(data.liMi[i], Ia_augmented);
      }
    }
  };

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
    template<typename T>
    class Holder>
  template<typename ConfigVectorType>
  void DelassusOperatorRigidBodyTpl<Scalar, Options, JointCollectionTpl, Holder>::compute(
    const Eigen::MatrixBase<ConfigVectorType> & q)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model().nq, "The joint configuration vector is not of right size");

    const Model & model_ref = model();

    typedef DelassusOperatorRigidBodyTplComputeForwardPass<
      DelassusOperatorRigidBodyTpl, ConfigVectorType>
      Pass1;
    for (JointIndex i = 1; i < (JointIndex)model_ref.njoints; ++i)
    {
      typename Pass1::ArgsType args(model_ref, this->m_custom_data, q.derived());
      Pass1::run(model_ref.joints[i], this->m_custom_data.joints[i], args);
    }

    compute();
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
    template<typename T>
    class Holder>
  void DelassusOperatorRigidBodyTpl<Scalar, Options, JointCollectionTpl, Holder>::compute()
  {
    const Model & model_ref = model();
    Data & data_ref = data();
    CustomData & custom_data = this->m_custom_data;
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    ConstraintDataVector & constraint_datas_ref = constraint_datas();
    typedef typename Data::Vector3 Vector3;

    for (JointIndex i = 1; i < (JointIndex)model_ref.njoints; ++i)
    {
      custom_data.Yaba[i] = custom_data.Yaba_augmented[i] = model_ref.inertias[i].matrix();
    }

    // Append constraint inertia to Yaba_augmented
    for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
    {
      const RigidConstraintModel & cmodel = constraint_models_ref[ee_id];
      RigidConstraintData & cdata = constraint_datas_ref[ee_id];

      const Vector3 constraint_diagonal_inertia =
        this->m_damping_inverse.template segment<3>(Eigen::DenseIndex(ee_id * 3));

      typedef typename CustomData::Matrix6Vector InertiaAlignedVector;
      typedef typename InertiaAlignedVector::vector_base InertiaStdVector;
      cmodel.appendConstraintDiagonalInertiaToJointInertias(
        model_ref, data_ref, cdata, constraint_diagonal_inertia,
        static_cast<InertiaStdVector &>(custom_data.Yaba_augmented));
    }

    typedef DelassusOperatorRigidBodyTplComputeBackwardPass<DelassusOperatorRigidBodyTpl> Pass2;
    for (JointIndex i = JointIndex(model_ref.njoints - 1); i > 0; --i)
    {
      typename Pass2::ArgsType args(model_ref, this->m_custom_data);
      Pass2::run(model_ref.joints[i], this->m_custom_data.joints[i], args);
    }

    // Make a pass over the whole set of constraints to update the content
    {
      // TODO(jcarpent): change data_ref for custom_data
      evalConstraints(model_ref, data_ref, constraint_models_ref, constraint_datas_ref);
    }

    compute_conclude();
  }

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodyTplApplyOnTheRightBackwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodyTplApplyOnTheRightBackwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    //    typedef typename DelassusOperator::Data Data;
    typedef typename DelassusOperator::CustomData Data;

    typedef boost::fusion::vector<
      const Model &,
      //    Data &,
      Data &>
      ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      //                     Data & data
      Data & data)
    {
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Force Force;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      jmodel.jointVelocitySelector(data.u) =
        jdata.S().transpose() * data.f[i]; // The sign is switched compare to ABA

      if (parent > 0)
      {
        Force & pa = data.f[i];
        pa.toVector().noalias() -=
          jdata.UDinv()
          * jmodel.jointVelocitySelector(
            data.u); // The sign is switched compare to ABA as the sign of data.f[i] is switched too
        data.f[parent] += data.liMi[i].act(pa);
      }
    }
  };

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodyTplApplyOnTheRightForwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodyTplApplyOnTheRightForwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    //    typedef typename DelassusOperator::Data Data;
    typedef typename DelassusOperator::CustomData Data;

    typedef boost::fusion::vector<
      const Model &,
      //    Data &,
      Data &>
      ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      //                     Data & data
      Data & data)
    {
      typedef typename Model::JointIndex JointIndex;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      //      typename JointData::TangentVector_t ddq_joint;
      auto ddq_joint = jmodel.jointVelocitySelector(data.ddq);
      if (parent > 0)
      {
        data.a[i] += data.liMi[i].actInv(data.a[parent]);
        ddq_joint = jdata.Dinv() * jmodel.jointVelocitySelector(data.u)
                    - jdata.UDinv().transpose() * data.a[i].toVector();
        data.a[i] += jdata.S() * ddq_joint;
      }
      else
      {
        ddq_joint.noalias() = jdata.Dinv() * jmodel.jointVelocitySelector(data.u);
        data.a[i] = jdata.S() * ddq_joint;
      }
    }

  }; // struct DelassusOperatorRigidBodyTplApplyOnTheRightForwardPass

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
    template<typename T>
    class Holder>
  template<typename MatrixIn, typename MatrixOut>
  void DelassusOperatorRigidBodyTpl<Scalar, Options, JointCollectionTpl, Holder>::applyOnTheRight(
    const Eigen::MatrixBase<MatrixIn> & rhs, const Eigen::MatrixBase<MatrixOut> & res_) const
  {
    MatrixOut & res = res_.const_cast_derived();
    PINOCCHIO_CHECK_SAME_MATRIX_SIZE(rhs, res);

    const Model & model_ref = model();
    const Data & data_ref = data();
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    const ConstraintDataVector & constraint_datas_ref = constraint_datas();

    // Make a pass over the whole set of constraints to add the contributions of constraint forces
    mapConstraintForcesToJointForces(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, rhs, m_custom_data.f);

    // Backward sweep: propagate joint force contributions
    {
      typedef DelassusOperatorRigidBodyTplApplyOnTheRightBackwardPass<DelassusOperatorRigidBodyTpl>
        Pass1;
      typename Pass1::ArgsType args1(model_ref, this->m_custom_data);
      for (JointIndex i = JointIndex(model_ref.njoints - 1); i > 0; --i)
      {
        Pass1::run(model_ref.joints[i], this->m_custom_data.joints[i], args1);
      }
    }

    // Forward sweep: compute joint accelerations
    {
      typedef DelassusOperatorRigidBodyTplApplyOnTheRightForwardPass<DelassusOperatorRigidBodyTpl>
        Pass2;
      for (auto & motion : m_custom_data.a)
        motion.setZero();
      typename Pass2::ArgsType args2(model_ref, this->m_custom_data);
      for (JointIndex i = 1; i < JointIndex(model_ref.njoints); ++i)
      {
        Pass2::run(model_ref.joints[i], this->m_custom_data.joints[i], args2);
      }
    }

    // Make a pass over the whole set of constraints to project back the accelerations onto the
    // joint
    mapJointMotionsToConstraintMotions(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, this->m_custom_data.a, res);

    // Add damping contribution
    res.array() += m_damping.array() * rhs.array();
  }

  //  template<typename DelassusOperator>
  //  struct DelassusOperatorRigidBodyTplSolveInPlaceBackwardPass
  //  : public fusion::JointUnaryVisitorBase<
  //  DelassusOperatorRigidBodyTplSolveInPlaceBackwardPass<DelassusOperator> >
  //  {
  //    typedef typename DelassusOperator::Model Model;
  //      //    typedef typename DelassusOperator::Data Data;
  //    typedef typename DelassusOperator::CustomData Data;
  //
  //    typedef boost::fusion::vector<const Model &,
  //      //    Data &,
  //    Data &
  //    > ArgsType;
  //
  //    template<typename JointModel>
  //    static void algo(const pinocchio::JointModelBase<JointModel> & jmodel,
  //                     pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
  //                     const Model & model,
  //                     //                     Data & data
  //                     Data & data)
  //    {
  //      typedef typename Model::JointIndex JointIndex;
  //      typedef typename Data::Force Force;
  //
  //      const JointIndex i = jmodel.id();
  //      const JointIndex parent = model.parents[i];
  //
  //      jmodel.jointVelocitySelector(data.u) = jdata.S().transpose()*data.f[i]; // The sign is
  //      switched compare to ABA
  //
  //      if (parent > 0)
  //      {
  //        Force & pa = data.f[i];
  //        pa.toVector().noalias() -= jdata.UDinv() * jmodel.jointVelocitySelector(data.u); // The
  //        sign is switched compare to ABA as the sign of data.f[i] is switched too data.f[parent]
  //        += data.liMi[i].act(pa);
  //      }
  //    }
  //
  //  };

  template<
    typename Scalar,
    int Options,
    template<typename, int>
    class JointCollectionTpl,
    template<typename T>
    class Holder>
  template<typename MatrixLike>
  void DelassusOperatorRigidBodyTpl<Scalar, Options, JointCollectionTpl, Holder>::solveInPlace(
    const Eigen::MatrixBase<MatrixLike> & mat_) const
  {
    MatrixLike & mat = mat_.const_cast_derived();
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      mat.rows(), size(), "The input matrix does not match the size of the Delassus.");

    const Model & model_ref = model();
    const Data & data_ref = data();
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    const ConstraintDataVector & constraint_datas_ref = constraint_datas();

    mat.array() *= m_damping_inverse.array();

    // Make a pass over the whole set of constraints to add the contributions of constraint forces
    mapConstraintForcesToJointForces(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, mat, m_custom_data.f);

    // Backward sweep: propagate joint force contributions
    {
      typedef DelassusOperatorRigidBodyTplApplyOnTheRightBackwardPass<DelassusOperatorRigidBodyTpl>
        Pass1;
      typename Pass1::ArgsType args1(model_ref, this->m_custom_data);
      for (JointIndex i = JointIndex(model_ref.njoints - 1); i > 0; --i)
      {
        Pass1::run(model_ref.joints[i], this->m_custom_data.joints_augmented[i], args1);
      }
    }

    // Forward sweep: compute joint accelerations
    {
      typedef DelassusOperatorRigidBodyTplApplyOnTheRightForwardPass<DelassusOperatorRigidBodyTpl>
        Pass2;
      for (auto & motion : m_custom_data.a)
        motion.setZero();
      typename Pass2::ArgsType args2(model_ref, this->m_custom_data);
      for (JointIndex i = 1; i < JointIndex(model_ref.njoints); ++i)
      {
        Pass2::run(model_ref.joints[i], this->m_custom_data.joints_augmented[i], args2);
      }
    }

    // Make a pass over the whole set of constraints to project back the joint accelerations onto
    // the constraints
    mapJointMotionsToConstraintMotions(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, this->m_custom_data.a,
      this->m_custom_data.tmp_vec);

    mat.noalias() -= m_damping_inverse.asDiagonal() * this->m_custom_data.tmp_vec;
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__
