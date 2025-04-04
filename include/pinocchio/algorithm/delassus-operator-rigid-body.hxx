//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__

#include "pinocchio/algorithm/check.hpp"
#include "pinocchio/algorithm/constraints/constraints.hpp"

#include "pinocchio/algorithm/loop-constrained-aba.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/contact-jacobian.hpp"

namespace pinocchio
{

  template<typename DelassusOperator, typename ConfigVectorType>
  struct DelassusOperatorRigidBodySystemsComputeForwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodySystemsComputeForwardPass<DelassusOperator, ConfigVectorType>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;

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
      auto & oMi = data.oMi[i];
      if (parent > 0)
        oMi = data.oMi[parent] * data.liMi[i];
      else
        oMi = data.liMi[i];

      // ABA in WORLD frame requires these quantities
      jmodel.jointCols(data.J) = oMi.act(jdata.S());
    }
  };

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodySystemsComputeBackwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodySystemsComputeBackwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
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
      typedef std::pair<JointIndex, JointIndex> JointPair;

      const auto & neighbours = data.neighbour_links;
      auto & joint_cross_coupling = data.joint_cross_coupling;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      // ApplyOnTheRight
      {
        auto & Ia = data.Yaba[i];
        jmodel.calc_aba(
          jdata.derived(), jmodel.jointVelocitySelector(model.armature), Ia, parent > 0);
        if (parent > 0)
        {
          data.Yaba[parent] += impl::internal::SE3actOn<Scalar>::run(data.liMi[i], Ia);
        }
      }

      // SolveInPlace
      {
        JointData & _jdata_augmented = boost::get<JointData>(data.joints_augmented[i]);
        JointDataBase<JointData> & jdata_augmented =
          static_cast<JointDataBase<JointData> &>(_jdata_augmented);

        auto Jcols = jmodel.jointCols(data.J);
        auto & Ia_augmented = data.oYaba_augmented[i];

        jdata_augmented.U().noalias() = Ia_augmented * Jcols;
        jdata_augmented.StU().noalias() = Jcols.transpose() * jdata_augmented.U();

        // Account for the rotor inertia contribution
        jdata_augmented.StU().diagonal() += jmodel.jointVelocitySelector(model.armature);

        pinocchio::internal::PerformStYSInversion<Scalar>::run(
          jdata_augmented.StU(), jdata_augmented.Dinv());

        jdata_augmented.UDinv().noalias() = jdata_augmented.U() * jdata_augmented.Dinv();

        if (parent > 0)
        {
          Ia_augmented.noalias() -= jdata_augmented.UDinv() * jdata_augmented.U().transpose();
          data.oYaba_augmented[parent] += Ia_augmented;
        }
      }
    }
  };

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    template<typename T> class Holder>
  template<typename ConfigVectorType>
  void DelassusOperatorRigidBodySystemsTpl<
    Scalar,
    Options,
    JointCollectionTpl,
    ConstraintModel,
    Holder>::compute(const Eigen::MatrixBase<ConfigVectorType> & q)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model().nq, "The joint configuration vector is not of right size");

    const Model & model_ref = model();
    Data & data_ref = data();

    // Zero-th order forward kinematics
    typedef DelassusOperatorRigidBodySystemsComputeForwardPass<
      DelassusOperatorRigidBodySystemsTpl, ConfigVectorType>
      Pass1;
    for (JointIndex i = 1; i < (JointIndex)model_ref.njoints; ++i)
    {
      typename Pass1::ArgsType args(model_ref, data_ref, q.derived());
      Pass1::run(model_ref.joints[i], data_ref.joints[i], args);
    }

    compute();
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    template<typename T> class Holder>
  void DelassusOperatorRigidBodySystemsTpl<
    Scalar,
    Options,
    JointCollectionTpl,
    ConstraintModel,
    Holder>::compute()
  {
    typedef typename Data::Inertia Inertia;

    const Model & model_ref = model();
    Data & data_ref = data();
    //    CustomData & custom_data = this->m_custom_data;
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    ConstraintDataVector & constraint_datas_ref = constraint_datas();

    // Compute joint ordering for solveInPlace
    computeJointMinimalOrdering(model_ref, data_ref, constraint_models_ref);

    for (JointIndex i = 1; i < JointIndex(model_ref.njoints); ++i)
    {
      data_ref.Yaba[i] = model_ref.inertias[i].matrix();
      const Inertia oinertia = data_ref.oMi[i].act(model_ref.inertias[i]);
      data_ref.oYaba_augmented[i] = oinertia.matrix();
    }

    // Append constraint inertia to oYaba_augmented
    Eigen::Index row_id = 0;
    for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
    {
      const ConstraintModel & cmodel = constraint_models_ref[ee_id];
      ConstraintData & cdata = constraint_datas_ref[ee_id];

      const auto constraint_size = cmodel.size();

      const auto constraint_diagonal_inertia =
        this->m_sum_compliance_damping_inverse.segment(row_id, constraint_size);

      cmodel.calc(model_ref, data_ref, cdata);
      cmodel.appendCouplingConstraintInertias(
        model_ref, data_ref, cdata, constraint_diagonal_inertia, WorldFrameTag());

      row_id += constraint_size;
    }

    typedef DelassusOperatorRigidBodySystemsComputeBackwardPass<DelassusOperatorRigidBodySystemsTpl>
      Pass2;
    for (JointIndex i = JointIndex(model_ref.njoints - 1); i > 0; --i)
    {
      typename Pass2::ArgsType args(model_ref, data_ref);
      Pass2::run(model_ref.joints[i], data_ref.joints[i], args);
    }

    compute_conclude();
  }

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodySystemsTplApplyOnTheRightBackwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodySystemsTplApplyOnTheRightBackwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
    typedef typename DelassusOperator::CustomData CustomData;

    typedef boost::fusion::vector<const Model &, const Data &, CustomData &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      const pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      const Data & data,
      CustomData & custom_data)
    {
      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      jmodel.jointVelocitySelector(custom_data.u) -= jdata.S().transpose() * custom_data.f[i];

      if (parent > 0)
      {
        auto & pa = custom_data.f[i];
        pa.toVector().noalias() += jdata.UDinv() * jmodel.jointVelocitySelector(custom_data.u);
        custom_data.f[parent] += data.liMi[i].act(pa);
      }
    }
  };

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
    typedef typename DelassusOperator::CustomData CustomData;

    typedef boost::fusion::vector<const Model &, const Data &, CustomData &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      const pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      const Data & data,
      CustomData & custom_data)
    {
      typedef typename Model::JointIndex JointIndex;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      //      typename JointData::TangentVector_t ddq_joint;
      auto ddq_joint = jmodel.jointVelocitySelector(custom_data.ddq);
      if (parent > 0)
      {
        custom_data.a[i] += data.liMi[i].actInv(custom_data.a[parent]);
        ddq_joint = jdata.Dinv() * jmodel.jointVelocitySelector(custom_data.u)
                    - jdata.UDinv().transpose() * custom_data.a[i].toVector();
        custom_data.a[i] += jdata.S() * ddq_joint;
      }
      else
      {
        ddq_joint.noalias() = jdata.Dinv() * jmodel.jointVelocitySelector(custom_data.u);
        custom_data.a[i] = jdata.S() * ddq_joint;
      }
    }

  }; // struct DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    template<typename T> class Holder>
  template<typename MatrixIn, typename MatrixOut>
  void DelassusOperatorRigidBodySystemsTpl<
    Scalar,
    Options,
    JointCollectionTpl,
    ConstraintModel,
    Holder>::
    applyOnTheRight(
      const Eigen::MatrixBase<MatrixIn> & rhs, const Eigen::MatrixBase<MatrixOut> & res_) const
  {
    MatrixOut & res = res_.const_cast_derived();
    PINOCCHIO_CHECK_SAME_MATRIX_SIZE(rhs, res);

    const Model & model_ref = model();
    const Data & data_ref = data();
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    const ConstraintDataVector & constraint_datas_ref = constraint_datas();
    auto & custom_data = this->m_custom_data;

    // Make a pass over the whole set of constraints to add the contributions of constraint forces
    // mapConstraintForcesToJointForces(
    //   model_ref, data_ref, constraint_models_ref, constraint_datas_ref, rhs, m_custom_data.f);
    // TODO(jcarpent): extend the code to operator on matrices

    //    typedef Eigen::Map<VectorXs> MapVectorXs;
    //    MapVectorXs u = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, model_ref.nv, 1));
    {
      auto & u = custom_data.u;
      u.setZero();
      Eigen::Index row_id = 0;
      for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
      {
        const ConstraintModel & cmodel = constraint_models_ref[ee_id];
        const ConstraintData & cdata = constraint_datas_ref[ee_id];
        const auto csize = cmodel.size();

        cmodel.jacobianTransposeMatrixProduct(
          model_ref, data_ref, cdata, rhs.middleRows(row_id, csize), u, AddTo());

        row_id += csize;
      }
    }

    // Backward sweep: propagate joint force contributions
    for (auto & f : m_custom_data.f)
      f.setZero();
    {
      typedef DelassusOperatorRigidBodySystemsTplApplyOnTheRightBackwardPass<
        DelassusOperatorRigidBodySystemsTpl>
        Pass1;
      typename Pass1::ArgsType args1(model_ref, data_ref, custom_data);
      for (JointIndex i = JointIndex(model_ref.njoints - 1); i > 0; --i)
      {
        Pass1::run(model_ref.joints[i], data_ref.joints[i], args1);
      }
    }

    // Forward sweep: compute joint accelerations
    {
      typedef DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass<
        DelassusOperatorRigidBodySystemsTpl>
        Pass2;
      for (auto & motion : custom_data.a)
        motion.setZero();
      typename Pass2::ArgsType args2(model_ref, data_ref, custom_data);
      for (JointIndex i = 1; i < JointIndex(model_ref.njoints); ++i)
      {
        Pass2::run(model_ref.joints[i], data_ref.joints[i], args2);
      }
    }

    // Make a pass over the whole set of constraints to project back the accelerations onto the
    // joint
    //    mapJointMotionsToConstraintMotions(
    //      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, this->m_custom_data.a,
    //      res);

    // TODO(jcarpent): extend the code to operator on matrices
    {
      const auto & ddq = custom_data.ddq;
      Eigen::Index row_id = 0;
      for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
      {
        const ConstraintModel & cmodel = constraint_models_ref[ee_id];
        const ConstraintData & cdata = constraint_datas_ref[ee_id];
        const auto csize = cmodel.size();

        cmodel.jacobianMatrixProduct(
          model_ref, data_ref, cdata, ddq, res.middleRows(row_id, csize));

        row_id += csize;
      }
    }

    // Add damping contribution
    res.array() += m_sum_compliance_damping.array() * rhs.array();
  }

  //  template<typename DelassusOperator>
  //  struct DelassusOperatorRigidBodySystemsTplSolveInPlaceBackwardPass
  //  : public fusion::JointUnaryVisitorBase<
  //  DelassusOperatorRigidBodySystemsTplSolveInPlaceBackwardPass<DelassusOperator> >
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
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    template<typename T> class Holder>
  template<typename MatrixLike>
  void DelassusOperatorRigidBodySystemsTpl<
    Scalar,
    Options,
    JointCollectionTpl,
    ConstraintModel,
    Holder>::solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat_) const
  {
    MatrixLike & mat = mat_.const_cast_derived();
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      mat.rows(), size(), "The input matrix does not match the size of the Delassus.");

    const Model & model_ref = model();
    const Data & data_ref = data();
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    const ConstraintDataVector & constraint_datas_ref = constraint_datas();

    mat.array() *= m_sum_compliance_damping_inverse.array();

    // Make a pass over the whole set of constraints to add the contributions of constraint forces
    mapConstraintForcesToJointForces(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, mat, m_custom_data.f);

    // Backward sweep: propagate joint force contributions
    {
      typedef DelassusOperatorRigidBodySystemsTplApplyOnTheRightBackwardPass<
        DelassusOperatorRigidBodySystemsTpl>
        Pass1;
      typename Pass1::ArgsType args1(model_ref, this->m_custom_data);
      for (JointIndex i = JointIndex(model_ref.njoints - 1); i > 0; --i)
      {
        Pass1::run(model_ref.joints[i], this->m_custom_data.joints_augmented[i], args1);
      }
    }

    // Forward sweep: compute joint accelerations
    {
      typedef DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass<
        DelassusOperatorRigidBodySystemsTpl>
        Pass2;
      for (auto & motion : m_custom_data.a)
        motion.setZero();
      typename Pass2::ArgsType args2(model_ref, this->m_custom_data);
      for (JointIndex i = 1; i < JointIndex(model_ref.njoints); ++i)
      {
        Pass2::run(model_ref.joints[i], this->m_custom_data.joints_augmented[i], args2);
      }
    }

    typedef Eigen::Map<VectorXs> MapVectorXs;
    MapVectorXs tmp_vec = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, size(), 1));

    // Make a pass over the whole set of constraints to project back the joint accelerations onto
    // the constraints
    mapJointMotionsToConstraintMotions(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, this->m_custom_data.a,
      this->m_custom_data.tmp_vec);

    mat.noalias() -= m_sum_compliance_damping_inverse.asDiagonal() * this->m_custom_data.tmp_vec;
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__
