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

#include "pinocchio/algorithm/delassus-operator-rigid-body-visitors.hxx"

namespace pinocchio
{

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
    Holder>::compute(bool damping_compliance_update_only)
  {
    typedef typename Data::Inertia Inertia;

    const Model & model_ref = model();
    Data & data_ref = data();
    //    CustomData & custom_data = this->m_custom_data;
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    ConstraintDataVector & constraint_datas_ref = constraint_datas();

    // Compute joint ordering for solveInPlace
    if (!damping_compliance_update_only)
      computeJointMinimalOrdering(model_ref, data_ref, constraint_models_ref);

    for (JointIndex i = 1; i < JointIndex(model_ref.njoints); ++i)
    {
      const auto & joint_inertia = model_ref.inertias[i];
      if (!damping_compliance_update_only)
        data_ref.Yaba[i] = joint_inertia.matrix();
      const Inertia oinertia = data_ref.oMi[i].act(joint_inertia);
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

      if (!damping_compliance_update_only)
        cmodel.calc(model_ref, data_ref, cdata);
      cmodel.appendCouplingConstraintInertias(
        model_ref, data_ref, cdata, constraint_diagonal_inertia, WorldFrameTag());

      row_id += constraint_size;
    }

    if (damping_compliance_update_only)
    {
      typedef DelassusOperatorRigidBodySystemsComputeBackwardPass<
        DelassusOperatorRigidBodySystemsTpl, true>
        Pass2;
      for (const JointIndex i : data_ref.elimination_order)
      {
        typename Pass2::ArgsType args(model_ref, data_ref);
        Pass2::run(model_ref.joints[i], data_ref.joints[i], args);
      }
    }
    else
    {
      typedef DelassusOperatorRigidBodySystemsComputeBackwardPass<
        DelassusOperatorRigidBodySystemsTpl, false>
        Pass2;
      for (const JointIndex i : data_ref.elimination_order)
      {
        typename Pass2::ArgsType args(model_ref, data_ref);
        Pass2::run(model_ref.joints[i], data_ref.joints[i], args);
      }
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
        const auto rhs_rows = rhs.middleRows(row_id, csize);

        cmodel.jacobianTransposeMatrixProduct(model_ref, data_ref, cdata, rhs_rows, u, AddTo());

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

    //    PINOCCHIO_THROW_IF(
    //      m_dirty, std::logic_error,
    //      "The DelassusOperator has dirty quantities. Please call compute() method first.");
    if (isDirty())
      self_const_cast().compute(true);

    const Model & model_ref = model();
    const Data & data_ref = data();
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    const ConstraintDataVector & constraint_datas_ref = constraint_datas();

    mat.array() *= m_sum_compliance_damping_inverse.array();

    //    for(auto & of_augmented: m_custom_data.of_augmented)
    //      of_augmented.setZero();
    //
    //
    //    // Make a pass over the whole set of constraints to add the contributions of constraint
    //    forces mapConstraintForcesToJointForces(
    //      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, mat,
    //      m_custom_data.of_augmented);

    typedef Eigen::Map<VectorXs> MapVectorXs;
    MapVectorXs u = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, model_ref.nv, 1));

    {
      u.setZero();
      Eigen::Index row_id = 0;
      for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
      {
        const ConstraintModel & cmodel = constraint_models_ref[ee_id];
        const ConstraintData & cdata = constraint_datas_ref[ee_id];
        const auto csize = cmodel.size();
        const auto mat_rows = mat.middleRows(row_id, csize);

        cmodel.jacobianTransposeMatrixProduct(model_ref, data_ref, cdata, mat_rows, u, AddTo());

        row_id += csize;
      }
    }

    const auto & augmented_mass_matrix_operator = this->getAugmentedMassMatrixOperator();
    augmented_mass_matrix_operator.solveInPlace(u);

    typedef Eigen::Map<VectorXs> MapVectorXs;
    MapVectorXs tmp_vec = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, size(), 1));
    {
      Eigen::Index row_id = 0;
      for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
      {
        const ConstraintModel & cmodel = constraint_models_ref[ee_id];
        const ConstraintData & cdata = constraint_datas_ref[ee_id];
        const auto csize = cmodel.size();

        cmodel.jacobianMatrixProduct(
          model_ref, data_ref, cdata, u, tmp_vec.middleRows(row_id, csize));

        row_id += csize;
      }
    }

    // Make a pass over the whole set of constraints to project back the joint accelerations onto
    // the constraints
    //    mapJointMotionsToConstraintMotions(
    //      model_ref, data_ref, constraint_models_ref, constraint_datas_ref,
    //      this->m_custom_data.oa_augmented, tmp_vec);
    //

    mat.noalias() -= m_sum_compliance_damping_inverse.asDiagonal() * tmp_vec;
  }

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
    Holder>::AugmentedMassMatrixOperator::
    solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat_, bool reset_joint_force_vector) const
  {
    MatrixLike & mat = mat_.const_cast_derived();
    const auto & model_ref = m_self.model();
    const auto & data_ref = m_self.data();
    DelassusOperatorRigidBodySystemsTpl::CustomData & custom_data =
      const_cast<DelassusOperatorRigidBodySystemsTpl &>(m_self).getCustomData();
    const auto & elimination_order = data_ref.elimination_order;

    if (reset_joint_force_vector)
    {
      for (auto & of_augmented : custom_data.of_augmented)
        of_augmented.setZero();
    }

    // Backward sweep: propagate joint force contributions
    {
      custom_data.u = mat;
      typedef DelassusOperatorRigidBodySystemsTplSolveInPlaceBackwardPass<
        DelassusOperatorRigidBodySystemsTpl>
        Pass1;
      typename Pass1::ArgsType args1(model_ref, data_ref, custom_data);
      for (const JointIndex i : elimination_order)
      {
        Pass1::run(model_ref.joints[i], data_ref.joints_augmented[i], args1);
      }
    }

    // Forward sweep: compute joint accelerations
    {
      typedef DelassusOperatorRigidBodySystemsTplSolveInPlaceForwardPass<
        DelassusOperatorRigidBodySystemsTpl>
        Pass2;
      custom_data.oa_augmented[0].setZero();
      typename Pass2::ArgsType args2(model_ref, data_ref, custom_data);
      for (int it = int(elimination_order.size()) - 1; it >= 0; it--)
      {
        const JointIndex i = elimination_order[size_t(it)];
        Pass2::run(model_ref.joints[i], data_ref.joints_augmented[i], args2);
      }
    }

    mat = custom_data.ddq;
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__
