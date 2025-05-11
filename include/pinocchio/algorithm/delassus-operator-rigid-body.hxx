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
  void DelassusOperatorRigidBodySystemsTpl<
    Scalar,
    Options,
    JointCollectionTpl,
    ConstraintModel,
    Holder>::
    update(
      const ConstraintModelVectorHolder & constraint_models_ref,
      const ConstraintDataVectorHolder & constraint_datas_ref)
  {
    m_constraint_models_ref = constraint_models_ref;
    m_constraint_datas_ref = constraint_datas_ref;

    computeJointMinimalOrdering(model(), data(), helper::get_ref(constraint_models_ref));
    m_dirty = true;
  }

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
    Holder>::
    compute(
      const Eigen::MatrixBase<ConfigVectorType> & q, bool apply_on_the_right, bool solve_in_place)
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

    compute(apply_on_the_right, solve_in_place);
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
    Holder>::compute_or_update_decomposition(bool apply_on_the_right, bool solve_in_place)
  {
    typedef typename Data::Inertia Inertia;
    using Matrix6 = typename Inertia::Matrix6;

    const Model & model_ref = model();
    Data & data_ref = data();
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    ConstraintDataVector & constraint_datas_ref = constraint_datas();

    //      computeJointMinimalOrdering(model_ref, data_ref, constraint_models_ref);

    for (JointIndex i = 1; i < JointIndex(model_ref.njoints); ++i)
    {
      const auto & joint_inertia = model_ref.inertias[i];
      if (apply_on_the_right)
        data_ref.Yaba[i] = joint_inertia.matrix();
      if (solve_in_place)
      {
        const Inertia oinertia = data_ref.oMi[i].act(joint_inertia);
        data_ref.oYaba_augmented[i] = oinertia.matrix();
      }
    }

    if (solve_in_place)
    {
      data_ref.joint_apparent_inertia = model_ref.armature;
      data_ref.joint_cross_coupling.apply([](Matrix6 & v) { v.setZero(); });

      // Append constraint inertia to oYaba_augmented
      Eigen::Index row_id = 0;
      for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
      {
        const auto & cmodel = helper::get_ref(constraint_models_ref[ee_id]);
        auto & cdata = helper::get_ref(constraint_datas_ref[ee_id]);

        const auto constraint_size = cmodel.size();

        const auto constraint_diagonal_inertia =
          this->m_sum_compliance_damping_inverse.segment(row_id, constraint_size);

        cmodel.appendCouplingConstraintInertias(
          model_ref, data_ref, cdata, constraint_diagonal_inertia, WorldFrameTag());

        row_id += constraint_size;
      }
    }

#define DO_PASS(apply_on_the_right_v, solve_in_place_v)                                            \
  {                                                                                                \
    typedef DelassusOperatorRigidBodySystemsComputeBackwardPass<                                   \
      DelassusOperatorRigidBodySystemsTpl, apply_on_the_right_v, solve_in_place_v>                 \
      Pass2;                                                                                       \
    for (const JointIndex i : data_ref.elimination_order)                                          \
    {                                                                                              \
      typename Pass2::ArgsType args(model_ref, data_ref);                                          \
      Pass2::run(model_ref.joints[i], data_ref.joints[i], args);                                   \
    }                                                                                              \
  }

    if (apply_on_the_right)
    {
      if (solve_in_place)
      {
        DO_PASS(true, true);
      }
      else
      {
        DO_PASS(true, false);
      }
    }
    else
    {
      if (solve_in_place)
      {
        DO_PASS(false, true);
      }
      else
      {
        DO_PASS(false, false);
      }
    }
#undef DO_PASS

    compute_conclude();
  }

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
    auto & u = custom_data.u;

    // Make a pass over the whole set of constraints to add the contributions of constraint forces
    mapConstraintForcesToJointSpace(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, rhs, m_custom_data.f, u,
      LocalFrameTag());
    // TODO(jcarpent): extend the code to operator on matrices

    //    typedef Eigen::Map<VectorXs> MapVectorXs;
    //    MapVectorXs u = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, model_ref.nv, 1));
    //    {
    //      auto & u = custom_data.u;
    //      u.setZero();
    //      Eigen::Index row_id = 0;
    //      for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
    //      {
    //        const auto & cmodel =
    //          helper::get_ref(constraint_models_ref[ee_id]);
    //        const auto & cdata =
    //          helper::get_ref(constraint_datas_ref[ee_id]);
    //        const auto csize = cmodel.size();
    //        const auto rhs_rows = rhs.middleRows(row_id, csize);
    //
    //        cmodel.jacobianTransposeMatrixProduct(model_ref, data_ref, cdata, rhs_rows, u,
    //        AddTo());
    //
    //        row_id += csize;
    //      }
    //    }

    // Backward sweep: propagate joint force contributions
    {
      //      for (auto & f : m_custom_data.f)
      //        f.setZero();
      //      auto & u = custom_data.u;
      //      u.setZero();

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
    mapJointSpaceToConstraintMotions(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, custom_data.a,
      custom_data.ddq, res, LocalFrameTag());

    // TODO(jcarpent): extend the code to operator on matrices
    //    {
    //      const auto & ddq = custom_data.ddq;
    //      Eigen::Index row_id = 0;
    //      for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
    //      {
    //        const auto & cmodel =
    //          helper::get_ref(constraint_models_ref[ee_id]);
    //        const auto & cdata =
    //          helper::get_ref(constraint_datas_ref[ee_id]);
    //        const auto csize = cmodel.size();
    //
    //        cmodel.jacobianMatrixProduct(
    //          model_ref, data_ref, cdata, ddq, res.middleRows(row_id, csize));
    //
    //        row_id += csize;
    //      }
    //    }

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
      self_const_cast().updateDecomposition();

    const Model & model_ref = model();
    const Data & data_ref = data();
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    const ConstraintDataVector & constraint_datas_ref = constraint_datas();
    auto & custom_data = this->m_custom_data;

    mat.array() *= m_sum_compliance_damping_inverse.array();

    // Make a pass over the whole set of constraints to add the contributions of constraint

    mapConstraintForcesToJointForces(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, mat,
      custom_data.of_augmented, WorldFrameTag());

    typedef Eigen::Map<VectorXs> MapVectorXs;
    MapVectorXs u = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, model_ref.nv, 1));
    u.setZero();
    //    {
    //      u.setZero();
    //      Eigen::Index row_id = 0;
    //      for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
    //      {
    //        const auto & cmodel =
    //          helper::get_ref(constraint_models_ref[ee_id]);
    //        const auto & cdata =
    //          helper::get_ref(constraint_datas_ref[ee_id]);
    //        const auto csize = cmodel.size();
    //        const auto mat_rows = mat.middleRows(row_id, csize);
    //
    //        cmodel.jacobianTransposeMatrixProduct(model_ref, data_ref, cdata, mat_rows, u,
    //        AddTo());
    //
    //        row_id += csize;
    //      }
    //    }

    const auto & augmented_mass_matrix_operator = this->getAugmentedMassMatrixOperator();
    augmented_mass_matrix_operator.solveInPlace(u, false);

    typedef Eigen::Map<VectorXs> MapVectorXs;
    MapVectorXs tmp_vec = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, size(), 1));
    //    {
    //      Eigen::Index row_id = 0;
    //      for (size_t ee_id = 0; ee_id < constraint_models_ref.size(); ++ee_id)
    //      {
    //        const auto & cmodel =
    //          helper::get_ref(constraint_models_ref[ee_id]);
    //        const auto & cdata =
    //          helper::get_ref(constraint_datas_ref[ee_id]);
    //        const auto csize = cmodel.size();
    //
    //        cmodel.jacobianMatrixProduct(
    //          model_ref, data_ref, cdata, u, tmp_vec.middleRows(row_id, csize));
    //
    //        row_id += csize;
    //      }
    //    }

    // Make a pass over the whole set of constraints to project back the joint accelerations onto
    // the constraints
    mapJointMotionsToConstraintMotions(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, custom_data.oa_augmented,
      tmp_vec, WorldFrameTag());

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
    auto & custom_data = const_cast<DelassusOperatorRigidBodySystemsTpl &>(m_self).getCustomData();
    const auto & elimination_order = data_ref.elimination_order;

    if (reset_joint_force_vector)
    {
      for (auto & of_augmented : custom_data.of_augmented)
        of_augmented.setZero();
    }

    // Backward sweep: propagate joint force contributions
    {
      custom_data.u = mat;
      typedef AugmentedMassMatrixOperatorSolveInPlaceBackwardPass<
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
      typedef AugmentedMassMatrixOperatorSolveInPlaceForwardPass<
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
