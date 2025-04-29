//
// Copyright (c) 2016-2020 CNRS INRIA
//

#ifndef __pinocchio_algorithm_joint_configuration_hxx__
#define __pinocchio_algorithm_joint_configuration_hxx__

#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/multibody/liegroup/liegroup-algo.hpp"
#include "pinocchio/multibody/liegroup/liegroup-joint.hpp"

/* --- Details -------------------------------------------------------------------- */
namespace pinocchio
{

  // --------------- API with return value as argument ---------------------- //

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentVectorType,
    typename ReturnType>
  void integrate(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType> & v,
    const Eigen::MatrixBase<ReturnType> & qout)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      v.size(), model.nv, "The joint velocity vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      qout.size(), model.nq, "The output argument is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;
    ReturnType & res = PINOCCHIO_EIGEN_CONST_CAST(ReturnType, qout);

    typedef IntegrateStep<LieGroup_t, ConfigVectorType, TangentVectorType, ReturnType> Algo;
    typename Algo::ArgsType args(q.derived(), v.derived(), res);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2,
    typename ReturnType>
  void interpolate(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorIn1> & q0,
    const Eigen::MatrixBase<ConfigVectorIn2> & q1,
    const Scalar & u,
    const Eigen::MatrixBase<ReturnType> & qout)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q0.size(), model.nq, "The first configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q1.size(), model.nq, "The second configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      qout.size(), model.nq, "The output argument is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;
    ReturnType & res = PINOCCHIO_EIGEN_CONST_CAST(ReturnType, qout);

    typedef InterpolateStep<LieGroup_t, ConfigVectorIn1, ConfigVectorIn2, Scalar, ReturnType> Algo;
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(
        model.joints[i], typename Algo::ArgsType(q0.derived(), q1.derived(), u, res.derived()));
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2,
    typename ReturnType>
  void difference(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorIn1> & q0,
    const Eigen::MatrixBase<ConfigVectorIn2> & q1,
    const Eigen::MatrixBase<ReturnType> & dvout)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q0.size(), model.nq, "The first configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q1.size(), model.nq, "The second configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      dvout.size(), model.nv, "The output argument is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;
    ReturnType & res = PINOCCHIO_EIGEN_CONST_CAST(ReturnType, dvout);

    typedef DifferenceStep<LieGroup_t, ConfigVectorIn1, ConfigVectorIn2, ReturnType> Algo;
    typename Algo::ArgsType args(q0.derived(), q1.derived(), res);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2,
    typename ReturnType>
  void squaredDistance(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorIn1> & q0,
    const Eigen::MatrixBase<ConfigVectorIn2> & q1,
    const Eigen::MatrixBase<ReturnType> & out)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q0.size(), model.nq, "The first configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q1.size(), model.nq, "The second configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      out.size(), (model.njoints - 1), "The output argument is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;
    ReturnType & distances = PINOCCHIO_EIGEN_CONST_CAST(ReturnType, out);

    typedef SquaredDistanceStep<LieGroup_t, ConfigVectorIn1, ConfigVectorIn2, ReturnType> Algo;
    for (JointIndex i = 0; i < (JointIndex)model.njoints - 1; ++i)
    {
      typename Algo::ArgsType args(i, q0.derived(), q1.derived(), distances.derived());
      Algo::run(model.joints[i + 1], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2,
    typename ReturnType>
  void randomConfiguration(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorIn1> & lowerLimits,
    const Eigen::MatrixBase<ConfigVectorIn2> & upperLimits,
    const Eigen::MatrixBase<ReturnType> & qout)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      lowerLimits.size(), model.nq, "The lower limits vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      upperLimits.size(), model.nq, "The upper limits vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      qout.size(), model.nq, "The output argument is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;
    ReturnType & q = PINOCCHIO_EIGEN_CONST_CAST(ReturnType, qout);

    typedef RandomConfigurationStep<LieGroup_t, ReturnType, ConfigVectorIn1, ConfigVectorIn2> Algo;
    typename Algo::ArgsType args(
      PINOCCHIO_EIGEN_CONST_CAST(ReturnType, q), lowerLimits.derived(), upperLimits.derived());
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ReturnType>
  void neutral(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ReturnType> & qout)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      qout.size(), model.nq, "The output argument is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;
    ReturnType & neutral_elt = PINOCCHIO_EIGEN_CONST_CAST(ReturnType, qout);

    typename NeutralStep<LieGroup_t, ReturnType>::ArgsType args(neutral_elt.derived());
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      NeutralStep<LieGroup_t, ReturnType>::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentVectorType,
    typename JacobianMatrixType>
  void dIntegrate(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType> & v,
    const Eigen::MatrixBase<JacobianMatrixType> & J,
    const ArgumentPosition arg,
    const AssignmentOperatorType op)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      v.size(), model.nv, "The joint velocity vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      J.rows(), model.nv, "The output argument is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      J.cols(), model.nv, "The output argument is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;

    typedef dIntegrateStep<LieGroup_t, ConfigVectorType, TangentVectorType, JacobianMatrixType>
      Algo;
    typename Algo::ArgsType args(
      q.derived(), v.derived(), PINOCCHIO_EIGEN_CONST_CAST(JacobianMatrixType, J), arg, op);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentMapMatrixType>
  void tangentMap(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentMapMatrixType> & TM,
    const AssignmentOperatorType op)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      TM.rows(), model.nq, "The output argument is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      TM.cols(), model.nv, "The output argument is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;

    typedef TangentMapStep<LieGroup_t, ConfigVectorType, TangentMapMatrixType> Algo;
    typename Algo::ArgsType args(
      q.derived(), PINOCCHIO_EIGEN_CONST_CAST(TangentMapMatrixType, TM), op);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentMapMatrixType>
  void compactTangentMap(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const std::vector<typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex> &
      joint_selection,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentMapMatrixType> & TMc)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The configuration vector is not of the right size");
    // Assume TMc.rows() == SUM_(j in joint_selection) j.nq() --> assert at the end
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      TMc.cols(), MAX_JOINT_NV, "The output argument is not of the right size");

    typedef CompactSetTangentMapStep<LieGroup_t, ConfigVectorType, TangentMapMatrixType> Algo;

    int idx = 0;
    typename Algo::ArgsType args(
      q.derived(), PINOCCHIO_EIGEN_CONST_CAST(TangentMapMatrixType, TMc), idx);
    for (size_t i = 0; i < joint_selection.size(); ++i)
    {
      Algo::run(model.joints[joint_selection[i]], args);
    }
    assert(idx == TMc.rows());
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename MatrixInType,
    typename MatrixOutType>
  void tangentMapProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<MatrixInType> & mat_in,
    const Eigen::MatrixBase<MatrixOutType> & mat_out,
    const AssignmentOperatorType op)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      mat_in.rows(), model.nv, "The input matrix is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      mat_out.rows(), model.nq, "The output matrix is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      mat_in.cols(), mat_out.cols(), "The input/output matrix sized do not match");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;

    typedef TangentMapProductStep<LieGroup_t, ConfigVectorType, MatrixInType, MatrixOutType> Algo;
    typename Algo::ArgsType args(
      q.derived(), mat_in.derived(), PINOCCHIO_EIGEN_CONST_CAST(MatrixOutType, mat_out), op);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename MatrixInType,
    typename MatrixOutType>
  void tangentMapTransposeProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<MatrixInType> & mat_in,
    const Eigen::MatrixBase<MatrixOutType> & mat_out,
    const AssignmentOperatorType op)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      mat_in.rows(), model.nq, "The input matrix is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      mat_out.rows(), model.nv, "The output matrix is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      mat_in.cols(), mat_out.cols(), "The input/output matrix sized do not match");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;

    typedef TangentMapTransposeProductStep<
      LieGroup_t, ConfigVectorType, MatrixInType, MatrixOutType>
      Algo;
    typename Algo::ArgsType args(
      q.derived(), mat_in.derived(), PINOCCHIO_EIGEN_CONST_CAST(MatrixOutType, mat_out), op);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentVectorType,
    typename JacobianMatrixType1,
    typename JacobianMatrixType2>
  void dIntegrateTransport(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType> & v,
    const Eigen::MatrixBase<JacobianMatrixType1> & Jin,
    const Eigen::MatrixBase<JacobianMatrixType2> & Jout,
    const ArgumentPosition arg)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      v.size(), model.nv, "The joint velocity vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      Jin.rows(), model.nv, "The input matrix is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      Jout.rows(), Jin.rows(), "The output argument should be the same size as input matrix");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      Jout.cols(), Jin.cols(), "The output argument should be the same size as input matrix");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;

    typedef dIntegrateTransportStep<
      LieGroup_t, ConfigVectorType, TangentVectorType, JacobianMatrixType1, JacobianMatrixType2>
      Algo;
    typename Algo::ArgsType args(
      q.derived(), v.derived(), Jin.derived(),
      PINOCCHIO_EIGEN_CONST_CAST(JacobianMatrixType2, Jout), arg);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentVectorType,
    typename JacobianMatrixType>
  void dIntegrateTransport(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType> & v,
    const Eigen::MatrixBase<JacobianMatrixType> & J,
    const ArgumentPosition arg)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      v.size(), model.nv, "The joint velocity vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J.rows(), model.nv, "The input matrix is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;

    typedef dIntegrateTransportInPlaceStep<
      LieGroup_t, ConfigVectorType, TangentVectorType, JacobianMatrixType>
      Algo;
    typename Algo::ArgsType args(
      q.derived(), v.derived(), PINOCCHIO_EIGEN_CONST_CAST(JacobianMatrixType, J), arg);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVector1,
    typename ConfigVector2,
    typename JacobianMatrix>
  void dDifference(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVector1> & q0,
    const Eigen::MatrixBase<ConfigVector2> & q1,
    const Eigen::MatrixBase<JacobianMatrix> & J,
    const ArgumentPosition arg)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q0.size(), model.nq, "The configuration vector q0 is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q1.size(), model.nq, "The configuration vector q1 is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      J.rows(), model.nv, "The output argument is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      J.cols(), model.nv, "The output argument is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;

    typedef dDifferenceStep<LieGroup_t, ConfigVector1, ConfigVector2, JacobianMatrix> Algo;
    typename Algo::ArgsType args(
      q0.derived(), q1.derived(), PINOCCHIO_EIGEN_CONST_CAST(JacobianMatrix, J), arg);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2>
  Scalar squaredDistanceSum(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorIn1> & q0,
    const Eigen::MatrixBase<ConfigVectorIn2> & q1)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q0.size(), model.nq, "The first configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q1.size(), model.nq, "The second configuration vector is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;
    typename ConfigVectorIn1::Scalar squaredDistance = Scalar(0.0);

    typedef SquaredDistanceSumStep<LieGroup_t, ConfigVectorIn1, ConfigVectorIn2, Scalar> Algo;
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      typename Algo::ArgsType args(q0.derived(), q1.derived(), squaredDistance);
      Algo::run(model.joints[i], args);
    }

    return squaredDistance;
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2>
  Scalar distance(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorIn1> & q0,
    const Eigen::MatrixBase<ConfigVectorIn2> & q1)
  {
    const Scalar & squaredDistance = squaredDistanceSum<
      LieGroup_t, Scalar, Options, JointCollectionTpl, ConfigVectorIn1, ConfigVectorIn2>(
      model, q0.derived(), q1.derived());
    return math::sqrt(squaredDistance);
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType>
  void normalize(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorType> & qout)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      qout.size(), model.nq, "The output argument is not of the right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;

    typedef NormalizeStep<LieGroup_t, ConfigVectorType> Algo;
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(
        model.joints[i],
        typename Algo::ArgsType(PINOCCHIO_EIGEN_CONST_CAST(ConfigVectorType, qout)));
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn>
  bool isNormalized(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorIn> & q,
    const Scalar & prec)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The configuration vector is not of the right size");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      (check_expression_if_real<Scalar, true>(prec >= 0)), "The precision should be positive");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;

    bool result = true;
    typedef IsNormalizedStep<LieGroup_t, ConfigVectorIn, Scalar> Algo;
    typename Algo::ArgsType args(q.derived(), prec, result);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
      if (!result)
        return false;
    }

    return true;
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2>
  bool isSameConfiguration(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorIn1> & q1,
    const Eigen::MatrixBase<ConfigVectorIn2> & q2,
    const Scalar & prec)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q1.size(), model.nq, "The first configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q2.size(), model.nq, "The second configuration vector is not of the right size");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      (check_expression_if_real<Scalar, true>(prec >= 0)), "The precision should be positive");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;

    bool result = true;
    typedef IsSameConfigurationStep<LieGroup_t, ConfigVectorIn1, ConfigVectorIn2, Scalar> Algo;
    typename Algo::ArgsType args(result, q1.derived(), q2.derived(), prec);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
      if (!result)
        return false;
    }

    return true;
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVector,
    typename JacobianMatrix>
  void integrateCoeffWiseJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVector> & q,
    const Eigen::MatrixBase<JacobianMatrix> & jacobian)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The configuration vector is not of the right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(jacobian.rows(), model.nq);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian.cols(), model.nv, "The jacobian does not have the right dimension");

    typedef IntegrateCoeffWiseJacobianStep<LieGroup_t, ConfigVector, JacobianMatrix> Algo;
    typename Algo::ArgsType args(q.derived(), PINOCCHIO_EIGEN_CONST_CAST(JacobianMatrix, jacobian));
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl>
  void lieGroup(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    typename LieGroup_t::template operationProduct<Scalar, Options>::type & lgo)
  {

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef LieGroupInstanceStep<LieGroup_t, Scalar, Options> Algo;
    typedef typename Model::JointIndex JointIndex;

    typename Algo::ArgsType args(lgo);
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Algo::run(model.joints[i], args);
    }
  }

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  void indexvInfo(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const std::vector<typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex> &
      joint_selection,
    std::vector<int> & nvs,
    std::vector<int> & idx_vs)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(nvs.size(), 0, "The nvs vector must empty");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(idx_vs.size(), 0, "The nvs vector must empty");

    typename IndexvInfoStep::ArgsType args(nvs, idx_vs);

    for (size_t i = 0; i < joint_selection.size(); ++i)
    {
      IndexvInfoStep::run(model.joints[joint_selection[i]], args);
    }
  }

  // ----------------- API that allocates memory ---------------------------- //

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentVectorType>
  typename PINOCCHIO_EIGEN_PLAIN_TYPE(ConfigVectorType) integrate(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType> & v)
  {
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(ConfigVectorType) ReturnType;
    ReturnType res(model.nq);
    integrate<
      LieGroup_t, Scalar, Options, JointCollectionTpl, ConfigVectorType, TangentVectorType,
      ReturnType>(model, q.derived(), v.derived(), res);
    return res;
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2>
  typename PINOCCHIO_EIGEN_PLAIN_TYPE(ConfigVectorIn1) interpolate(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorIn1> & q0,
    const Eigen::MatrixBase<ConfigVectorIn2> & q1,
    const Scalar & u)
  {
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(ConfigVectorIn1) ReturnType;
    ReturnType res(model.nq);
    interpolate<
      LieGroup_t, Scalar, Options, JointCollectionTpl, ConfigVectorIn1, ConfigVectorIn2,
      ReturnType>(model, q0.derived(), q1.derived(), u, res);
    return res;
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2>
  typename PINOCCHIO_EIGEN_PLAIN_TYPE(ConfigVectorIn1) difference(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorIn1> & q0,
    const Eigen::MatrixBase<ConfigVectorIn2> & q1)
  {
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(ConfigVectorIn1) ReturnType;
    ReturnType res(model.nv);
    difference<
      LieGroup_t, Scalar, Options, JointCollectionTpl, ConfigVectorIn1, ConfigVectorIn2,
      ReturnType>(model, q0.derived(), q1.derived(), res);
    return res;
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2>
  typename PINOCCHIO_EIGEN_PLAIN_TYPE(ConfigVectorIn1) squaredDistance(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const Eigen::MatrixBase<ConfigVectorIn1> & q0,
    const Eigen::MatrixBase<ConfigVectorIn2> & q1)
  {
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(ConfigVectorIn1) ReturnType;
    ReturnType distances(ReturnType::Zero(model.njoints - 1));
    squaredDistance<
      LieGroup_t, Scalar, Options, JointCollectionTpl, ConfigVectorIn1, ConfigVectorIn2,
      ReturnType>(model, q0.derived(), q1.derived(), distances);
    return distances;
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2>
  typename PINOCCHIO_EIGEN_PLAIN_TYPE_NO_PARENS(
    (typename ModelTpl<Scalar, Options, JointCollectionTpl>::ConfigVectorType))
    randomConfiguration(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const Eigen::MatrixBase<ConfigVectorIn1> & lowerLimits,
      const Eigen::MatrixBase<ConfigVectorIn2> & upperLimits)
  {
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE_NO_PARENS(
      (typename ModelTpl<Scalar, Options, JointCollectionTpl>::ConfigVectorType)) ReturnType;
    ReturnType q(model.nq);
    randomConfiguration<
      LieGroup_t, Scalar, Options, JointCollectionTpl, ConfigVectorIn1, ConfigVectorIn2,
      ReturnType>(model, lowerLimits.derived(), upperLimits.derived(), q);
    return q;
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl>
  typename PINOCCHIO_EIGEN_PLAIN_TYPE_NO_PARENS(
    (typename ModelTpl<Scalar, Options, JointCollectionTpl>::ConfigVectorType))
    randomConfiguration(const ModelTpl<Scalar, Options, JointCollectionTpl> & model)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::ConfigVectorType ConfigVectorType;
    return randomConfiguration<
      LieGroup_t, Scalar, Options, JointCollectionTpl, ConfigVectorType, ConfigVectorType>(
      model, model.lowerPositionLimit, model.upperPositionLimit);
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl>
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options>
  neutral(const ModelTpl<Scalar, Options, JointCollectionTpl> & model)
  {
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> ReturnType;
    ReturnType q(model.nq);
    neutral<LieGroup_t, Scalar, Options, JointCollectionTpl, ReturnType>(model, q);
    return q;
  }

  template<
    typename LieGroup_t,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl>
  typename LieGroup_t::template operationProduct<Scalar, Options>::type
  lieGroup(const ModelTpl<Scalar, Options, JointCollectionTpl> & model)
  {
    typedef typename LieGroup_t::template operationProduct<Scalar, Options>::type LGO;
    LGO lgo;
    lieGroup<LieGroup_t, Scalar, Options, JointCollectionTpl>(model, lgo);
    return lgo;
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_joint_configuration_hxx__
