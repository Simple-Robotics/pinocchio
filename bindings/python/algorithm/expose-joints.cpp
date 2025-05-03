//
// Copyright (c) 2015-2021 CNRS INRIA
//

#include "pinocchio/bindings/python/algorithm/algorithms.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

namespace pinocchio
{
  namespace python
  {

    typedef typename LieGroupMap::template operationProduct<context::Scalar, context::Options>::type
      LgType;

    static context::VectorXs
    normalize_proxy(const context::Model & model, const context::VectorXs & config)
    {
      context::VectorXs q(config);
      normalize(model, q);
      return q;
    }

    static context::VectorXs randomConfiguration_proxy(const context::Model & model)
    {
      return randomConfiguration(model);
    }

    bp::tuple dIntegrate_proxy(
      const context::Model & model, const context::VectorXs & q, const context::VectorXs & v)
    {
      context::MatrixXs J0(context::MatrixXs::Zero(model.nv, model.nv));
      context::MatrixXs J1(context::MatrixXs::Zero(model.nv, model.nv));

      dIntegrate(model, q, v, J0, ARG0);
      dIntegrate(model, q, v, J1, ARG1);

      return bp::make_tuple(J0, J1);
    }

    context::MatrixXs dIntegrate_arg_proxy(
      const context::Model & model,
      const context::VectorXs & q,
      const context::VectorXs & v,
      const ArgumentPosition arg)
    {
      context::MatrixXs J(context::MatrixXs::Zero(model.nv, model.nv));

      dIntegrate(model, q, v, J, arg, SETTO);

      return J;
    }

    context::MatrixXs dIntegrateTransport_proxy(
      const context::Model & model,
      const context::VectorXs & q,
      const context::VectorXs & v,
      const context::MatrixXs & Jin,
      const ArgumentPosition arg)
    {
      int ncols = Jin.cols();
      context::MatrixXs Jout(context::MatrixXs::Zero(model.nv, ncols));
      dIntegrateTransport(model, q, v, Jin, Jout, arg);
      return Jout;
    }

    bp::tuple dDifference_proxy(
      const context::Model & model, const context::VectorXs & q1, const context::VectorXs & q2)
    {
      context::MatrixXs J0(context::MatrixXs::Zero(model.nv, model.nv));
      context::MatrixXs J1(context::MatrixXs::Zero(model.nv, model.nv));

      dDifference(model, q1, q2, J0, ARG0);
      dDifference(model, q1, q2, J1, ARG1);

      return bp::make_tuple(J0, J1);
    }

    context::MatrixXs dDifference_arg_proxy(
      const context::Model & model,
      const context::VectorXs & q1,
      const context::VectorXs & q2,
      const ArgumentPosition arg)
    {
      context::MatrixXs J(context::MatrixXs::Zero(model.nv, model.nv));

      dDifference(model, q1, q2, J, arg);

      return J;
    }

    context::MatrixXs tangentMap_proxy(const context::Model & model, const context::VectorXs & q)
    {
      context::MatrixXs TM(context::MatrixXs::Zero(model.nq, model.nv));

      tangentMap(model, q, TM, SETTO);

      return TM;
    }

    context::MatrixXs
    compactTangentMap_proxy(const context::Model & model, const context::VectorXs & q)
    {
      context::MatrixXs TMc(context::MatrixXs::Zero(model.nq, MAX_JOINT_NV));
      typedef typename context::Model::JointIndex JointIndex;
      std::vector<JointIndex> joint_selection;
      for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
      {
        joint_selection.push_back(i);
      }

      compactTangentMap(model, joint_selection, q, TMc);

      return TMc;
    }

    context::MatrixXs tangentMapProduct_proxy(
      const context::Model & model, const context::VectorXs & q, const context::MatrixXs & mat_in)
    {
      context::MatrixXs mat_out(context::MatrixXs::Zero(model.nq, mat_in.cols()));

      tangentMapProduct(model, q, mat_in, mat_out, SETTO);

      return mat_out;
    }

    context::MatrixXs tangentMapTransposeProduct_proxy(
      const context::Model & model, const context::VectorXs & q, const context::MatrixXs & mat_in)
    {
      context::MatrixXs mat_out(context::MatrixXs::Zero(model.nv, mat_in.cols()));

      tangentMapTransposeProduct(model, q, mat_in, mat_out, SETTO);

      return mat_out;
    }

    LgType lieGroup_proxy(const context::Model & model)
    {
      LgType res;

      lieGroup(model, res);

      return res;
    }

    bp::tuple indexvInfo_proxy(const context::Model & model)
    {
      std::vector<int> nvs;
      std::vector<int> idx_vs;

      typedef typename context::Model::JointIndex JointIndex;
      std::vector<JointIndex> joint_selection;
      for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
      {
        joint_selection.push_back(i);
      }

      indexvInfo(model, joint_selection, nvs, idx_vs);

      return bp::make_tuple(nvs, idx_vs);
    }

    void exposeJointsAlgo()
    {
      typedef context::Scalar Scalar;
      typedef context::VectorXs VectorXs;
      enum
      {
        Options = context::Options
      };

      bp::def(
        "integrate", &integrate<Scalar, Options, JointCollectionDefaultTpl, VectorXs, VectorXs>,
        bp::args("model", "q", "v"),
        "Integrate the joint configuration vector q with a tangent vector v during one unit time.\n"
        "This is the canonical integrator of a Configuration Space composed of Lie groups, such as "
        "most robots.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq: the joint configuration vector (size model.nq)\n"
        "\tv: the joint velocity vector (size model.nv)\n");

      bp::def(
        "dIntegrate", &dIntegrate_proxy, bp::args("model", "q", "v"),
        "Computes the partial derivatives of the integrate function with respect to the first "
        "and the second argument, and returns the two Jacobians as a tuple.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq: the joint configuration vector (size model.nq)\n"
        "\tv: the joint velocity vector (size model.nv)\n");

      bp::def(
        "dIntegrate", &dIntegrate_arg_proxy, bp::args("model", "q", "v", "argument_position"),
        "Computes the partial derivatives of the integrate function with respect to the first (arg "
        "== ARG0) "
        "or the second argument (arg == ARG1).\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq: the joint configuration vector (size model.nq)\n"
        "\tv: the joint velocity vector (size model.nv)\n"
        "\targument_position: either pinocchio.ArgumentPosition.ARG0 or "
        "pinocchio.ArgumentPosition.ARG1, depending on the desired Jacobian value.\n");

      bp::def(
        "dIntegrateTransport", &dIntegrateTransport_proxy,
        bp::args("model", "q", "v", "Jin", "argument_position"),
        "Takes a matrix expressed at q (+) v and uses parallel transport to express it in the "
        "tangent space at q."
        "\tThis operation does the product of the matrix by the Jacobian of the integration "
        "operation, but more efficiently."
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq: the joint configuration vector (size model.nq)\n"
        "\tv: the joint velocity vector (size model.nv)\n"
        "\tJin: the input matrix (row size model.nv)"
        "\targument_position: either pinocchio.ArgumentPosition.ARG0 (q) or "
        "pinocchio.ArgumentPosition.ARG1 (v), depending on the desired Jacobian value.\n");

      bp::def(
        "interpolate", &interpolate<Scalar, Options, JointCollectionDefaultTpl, VectorXs, VectorXs>,
        bp::args("model", "q1", "q2", "alpha"),
        "Interpolate between two given joint configuration vectors q1 and q2.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq1: the initial joint configuration vector (size model.nq)\n"
        "\tq2: the terminal joint configuration vector (size model.nq)\n"
        "\talpha: the interpolation coefficient in [0,1]\n");

      bp::def(
        "difference", &difference<Scalar, Options, JointCollectionDefaultTpl, VectorXs, VectorXs>,
        bp::args("model", "q1", "q2"),
        "Difference between two joint configuration vectors, i.e. the tangent vector that must be "
        "integrated during one unit time"
        "to go from q1 to q2.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq1: the initial joint configuration vector (size model.nq)\n"
        "\tq2: the terminal joint configuration vector (size model.nq)\n");

      bp::def(
        "squaredDistance",
        &squaredDistance<Scalar, Options, JointCollectionDefaultTpl, VectorXs, VectorXs>,
        bp::args("model", "q1", "q2"),
        "Squared distance vector between two joint configuration vectors.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq1: the initial joint configuration vector (size model.nq)\n"
        "\tq2: the terminal joint configuration vector (size model.nq)\n");

      bp::def(
        "distance", &distance<Scalar, Options, JointCollectionDefaultTpl, VectorXs, VectorXs>,
        bp::args("model", "q1", "q2"),
        "Distance between two joint configuration vectors.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq1: the initial joint configuration vector (size model.nq)\n"
        "\tq2: the terminal joint configuration vector (size model.nq)\n");

      bp::def(
        "dDifference", &dDifference_proxy, bp::args("model", "q1", "q2"),
        "Computes the partial derivatives of the difference function with respect to the first "
        "and the second argument, and returns the two Jacobians as a tuple.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq1: the initial joint configuration vector (size model.nq)\n"
        "\tq2: the terminal joint configuration vector (size model.nq)\n");

      bp::def(
        "dDifference", &dDifference_arg_proxy, bp::args("model", "q1", "q2", "argument_position"),
        "Computes the partial derivatives of the difference function with respect to the first "
        "(arg == ARG0) "
        "or the second argument (arg == ARG1).\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq1: the initial joint configuration vector (size model.nq)\n"
        "\tq2: the terminal joint configuration vector (size model.nq)\n"
        "\targument_position: either pinocchio.ArgumentPosition.ARG0 or "
        "pinocchio.ArgumentPosition.ARG1, depending on the desired Jacobian value.\n");

      bp::def(
        "tangentMap", &tangentMap_proxy, bp::args("model", "q"),
        "Computes the tangent map in configuration q that map of a small variation express in the "
        "Lie algebra as a small variation in the parametric space.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq: the joint configuration vector (size model.nq)\n");

      bp::def(
        "compactTangentMap", &compactTangentMap_proxy, bp::args("model", "q"),
        "Computes the tangent map in configuration q that map of a small variation express in the "
        "Lie algebra as a small variation in the parametric space. Store the result in a compact "
        "manner that can be exploited using indexvInfo.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq: the joint configuration vector (size model.nq)\n");

      bp::def(
        "tangentMapProduct", &tangentMapProduct_proxy, bp::args("model", "q", "mat_in"),
        "Apply the tangent map to a matrix mat_in.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq: the joint configuration vector (size model.nq)\n"
        "\tmat_in: a matrix (size model.nq, ncols)");

      bp::def(
        "tangentMapTransposeProduct", &tangentMapTransposeProduct_proxy,
        bp::args("model", "q", "mat_in"),
        "Apply the tangent map to a matrix mat_in.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq: the joint configuration vector (size model.nq)\n"
        "\tmat_in: a matrix (size model.nv, ncols)");

      bp::def(
        "randomConfiguration", &randomConfiguration_proxy, bp::arg("model"),
        "Generate a random configuration in the bounds given by the lower and upper limits "
        "contained in model.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n");

      bp::def(
        "randomConfiguration",
        &randomConfiguration<Scalar, Options, JointCollectionDefaultTpl, VectorXs, VectorXs>,
        bp::args("model", "lower_bound", "upper_bound"),
        "Generate a random configuration in the bounds given by the Joint lower and upper limits "
        "arguments.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tlower_bound: the lower bound on the joint configuration vectors (size model.nq)\n"
        "\tupper_bound: the upper bound on the joint configuration vectors (size model.nq)\n");

      bp::def(
        "neutral", &neutral<Scalar, Options, JointCollectionDefaultTpl>, bp::arg("model"),
        "Returns the neutral configuration vector associated to the model.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n");

      bp::def(
        "normalize", normalize_proxy, bp::args("model", "q"),
        "Returns the configuration normalized.\n"
        "For instance, when the configuration vectors contains some quaternion values, it must be "
        "required to renormalize these components to keep orthonormal rotation values.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq: a joint configuration vector to normalize (size model.nq)\n");

      bp::def(
        "lieGroup", lieGroup_proxy, bp::args("model"),
        "Returns the Lie group associated to the model. It is the cartesian product of the lie "
        "groups of all its joints.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n");

      bp::def(
        "indexvInfo", indexvInfo_proxy, bp::args("model"),
        "Returns two vectors of size model.nq that gives for each q_i, the associated idx_v and nv "
        "of the joint for which q_i is a configuration component.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n");

#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS

      static const Scalar dummy_precision = Eigen::NumTraits<Scalar>::dummy_precision();
      bp::def(
        "isSameConfiguration",
        &isSameConfiguration<Scalar, Options, JointCollectionDefaultTpl, VectorXs, VectorXs>,
        bp::args("model", "q1", "q2", "prec"),
        "Return true if two configurations are equivalent within the given precision provided by "
        "prec.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq1: a joint configuration vector (size model.nq)\n"
        "\tq2: a joint configuration vector (size model.nq)\n"
        "\tprec: requested accuracy for the comparison\n");

      bp::def(
        "isNormalized", &isNormalized<Scalar, Options, JointCollectionDefaultTpl, VectorXs>,
        (bp::arg("model"), bp::arg("q"), bp::arg("prec") = dummy_precision),
        "Check whether a configuration vector is normalized within the given precision provided by "
        "prec.\n\n"
        "Parameters:\n"
        "\tmodel: model of the kinematic tree\n"
        "\tq: a joint configuration vector (size model.nq)\n"
        "\tprec: requested accuracy for the check\n");
#endif
    }

  } // namespace python
} // namespace pinocchio
