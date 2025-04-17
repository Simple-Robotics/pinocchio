//
// Copyright (c) 2018-2020 CNRS INRIA
//

#ifndef __pinocchio_multibody_liegroup_liegroup_algo_hxx__
#define __pinocchio_multibody_liegroup_liegroup_algo_hxx__

#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/multibody/joint/joint-composite.hpp"
#include "pinocchio/multibody/joint/joint-mimic.hpp"

namespace pinocchio
{

  namespace details
  {
    template<typename Visitor>
    struct Dispatch
    {
      template<
        typename Scalar,
        int Options,
        template<typename, int> class JointCollectionTpl,
        typename ArgsType>
      static void
      run(const JointModelCompositeTpl<Scalar, Options, JointCollectionTpl> & jmodel, ArgsType args)
      {
        for (size_t i = 0; i < jmodel.joints.size(); ++i)
          Visitor::run(jmodel.joints[i], args);
      }
    };

#define PINOCCHIO_DETAILS_WRITE_ARGS_0(JM) const JointModelBase<JM> & jmodel
#define PINOCCHIO_DETAILS_WRITE_ARGS_1(JM)                                                         \
  PINOCCHIO_DETAILS_WRITE_ARGS_0(JM), typename boost::fusion::result_of::at_c<ArgsType, 0>::type a0
#define PINOCCHIO_DETAILS_WRITE_ARGS_2(JM)                                                         \
  PINOCCHIO_DETAILS_WRITE_ARGS_1(JM), typename boost::fusion::result_of::at_c<ArgsType, 1>::type a1
#define PINOCCHIO_DETAILS_WRITE_ARGS_3(JM)                                                         \
  PINOCCHIO_DETAILS_WRITE_ARGS_2(JM), typename boost::fusion::result_of::at_c<ArgsType, 2>::type a2
#define PINOCCHIO_DETAILS_WRITE_ARGS_4(JM)                                                         \
  PINOCCHIO_DETAILS_WRITE_ARGS_3(JM), typename boost::fusion::result_of::at_c<ArgsType, 3>::type a3
#define PINOCCHIO_DETAILS_WRITE_ARGS_5(JM)                                                         \
  PINOCCHIO_DETAILS_WRITE_ARGS_4(JM), typename boost::fusion::result_of::at_c<ArgsType, 4>::type a4

#define PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_1(Algo)                                         \
  template<                                                                                        \
    typename Visitor, typename _Scalar, int _Options,                                              \
    template<typename, int> class JointCollection>                                                 \
  struct Algo<Visitor, JointModelCompositeTpl<_Scalar, _Options, JointCollection>>                 \
  {                                                                                                \
    typedef typename Visitor::ArgsType ArgsType;                                                   \
    typedef JointModelCompositeTpl<_Scalar, _Options, JointCollection> JointModelComposite;        \
    static void run(PINOCCHIO_DETAILS_WRITE_ARGS_1(JointModelComposite))                           \
    {                                                                                              \
      ::pinocchio::details::Dispatch<Visitor>::run(jmodel.derived(), ArgsType(a0));                \
    }                                                                                              \
  }

#define PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_2(Algo)                                         \
  template<                                                                                        \
    typename Visitor, typename _Scalar, int _Options,                                              \
    template<typename, int> class JointCollection>                                                 \
  struct Algo<Visitor, JointModelCompositeTpl<_Scalar, _Options, JointCollection>>                 \
  {                                                                                                \
    typedef typename Visitor::ArgsType ArgsType;                                                   \
    typedef JointModelCompositeTpl<_Scalar, _Options, JointCollection> JointModelComposite;        \
    static void run(PINOCCHIO_DETAILS_WRITE_ARGS_2(JointModelComposite))                           \
    {                                                                                              \
      ::pinocchio::details::Dispatch<Visitor>::run(jmodel.derived(), ArgsType(a0, a1));            \
    }                                                                                              \
  }

#define PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_3(Algo)                                         \
  template<                                                                                        \
    typename Visitor, typename _Scalar, int _Options,                                              \
    template<typename, int> class JointCollection>                                                 \
  struct Algo<Visitor, JointModelCompositeTpl<_Scalar, _Options, JointCollection>>                 \
  {                                                                                                \
    typedef typename Visitor::ArgsType ArgsType;                                                   \
    typedef JointModelCompositeTpl<_Scalar, _Options, JointCollection> JointModelComposite;        \
    static void run(PINOCCHIO_DETAILS_WRITE_ARGS_3(JointModelComposite))                           \
    {                                                                                              \
      ::pinocchio::details::Dispatch<Visitor>::run(jmodel.derived(), ArgsType(a0, a1, a2));        \
    }                                                                                              \
  }

#define PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_4(Algo)                                         \
  template<                                                                                        \
    typename Visitor, typename _Scalar, int _Options,                                              \
    template<typename, int> class JointCollection>                                                 \
  struct Algo<Visitor, JointModelCompositeTpl<_Scalar, _Options, JointCollection>>                 \
  {                                                                                                \
    typedef typename Visitor::ArgsType ArgsType;                                                   \
    typedef JointModelCompositeTpl<_Scalar, _Options, JointCollection> JointModelComposite;        \
    static void run(PINOCCHIO_DETAILS_WRITE_ARGS_4(JointModelComposite))                           \
    {                                                                                              \
      ::pinocchio::details::Dispatch<Visitor>::run(jmodel.derived(), ArgsType(a0, a1, a2, a3));    \
    }                                                                                              \
  }

#define PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_5(Algo)                                         \
  template<                                                                                        \
    typename Visitor, typename _Scalar, int _Options,                                              \
    template<typename, int> class JointCollection>                                                 \
  struct Algo<Visitor, JointModelCompositeTpl<_Scalar, _Options, JointCollection>>                 \
  {                                                                                                \
    typedef typename Visitor::ArgsType ArgsType;                                                   \
    typedef JointModelCompositeTpl<_Scalar, _Options, JointCollection> JointModelComposite;        \
    static void run(PINOCCHIO_DETAILS_WRITE_ARGS_5(JointModelComposite))                           \
    {                                                                                              \
      ::pinocchio::details::Dispatch<Visitor>::run(                                                \
        jmodel.derived(), ArgsType(a0, a1, a2, a3, a4));                                           \
    }                                                                                              \
  }

#define PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_1(Algo)                                               \
  template<                                                                                        \
    typename Visitor, typename _Scalar, int _Options,                                              \
    template<typename, int> class JointCollection>                                                 \
  struct Algo<Visitor, JointModelMimicTpl<_Scalar, _Options, JointCollection>>                     \
  {                                                                                                \
    typedef typename Visitor::ArgsType ArgsType;                                                   \
    typedef JointModelMimicTpl<_Scalar, _Options, JointCollection> JointModelMimic;                \
    static void run(PINOCCHIO_DETAILS_WRITE_ARGS_1(JointModelMimic))                               \
    {                                                                                              \
      PINOCCHIO_UNUSED_VARIABLE(jmodel);                                                           \
      PINOCCHIO_UNUSED_VARIABLE(a0);                                                               \
    }                                                                                              \
  }

#define PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_2(Algo)                                               \
  template<                                                                                        \
    typename Visitor, typename _Scalar, int _Options,                                              \
    template<typename, int> class JointCollection>                                                 \
  struct Algo<Visitor, JointModelMimicTpl<_Scalar, _Options, JointCollection>>                     \
  {                                                                                                \
    typedef typename Visitor::ArgsType ArgsType;                                                   \
    typedef JointModelMimicTpl<_Scalar, _Options, JointCollection> JointModelMimic;                \
    static void run(PINOCCHIO_DETAILS_WRITE_ARGS_2(JointModelMimic))                               \
    {                                                                                              \
      PINOCCHIO_UNUSED_VARIABLE(jmodel);                                                           \
      PINOCCHIO_UNUSED_VARIABLE(a0);                                                               \
      PINOCCHIO_UNUSED_VARIABLE(a1);                                                               \
    }                                                                                              \
  }

#define PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_3(Algo)                                               \
  template<                                                                                        \
    typename Visitor, typename _Scalar, int _Options,                                              \
    template<typename, int> class JointCollection>                                                 \
  struct Algo<Visitor, JointModelMimicTpl<_Scalar, _Options, JointCollection>>                     \
  {                                                                                                \
    typedef typename Visitor::ArgsType ArgsType;                                                   \
    typedef JointModelMimicTpl<_Scalar, _Options, JointCollection> JointModelMimic;                \
    static void run(PINOCCHIO_DETAILS_WRITE_ARGS_3(JointModelMimic))                               \
    {                                                                                              \
      PINOCCHIO_UNUSED_VARIABLE(jmodel);                                                           \
      PINOCCHIO_UNUSED_VARIABLE(a0);                                                               \
      PINOCCHIO_UNUSED_VARIABLE(a1);                                                               \
      PINOCCHIO_UNUSED_VARIABLE(a2);                                                               \
    }                                                                                              \
  }

#define PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_4(Algo)                                               \
  template<                                                                                        \
    typename Visitor, typename _Scalar, int _Options,                                              \
    template<typename, int> class JointCollection>                                                 \
  struct Algo<Visitor, JointModelMimicTpl<_Scalar, _Options, JointCollection>>                     \
  {                                                                                                \
    typedef typename Visitor::ArgsType ArgsType;                                                   \
    typedef JointModelMimicTpl<_Scalar, _Options, JointCollection> JointModelMimic;                \
    static void run(PINOCCHIO_DETAILS_WRITE_ARGS_4(JointModelMimic))                               \
    {                                                                                              \
      PINOCCHIO_UNUSED_VARIABLE(jmodel);                                                           \
      PINOCCHIO_UNUSED_VARIABLE(a0);                                                               \
      PINOCCHIO_UNUSED_VARIABLE(a1);                                                               \
      PINOCCHIO_UNUSED_VARIABLE(a2);                                                               \
      PINOCCHIO_UNUSED_VARIABLE(a3);                                                               \
    }                                                                                              \
  }

#define PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_5(Algo)                                               \
  template<                                                                                        \
    typename Visitor, typename _Scalar, int _Options,                                              \
    template<typename, int> class JointCollection>                                                 \
  struct Algo<Visitor, JointModelMimicTpl<_Scalar, _Options, JointCollection>>                     \
  {                                                                                                \
    typedef typename Visitor::ArgsType ArgsType;                                                   \
    typedef JointModelMimicTpl<_Scalar, _Options, JointCollection> JointModelMimic;                \
    static void run(PINOCCHIO_DETAILS_WRITE_ARGS_5(JointModelMimic))                               \
    {                                                                                              \
      PINOCCHIO_UNUSED_VARIABLE(jmodel);                                                           \
      PINOCCHIO_UNUSED_VARIABLE(a0);                                                               \
      PINOCCHIO_UNUSED_VARIABLE(a1);                                                               \
      PINOCCHIO_UNUSED_VARIABLE(a2);                                                               \
      PINOCCHIO_UNUSED_VARIABLE(a3);                                                               \
      PINOCCHIO_UNUSED_VARIABLE(a4);                                                               \
    }                                                                                              \
  }

#define PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_1(Algo, Visitor)                                     \
  typedef LieGroup_t LieGroupMap;                                                                  \
  template<typename JointModel>                                                                    \
  static void algo(PINOCCHIO_DETAILS_WRITE_ARGS_1(JointModel))                                     \
  {                                                                                                \
    AlgoDispatch<JointModel>::run(jmodel, a0);                                                     \
  }                                                                                                \
  template<typename JointModel>                                                                    \
  struct AlgoDispatch : Algo<Visitor, JointModel>                                                  \
  {                                                                                                \
    using Algo<Visitor, JointModel>::run;                                                          \
  };

#define PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_2(Algo, Visitor)                                     \
  typedef LieGroup_t LieGroupMap;                                                                  \
  template<typename JointModel>                                                                    \
  static void algo(PINOCCHIO_DETAILS_WRITE_ARGS_2(JointModel))                                     \
  {                                                                                                \
    AlgoDispatch<JointModel>::run(jmodel, a0, a1);                                                 \
  }                                                                                                \
  template<typename JointModel>                                                                    \
  struct AlgoDispatch : Algo<Visitor, JointModel>                                                  \
  {                                                                                                \
    using Algo<Visitor, JointModel>::run;                                                          \
  };

#define PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_3(Algo, Visitor)                                     \
  typedef LieGroup_t LieGroupMap;                                                                  \
  template<typename JointModel>                                                                    \
  static void algo(PINOCCHIO_DETAILS_WRITE_ARGS_3(JointModel))                                     \
  {                                                                                                \
    AlgoDispatch<JointModel>::run(jmodel, a0, a1, a2);                                             \
  }                                                                                                \
  template<typename JointModel>                                                                    \
  struct AlgoDispatch : Algo<Visitor, JointModel>                                                  \
  {                                                                                                \
    using Algo<Visitor, JointModel>::run;                                                          \
  };

#define PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_4(Algo, Visitor)                                     \
  typedef LieGroup_t LieGroupMap;                                                                  \
  template<typename JointModel>                                                                    \
  static void algo(PINOCCHIO_DETAILS_WRITE_ARGS_4(JointModel))                                     \
  {                                                                                                \
    AlgoDispatch<JointModel>::run(jmodel, a0, a1, a2, a3);                                         \
  }                                                                                                \
  template<typename JointModel>                                                                    \
  struct AlgoDispatch : Algo<Visitor, JointModel>                                                  \
  {                                                                                                \
    using Algo<Visitor, JointModel>::run;                                                          \
  };

#define PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_5(Algo, Visitor)                                     \
  typedef LieGroup_t LieGroupMap;                                                                  \
  template<typename JointModel>                                                                    \
  static void algo(PINOCCHIO_DETAILS_WRITE_ARGS_5(JointModel))                                     \
  {                                                                                                \
    AlgoDispatch<JointModel>::run(jmodel, a0, a1, a2, a3, a4);                                     \
  }                                                                                                \
  template<typename JointModel>                                                                    \
  struct AlgoDispatch : Algo<Visitor, JointModel>                                                  \
  {                                                                                                \
    using Algo<Visitor, JointModel>::run;                                                          \
  };

  } // namespace details

  template<typename Visitor, typename JointModel>
  struct IntegrateStepAlgo;

  template<
    typename LieGroup_t,
    typename ConfigVectorIn,
    typename TangentVectorIn,
    typename ConfigVectorOut>
  struct IntegrateStep
  : public fusion::JointUnaryVisitorBase<
      IntegrateStep<LieGroup_t, ConfigVectorIn, TangentVectorIn, ConfigVectorOut>>
  {
    typedef boost::fusion::
      vector<const ConfigVectorIn &, const TangentVectorIn &, ConfigVectorOut &>
        ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_3(IntegrateStepAlgo, IntegrateStep)
  };

  template<typename Visitor, typename JointModel>
  struct IntegrateStepAlgo
  {
    template<typename ConfigVectorIn, typename TangentVector, typename ConfigVectorOut>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn> & q,
      const Eigen::MatrixBase<TangentVector> & v,
      const Eigen::MatrixBase<ConfigVectorOut> & result)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.integrate(
        jmodel.jointConfigSelector(q.derived()), jmodel.jointVelocitySelector(v.derived()),
        jmodel.jointConfigSelector(PINOCCHIO_EIGEN_CONST_CAST(ConfigVectorOut, result)));
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_3(IntegrateStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_3(IntegrateStepAlgo);

  template<typename Visitor, typename JointModel>
  struct dIntegrateStepAlgo;

  template<
    typename LieGroup_t,
    typename ConfigVectorIn,
    typename TangentVectorIn,
    typename JacobianMatrixType>
  struct dIntegrateStep
  : public fusion::JointUnaryVisitorBase<
      dIntegrateStep<LieGroup_t, ConfigVectorIn, TangentVectorIn, JacobianMatrixType>>
  {
    typedef boost::fusion::vector<
      const ConfigVectorIn &,
      const TangentVectorIn &,
      JacobianMatrixType &,
      const ArgumentPosition &,
      const AssignmentOperatorType &>
      ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_5(dIntegrateStepAlgo, dIntegrateStep)
  };

  template<typename Visitor, typename JointModel>
  struct dIntegrateStepAlgo
  {
    template<typename ConfigVectorIn, typename TangentVector, typename JacobianMatrixType>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn> & q,
      const Eigen::MatrixBase<TangentVector> & v,
      const Eigen::MatrixBase<JacobianMatrixType> & mat,
      const ArgumentPosition & arg,
      const AssignmentOperatorType & op)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.dIntegrate(
        jmodel.jointConfigSelector(q.derived()), jmodel.jointVelocitySelector(v.derived()),
        jmodel.jointBlock(PINOCCHIO_EIGEN_CONST_CAST(JacobianMatrixType, mat)), arg, op);
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_5(dIntegrateStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_5(dIntegrateStepAlgo);

  template<typename Visitor, typename JointModel>
  struct TangentMapStepAlgo;

  template<
    typename LieGroup_t,
    typename IsCompact_t,
    typename ConfigVectorIn,
    typename TangentMapMatrixType>
  struct TangentMapStep
  : public fusion::JointUnaryVisitorBase<
      TangentMapStep<LieGroup_t, IsCompact_t, ConfigVectorIn, TangentMapMatrixType>>
  {
    typedef IsCompact_t IsCompact;
    typedef boost::fusion::
      vector<const ConfigVectorIn &, TangentMapMatrixType &, const AssignmentOperatorType &>
        ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_3(TangentMapStepAlgo, TangentMapStep)
  };

  template<typename Visitor, typename JointModel>
  struct TangentMapStepAlgo
  {
    template<typename ConfigVectorIn, typename TangentMapMatrixType>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn> & q,
      const Eigen::MatrixBase<TangentMapMatrixType> & TM,
      const AssignmentOperatorType & op)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.tangentMap(
        jmodel.jointConfigSelector(q.derived()),
        jmodel.jointQVBlock(
          PINOCCHIO_EIGEN_CONST_CAST(TangentMapMatrixType, TM), jmodel.idx_q(),
          Visitor::IsCompact::value ? 0 : jmodel.idx_v()),
        Visitor::IsCompact::value ? SETTO : op);
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_3(TangentMapStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_3(TangentMapStepAlgo);

  template<typename Visitor, typename JointModel>
  struct CompactSetTangentMapStepAlgo;

  template<typename LieGroup_t, typename ConfigVectorIn, typename CompactSetTangentMapMatrixType>
  struct CompactSetTangentMapStep
  : public fusion::JointUnaryVisitorBase<
      CompactSetTangentMapStep<LieGroup_t, ConfigVectorIn, CompactSetTangentMapMatrixType>>
  {
    typedef boost::fusion::vector<const ConfigVectorIn &, CompactSetTangentMapMatrixType &>
      ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_2(CompactSetTangentMapStepAlgo, CompactSetTangentMapStep)
  };

  template<typename Visitor, typename JointModel>
  struct CompactSetTangentMapStepAlgo
  {
    template<typename ConfigVectorIn, typename CompactSetTangentMapMatrixType>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn> & q,
      const Eigen::MatrixBase<CompactSetTangentMapMatrixType> & TMc)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.tangentMap(
        jmodel.jointConfigSelector(q.derived()),
        jmodel.jointQVBlock(
          PINOCCHIO_EIGEN_CONST_CAST(CompactSetTangentMapMatrixType, TMc), jmodel.idx_q(), 0),
        SETTO);
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_2(CompactSetTangentMapStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_2(CompactSetTangentMapStepAlgo);

  template<typename Visitor, typename JointModel>
  struct TangentMapProductStepAlgo;

  template<
    typename LieGroup_t,
    typename ConfigVectorIn,
    typename MatrixInType,
    typename MatrixOutType>
  struct TangentMapProductStep
  : public fusion::JointUnaryVisitorBase<
      TangentMapProductStep<LieGroup_t, ConfigVectorIn, MatrixInType, MatrixOutType>>
  {
    typedef boost::fusion::vector<
      const ConfigVectorIn &,
      const MatrixInType &,
      MatrixOutType &,
      const AssignmentOperatorType &>
      ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_4(TangentMapProductStepAlgo, TangentMapProductStep)
  };

  template<typename Visitor, typename JointModel>
  struct TangentMapProductStepAlgo
  {
    template<typename ConfigVectorIn, typename MatrixInType, typename MatrixOutType>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn> & q,
      const Eigen::MatrixBase<MatrixInType> & mat_in,
      const Eigen::MatrixBase<MatrixOutType> & mat_out,
      const AssignmentOperatorType & op)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.tangentMapProduct(
        jmodel.jointConfigSelector(q.derived()), jmodel.jointRows(mat_in.derived()),
        jmodel.jointQRows(PINOCCHIO_EIGEN_CONST_CAST(MatrixOutType, mat_out)), op);
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_4(TangentMapProductStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_4(TangentMapProductStepAlgo);

  template<typename Visitor, typename JointModel>
  struct CoTangentMapProductStepAlgo;

  template<
    typename LieGroup_t,
    typename ConfigVectorIn,
    typename MatrixInType,
    typename MatrixOutType>
  struct CoTangentMapProductStep
  : public fusion::JointUnaryVisitorBase<
      CoTangentMapProductStep<LieGroup_t, ConfigVectorIn, MatrixInType, MatrixOutType>>
  {
    typedef boost::fusion::vector<
      const ConfigVectorIn &,
      const MatrixInType &,
      MatrixOutType &,
      const AssignmentOperatorType &>
      ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_4(CoTangentMapProductStepAlgo, CoTangentMapProductStep)
  };

  template<typename Visitor, typename JointModel>
  struct CoTangentMapProductStepAlgo
  {
    template<typename ConfigVectorIn, typename MatrixInType, typename MatrixOutType>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn> & q,
      const Eigen::MatrixBase<MatrixInType> & mat_in,
      const Eigen::MatrixBase<MatrixOutType> & mat_out,
      const AssignmentOperatorType & op)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.coTangentMapProduct(
        jmodel.jointConfigSelector(q.derived()), jmodel.jointQRows(mat_in.derived()),
        jmodel.jointRows(PINOCCHIO_EIGEN_CONST_CAST(MatrixOutType, mat_out)), op);
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_4(CoTangentMapProductStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_4(CoTangentMapProductStepAlgo);

  template<typename Visitor, typename JointModel>
  struct dIntegrateTransportStepAlgo;

  template<
    typename LieGroup_t,
    typename ConfigVectorIn,
    typename TangentVectorIn,
    typename JacobianMatrixInType,
    typename JacobianMatrixOutType>
  struct dIntegrateTransportStep
  : public fusion::JointUnaryVisitorBase<dIntegrateTransportStep<
      LieGroup_t,
      ConfigVectorIn,
      TangentVectorIn,
      JacobianMatrixInType,
      JacobianMatrixOutType>>
  {
    typedef boost::fusion::vector<
      const ConfigVectorIn &,
      const TangentVectorIn &,
      const JacobianMatrixInType &,
      JacobianMatrixOutType &,
      const ArgumentPosition &>
      ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_5(dIntegrateTransportStepAlgo, dIntegrateTransportStep)
  };

  template<typename Visitor, typename JointModel>
  struct dIntegrateTransportStepAlgo
  {
    template<
      typename ConfigVectorIn,
      typename TangentVector,
      typename JacobianMatrixInType,
      typename JacobianMatrixOutType>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn> & q,
      const Eigen::MatrixBase<TangentVector> & v,
      const Eigen::MatrixBase<JacobianMatrixInType> & mat_in,
      const Eigen::MatrixBase<JacobianMatrixOutType> & mat_out,
      const ArgumentPosition & arg)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.dIntegrateTransport(
        jmodel.jointConfigSelector(q.derived()), jmodel.jointVelocitySelector(v.derived()),
        jmodel.jointRows(mat_in.derived()),
        jmodel.jointRows(PINOCCHIO_EIGEN_CONST_CAST(JacobianMatrixOutType, mat_out)), arg);
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_5(dIntegrateTransportStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_5(dIntegrateTransportStepAlgo);

  template<typename Visitor, typename JointModel>
  struct dIntegrateTransportInPlaceStepAlgo;

  template<
    typename LieGroup_t,
    typename ConfigVectorIn,
    typename TangentVectorIn,
    typename JacobianMatrixType>
  struct dIntegrateTransportInPlaceStep
  : public fusion::JointUnaryVisitorBase<dIntegrateTransportInPlaceStep<
      LieGroup_t,
      ConfigVectorIn,
      TangentVectorIn,
      JacobianMatrixType>>
  {
    typedef boost::fusion::vector<
      const ConfigVectorIn &,
      const TangentVectorIn &,
      JacobianMatrixType &,
      const ArgumentPosition &>
      ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_4(
      dIntegrateTransportInPlaceStepAlgo, dIntegrateTransportInPlaceStep)
  };

  template<typename Visitor, typename JointModel>
  struct dIntegrateTransportInPlaceStepAlgo
  {
    template<typename ConfigVectorIn, typename TangentVector, typename JacobianMatrixType>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn> & q,
      const Eigen::MatrixBase<TangentVector> & v,
      const Eigen::MatrixBase<JacobianMatrixType> & mat,
      const ArgumentPosition & arg)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.dIntegrateTransport(
        jmodel.jointConfigSelector(q.derived()), jmodel.jointVelocitySelector(v.derived()),
        jmodel.jointRows(PINOCCHIO_EIGEN_CONST_CAST(JacobianMatrixType, mat)), arg);
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_4(dIntegrateTransportInPlaceStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_4(dIntegrateTransportInPlaceStepAlgo);

  template<typename Visitor, typename JointModel>
  struct dDifferenceStepAlgo;

  template<
    typename LieGroup_t,
    typename ConfigVector1,
    typename ConfigVector2,
    typename JacobianMatrix>
  struct dDifferenceStep
  : public fusion::JointUnaryVisitorBase<
      dDifferenceStep<LieGroup_t, ConfigVector1, ConfigVector2, JacobianMatrix>>
  {
    typedef boost::fusion::vector<
      const ConfigVector1 &,
      const ConfigVector2 &,
      JacobianMatrix &,
      const ArgumentPosition &>
      ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_4(dDifferenceStepAlgo, dDifferenceStep)
  };

  template<typename Visitor, typename JointModel>
  struct dDifferenceStepAlgo
  {
    template<typename ConfigVector1, typename ConfigVector2, typename JacobianMatrix>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVector1> & q0,
      const Eigen::MatrixBase<ConfigVector2> & q1,
      const Eigen::MatrixBase<JacobianMatrix> & mat,
      const ArgumentPosition & arg)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.dDifference(
        jmodel.jointConfigSelector(q0.derived()), jmodel.jointConfigSelector(q1.derived()),
        jmodel.jointBlock(PINOCCHIO_EIGEN_CONST_CAST(JacobianMatrix, mat)), arg);
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_4(dDifferenceStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_4(dDifferenceStepAlgo);

  template<typename Visitor, typename JointModel>
  struct InterpolateStepAlgo;

  template<
    typename LieGroup_t,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2,
    typename Scalar,
    typename ConfigVectorOut>
  struct InterpolateStep
  : public fusion::JointUnaryVisitorBase<
      InterpolateStep<LieGroup_t, ConfigVectorIn1, ConfigVectorIn2, Scalar, ConfigVectorOut>>
  {
    typedef boost::fusion::
      vector<const ConfigVectorIn1 &, const ConfigVectorIn2 &, const Scalar &, ConfigVectorOut &>
        ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_4(InterpolateStepAlgo, InterpolateStep)
  };

  template<typename Visitor, typename JointModel>
  struct InterpolateStepAlgo
  {
    template<
      typename ConfigVectorIn1,
      typename ConfigVectorIn2,
      typename Scalar,
      typename ConfigVectorOut>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn1> & q0,
      const Eigen::MatrixBase<ConfigVectorIn2> & q1,
      const Scalar & u,
      const Eigen::MatrixBase<ConfigVectorOut> & result)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.interpolate(
        jmodel.jointConfigSelector(q0.derived()), jmodel.jointConfigSelector(q1.derived()), u,
        jmodel.jointConfigSelector(PINOCCHIO_EIGEN_CONST_CAST(ConfigVectorOut, result)));
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_4(InterpolateStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_4(InterpolateStepAlgo);

  template<typename Visitor, typename JointModel>
  struct DifferenceStepAlgo;

  template<
    typename LieGroup_t,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2,
    typename TangentVectorOut>
  struct DifferenceStep
  : public fusion::JointUnaryVisitorBase<
      DifferenceStep<LieGroup_t, ConfigVectorIn1, ConfigVectorIn2, TangentVectorOut>>
  {
    typedef boost::fusion::
      vector<const ConfigVectorIn1 &, const ConfigVectorIn2 &, TangentVectorOut &>
        ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_3(DifferenceStepAlgo, DifferenceStep)
  };

  template<typename Visitor, typename JointModel>
  struct DifferenceStepAlgo
  {
    template<typename ConfigVectorIn1, typename ConfigVectorIn2, typename TangentVectorOut>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn1> & q0,
      const Eigen::MatrixBase<ConfigVectorIn2> & q1,
      const Eigen::MatrixBase<TangentVectorOut> & result)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.difference(
        jmodel.jointConfigSelector(q0.derived()), jmodel.jointConfigSelector(q1.derived()),
        jmodel.jointVelocitySelector(PINOCCHIO_EIGEN_CONST_CAST(TangentVectorOut, result)));
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_3(DifferenceStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_3(DifferenceStepAlgo);

  template<typename Visitor, typename JointModel>
  struct SquaredDistanceStepAlgo;

  template<
    typename LieGroup_t,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2,
    typename DistanceVectorOut>
  struct SquaredDistanceStep
  : public fusion::JointUnaryVisitorBase<
      SquaredDistanceStep<LieGroup_t, ConfigVectorIn1, ConfigVectorIn2, DistanceVectorOut>>
  {
    typedef boost::fusion::vector<
      const JointIndex &,
      const ConfigVectorIn1 &,
      const ConfigVectorIn2 &,
      DistanceVectorOut &>
      ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_4(SquaredDistanceStepAlgo, SquaredDistanceStep)
  };

  template<typename Visitor, typename JointModel>
  struct SquaredDistanceStepAlgo
  {
    template<typename ConfigVectorIn1, typename ConfigVectorIn2, typename DistanceVectorOut>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const JointIndex & i,
      const Eigen::MatrixBase<ConfigVectorIn1> & q0,
      const Eigen::MatrixBase<ConfigVectorIn2> & q1,
      const Eigen::MatrixBase<DistanceVectorOut> & distances)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;
      typename LieGroupMap::template operation<JointModel>::type lgo;
      DistanceVectorOut & distances_ = PINOCCHIO_EIGEN_CONST_CAST(DistanceVectorOut, distances);
      distances_[(Eigen::DenseIndex)i] += lgo.squaredDistance(
        jmodel.jointConfigSelector(q0.derived()), jmodel.jointConfigSelector(q1.derived()));
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_4(SquaredDistanceStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_4(SquaredDistanceStepAlgo);

  template<typename Visitor, typename JointModel>
  struct SquaredDistanceSumStepAlgo;

  template<typename LieGroup_t, typename ConfigVectorIn1, typename ConfigVectorIn2, typename Scalar>
  struct SquaredDistanceSumStep
  : public fusion::JointUnaryVisitorBase<
      SquaredDistanceSumStep<LieGroup_t, ConfigVectorIn1, ConfigVectorIn2, Scalar>>
  {
    typedef boost::fusion::vector<const ConfigVectorIn1 &, const ConfigVectorIn2 &, Scalar &>
      ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_3(SquaredDistanceSumStepAlgo, SquaredDistanceSumStep)
  };

  template<typename Visitor, typename JointModel>
  struct SquaredDistanceSumStepAlgo
  {
    template<typename ConfigVectorIn1, typename ConfigVectorIn2>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn1> & q0,
      const Eigen::MatrixBase<ConfigVectorIn2> & q1,
      typename ConfigVectorIn1::Scalar & squaredDistance)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;
      typename LieGroupMap::template operation<JointModel>::type lgo;
      squaredDistance += lgo.squaredDistance(
        jmodel.jointConfigSelector(q0.derived()), jmodel.jointConfigSelector(q1.derived()));
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_3(SquaredDistanceSumStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_3(SquaredDistanceSumStepAlgo);

  template<typename Visitor, typename JointModel>
  struct RandomConfigurationStepAlgo;

  template<
    typename LieGroup_t,
    typename ConfigVectorOut,
    typename ConfigVectorIn1,
    typename ConfigVectorIn2>
  struct RandomConfigurationStep
  : public fusion::JointUnaryVisitorBase<
      RandomConfigurationStep<LieGroup_t, ConfigVectorOut, ConfigVectorIn1, ConfigVectorIn2>>
  {
    typedef boost::fusion::
      vector<ConfigVectorOut &, const ConfigVectorIn1 &, const ConfigVectorIn2 &>
        ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_3(RandomConfigurationStepAlgo, RandomConfigurationStep)
  };

  template<typename Visitor, typename JointModel>
  struct RandomConfigurationStepAlgo
  {
    template<typename ConfigVectorOut, typename ConfigVectorIn1, typename ConfigVectorIn2>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorOut> & q,
      const Eigen::MatrixBase<ConfigVectorIn1> & lowerLimits,
      const Eigen::MatrixBase<ConfigVectorIn2> & upperLimits)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.randomConfiguration(
        jmodel.jointConfigSelector(lowerLimits.derived()),
        jmodel.jointConfigSelector(upperLimits.derived()),
        jmodel.jointConfigSelector(PINOCCHIO_EIGEN_CONST_CAST(ConfigVectorOut, q)));
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_3(RandomConfigurationStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_3(RandomConfigurationStepAlgo);

  template<typename Visitor, typename JointModel>
  struct NormalizeStepAlgo;

  template<typename LieGroup_t, typename ConfigVectorType>
  struct NormalizeStep
  : public fusion::JointUnaryVisitorBase<NormalizeStep<LieGroup_t, ConfigVectorType>>
  {
    typedef boost::fusion::vector<ConfigVectorType &> ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_1(NormalizeStepAlgo, NormalizeStep)
  };

  template<typename Visitor, typename JointModel>
  struct NormalizeStepAlgo
  {
    template<typename ConfigVectorType>
    static void
    run(const JointModelBase<JointModel> & jmodel, const Eigen::MatrixBase<ConfigVectorType> & qout)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      lgo.normalize(jmodel.jointConfigSelector(PINOCCHIO_EIGEN_CONST_CAST(ConfigVectorType, qout)));
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_1(NormalizeStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_1(NormalizeStepAlgo);

  template<typename Visitor, typename JointModel>
  struct IsNormalizedStepAlgo;

  template<typename LieGroup_t, typename ConfigVectorIn, typename Scalar>
  struct IsNormalizedStep
  : public fusion::JointUnaryVisitorBase<IsNormalizedStep<LieGroup_t, ConfigVectorIn, Scalar>>
  {
    typedef boost::fusion::vector<const ConfigVectorIn &, const Scalar &, bool &> ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_3(IsNormalizedStepAlgo, IsNormalizedStep)
  };

  template<typename Visitor, typename JointModel>
  struct IsNormalizedStepAlgo
  {
    template<typename ConfigVectorIn>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorIn> & q,
      const typename ConfigVectorIn::Scalar & prec,
      bool & res)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      res &= lgo.isNormalized(jmodel.jointConfigSelector(q.derived()), prec);
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_3(IsNormalizedStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_3(IsNormalizedStepAlgo);

  template<typename Visitor, typename JointModel>
  struct IsSameConfigurationStepAlgo;

  template<typename LieGroup_t, typename ConfigVectorIn1, typename ConfigVectorIn2, typename Scalar>
  struct IsSameConfigurationStep
  : public fusion::JointUnaryVisitorBase<
      IsSameConfigurationStep<LieGroup_t, ConfigVectorIn1, ConfigVectorIn2, Scalar>>
  {
    typedef boost::fusion::
      vector<bool &, const ConfigVectorIn1 &, const ConfigVectorIn2 &, const Scalar &>
        ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_4(IsSameConfigurationStepAlgo, IsSameConfigurationStep)
  };

  template<typename Visitor, typename JointModel>
  struct IsSameConfigurationStepAlgo
  {
    template<typename ConfigVectorIn1, typename ConfigVectorIn2>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      bool & isSame,
      const Eigen::MatrixBase<ConfigVectorIn1> & q1,
      const Eigen::MatrixBase<ConfigVectorIn2> & q2,
      const typename ConfigVectorIn1::Scalar & prec)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      isSame &= lgo.isSameConfiguration(
        jmodel.jointConfigSelector(q1.derived()), jmodel.jointConfigSelector(q2.derived()), prec);
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_4(IsSameConfigurationStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_4(IsSameConfigurationStepAlgo);

  template<typename Visitor, typename JointModel>
  struct NeutralStepAlgo;

  template<typename LieGroup_t, typename ConfigVectorType>
  struct NeutralStep
  : public fusion::JointUnaryVisitorBase<NeutralStep<LieGroup_t, ConfigVectorType>>
  {
    typedef boost::fusion::vector<ConfigVectorType &> ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_1(NeutralStepAlgo, NeutralStep)
  };

  template<typename Visitor, typename JointModel>
  struct NeutralStepAlgo
  {
    template<typename ConfigVectorType>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorType> & neutral_elt)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typename LieGroupMap::template operation<JointModel>::type lgo;
      jmodel.jointConfigSelector(PINOCCHIO_EIGEN_CONST_CAST(ConfigVectorType, neutral_elt)) =
        lgo.neutral();
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_1(NeutralStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_1(NeutralStepAlgo);

  template<typename Visitor, typename JointModel>
  struct IntegrateCoeffWiseJacobianStepAlgo;

  template<typename LieGroup_t, typename ConfigVectorType, typename JacobianMatrix>
  struct IntegrateCoeffWiseJacobianStep
  : public fusion::JointUnaryVisitorBase<
      IntegrateCoeffWiseJacobianStep<LieGroup_t, ConfigVectorType, JacobianMatrix>>
  {
    typedef boost::fusion::vector<const ConfigVectorType &, JacobianMatrix &> ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_2(
      IntegrateCoeffWiseJacobianStepAlgo, IntegrateCoeffWiseJacobianStep)
  };

  template<typename Visitor, typename JointModel>
  struct IntegrateCoeffWiseJacobianStepAlgo
  {
    template<typename ConfigVectorType, typename JacobianMatrix>
    static void run(
      const JointModelBase<JointModel> & jmodel,
      const Eigen::MatrixBase<ConfigVectorType> & q,
      const Eigen::MatrixBase<JacobianMatrix> & jacobian)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;

      typedef typename LieGroupMap::template operation<JointModel>::type LieGroup;
      LieGroup lgo;
      lgo.integrateCoeffWiseJacobian(
        jmodel.jointConfigSelector(q.derived()),
        PINOCCHIO_EIGEN_CONST_CAST(JacobianMatrix, jacobian)
          .template block<LieGroup::NQ, LieGroup::NV>(
            jmodel.idx_q(), jmodel.idx_v(), jmodel.nq(), jmodel.nv()));
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_2(IntegrateCoeffWiseJacobianStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_2(IntegrateCoeffWiseJacobianStepAlgo);

  template<typename Visitor, typename JointModel>
  struct LieGroupInstanceStepAlgo;

  template<typename LieGroup_t, typename Scalar, int Options>
  struct LieGroupInstanceStep
  : public fusion::JointUnaryVisitorBase<LieGroupInstanceStep<LieGroup_t, Scalar, Options>>
  {
    typedef typename LieGroup_t::template product_variant<Scalar, Options>::type LgType;

    typedef boost::fusion::vector<LgType &> ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_1(LieGroupInstanceStepAlgo, LieGroupInstanceStep)
  };

  template<typename Visitor, typename JointModel>
  struct LieGroupInstanceStepAlgo
  {
    typedef typename Visitor::LgType LgType;

    static void run(const JointModelBase<JointModel> & jmodel, LgType & res_lgo)
    {
      typedef typename Visitor::LieGroupMap LieGroupMap;
      res_lgo *= jmodel.template lie_group<LieGroupMap>();
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_1(LieGroupInstanceStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_1(LieGroupInstanceStepAlgo);

  template<typename Visitor, typename JointModel>
  struct IndexvInfoStepAlgo;

  struct IndexvInfoStep : public fusion::JointUnaryVisitorBase<IndexvInfoStep>
  {
    typedef boost::blank LieGroup_t;
    typedef boost::fusion::vector<std::vector<int> &, std::vector<int> &> ArgsType;

    PINOCCHIO_DETAILS_VISITOR_METHOD_ALGO_2(IndexvInfoStepAlgo, IndexvInfoStep)
  };

  template<typename Visitor, typename JointModel>
  struct IndexvInfoStepAlgo
  {
    static void run(
      const JointModelBase<JointModel> & jmodel, std::vector<int> & nvs, std::vector<int> & idx_vs)
    {
      int idx_v = jmodel.idx_v();
      int nv = jmodel.nv();
      size_t idx_q = static_cast<size_t>(jmodel.idx_q());

      for (size_t idx = idx_q; idx < idx_q + static_cast<size_t>(jmodel.nq()); ++idx)
      {
        nvs[idx] = nv;
        idx_vs[idx] = idx_v;
      }
    }
  };

  PINOCCHIO_DETAILS_DISPATCH_JOINT_COMPOSITE_2(IndexvInfoStepAlgo);
  PINOCCHIO_DETAILS_CANCEL_JOINT_MIMIC_2(IndexvInfoStepAlgo);

} // namespace pinocchio

#endif // ifndef __pinocchio_multibody_liegroup_liegroup_algo_hxx__
