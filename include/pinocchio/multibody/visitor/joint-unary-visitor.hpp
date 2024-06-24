//
// Copyright (c) 2015-2024 CNRS INRIA
//

#ifndef __pinocchio_multibody_visitor_joint_unary_visitor_hpp__
#define __pinocchio_multibody_visitor_joint_unary_visitor_hpp__

#include <type_traits>
#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/get.hpp>

#include "pinocchio/multibody/visitor/fusion.hpp"
#include "pinocchio/multibody/joint/joint-base.hpp"

namespace pinocchio
{
  namespace fusion
  {

    namespace helper
    {
      /// \brief Add constness to T2 if T1 is const
      template<typename T1, typename T2>
      struct add_const_if_const
      {
        typedef T2 type;
      };

      template<typename T1, typename T2>
      struct add_const_if_const<const T1, T2>
      {
        typedef const T2 type;
      };
    } // namespace helper

    ///
    /// \brief Base structure for \b Unary visitation of a JointModel.
    ///        This structure provides runners to call the right visitor according to the number of
    ///        arguments. This should be used when deriving new rigid body algorithms.
    ///
    template<typename JointVisitorDerived, typename ReturnType = void>
    struct JointUnaryVisitorBase
    {

      template<
        typename Scalar,
        int Options,
        template<typename, int>
        class JointCollectionTpl,
        typename ArgsTmp>
      static ReturnType run(
        const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel,
        JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata,
        ArgsTmp args)
      {
        typedef JointModelTpl<Scalar, Options, JointCollectionTpl> JointModel;
        typedef JointDataTpl<Scalar, Options, JointCollectionTpl> JointData;

        InternalVisitorModelAndData<JointModel, JointData, ArgsTmp> visitor(jdata, args);
        return boost::apply_visitor(visitor, jmodel);
      }

      template<
        typename Scalar,
        int Options,
        template<typename, int>
        class JointCollectionTpl,
        typename ArgsTmp>
      static ReturnType run(
        const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel,
        const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata,
        ArgsTmp args)
      {
        typedef JointModelTpl<Scalar, Options, JointCollectionTpl> JointModel;
        typedef JointDataTpl<Scalar, Options, JointCollectionTpl> JointData;

        InternalVisitorModelAndData<JointModel, const JointData, ArgsTmp> visitor(jdata, args);
        return boost::apply_visitor(visitor, jmodel);
      }

      template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
      static ReturnType run(
        const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel,
        JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata)
      {
        typedef JointModelTpl<Scalar, Options, JointCollectionTpl> JointModel;
        typedef JointDataTpl<Scalar, Options, JointCollectionTpl> JointData;

        InternalVisitorModelAndData<JointModel, JointData, NoArg> visitor(jdata);
        return boost::apply_visitor(visitor, jmodel);
      }

      template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
      static ReturnType run(
        const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel,
        const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata)
      {
        typedef JointModelTpl<Scalar, Options, JointCollectionTpl> JointModel;
        typedef JointDataTpl<Scalar, Options, JointCollectionTpl> JointData;

        InternalVisitorModelAndData<JointModel, const JointData, NoArg> visitor(jdata);
        return boost::apply_visitor(visitor, jmodel);
      }

      template<typename JointModelDerived, typename ArgsTmp>
      static ReturnType run(
        const JointModelBase<JointModelDerived> & jmodel,
        typename JointModelBase<JointModelDerived>::JointDataDerived & jdata,
        ArgsTmp args)
      {
        typedef JointModelBase<JointModelDerived> JointModel;
        typedef typename JointModel::JointDataDerived JointData;

        InternalVisitorModelAndData<JointModelDerived, JointData, ArgsTmp> visitor(jdata, args);
        return visitor(jmodel.derived());
      }

      template<typename JointModelDerived, typename ArgsTmp>
      static ReturnType run(
        const JointModelBase<JointModelDerived> & jmodel,
        const typename JointModelBase<JointModelDerived>::JointDataDerived & jdata,
        ArgsTmp args)
      {
        typedef JointModelBase<JointModelDerived> JointModel;
        typedef typename JointModel::JointDataDerived JointData;

        InternalVisitorModelAndData<JointModelDerived, const JointData, ArgsTmp> visitor(
          jdata, args);
        return visitor(jmodel.derived());
      }

      template<typename JointModelDerived>
      static ReturnType run(
        const JointModelBase<JointModelDerived> & jmodel,
        typename JointModelBase<JointModelDerived>::JointDataDerived & jdata)
      {
        typedef JointModelBase<JointModelDerived> JointModel;
        typedef typename JointModel::JointDataDerived JointData;

        InternalVisitorModelAndData<JointModelDerived, JointData, NoArg> visitor(jdata);
        return visitor(jmodel.derived());
      }

      template<typename JointModelDerived>
      static ReturnType run(
        const JointModelBase<JointModelDerived> & jmodel,
        const typename JointModelBase<JointModelDerived>::JointDataDerived & jdata)
      {
        typedef JointModelBase<JointModelDerived> JointModel;
        typedef typename JointModel::JointDataDerived JointData;

        InternalVisitorModelAndData<JointModelDerived, const JointData, NoArg> visitor(jdata);
        return visitor(jmodel.derived());
      }

      template<
        typename Scalar,
        int Options,
        template<typename, int>
        class JointCollectionTpl,
        typename ArgsTmp>
      static ReturnType
      run(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel, ArgsTmp args)
      {
        InternalVisitorModel<ArgsTmp> visitor(args);
        return boost::apply_visitor(visitor, jmodel);
      }

      template<
        typename Scalar,
        int Options,
        template<typename, int>
        class JointCollectionTpl,
        typename ArgsTmp>
      static ReturnType
      run(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata, ArgsTmp args)
      {
        InternalVisitorModel<ArgsTmp> visitor(args);
        return boost::apply_visitor(visitor, jdata);
      }

      template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
      static ReturnType run(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel)
      {
        InternalVisitorModel<NoArg> visitor;
        return boost::apply_visitor(visitor, jmodel);
      }

      template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
      static ReturnType run(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata)
      {
        InternalVisitorModel<NoArg> visitor;
        return boost::apply_visitor(visitor, jdata);
      }

      template<typename JointModelDerived, typename ArgsTmp>
      static ReturnType run(const JointModelBase<JointModelDerived> & jmodel, ArgsTmp args)
      {
        InternalVisitorModel<ArgsTmp> visitor(args);
        return visitor(jmodel.derived());
      }

      template<typename JointDataDerived, typename ArgsTmp>
      static ReturnType run(const JointDataBase<JointDataDerived> & jdata, ArgsTmp args)
      {
        InternalVisitorModel<ArgsTmp> visitor(args);
        return visitor(jdata.derived());
      }

      template<typename JointModelDerived>
      static ReturnType run(const JointModelBase<JointModelDerived> & jmodel)
      {
        InternalVisitorModel<NoArg> visitor;
        return visitor(jmodel.derived());
      }

      template<typename JointDataDerived>
      static ReturnType run(const JointDataBase<JointDataDerived> & jdata)
      {
        InternalVisitorModel<NoArg> visitor;
        return visitor(jdata.derived());
      }

    private:
      template<typename JointModel, typename JointData, typename ArgType>
      struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
      {
        InternalVisitorModelAndData(JointData & jdata, ArgType args)
        : jdata(jdata)
        , args(args)
        {
        }

        template<typename JointModelDerived>
        ReturnType operator()(const JointModelBase<JointModelDerived> & jmodel) const
        {
          typedef typename helper::add_const_if_const<
            JointData, typename JointModelBase<JointModelDerived>::JointDataDerived>::type
            JointDataDerived;
          return bf::invoke(
            &JointVisitorDerived::template algo<JointModelDerived>,
            bf::append(
              boost::ref(jmodel.derived()), boost::ref(boost::get<JointDataDerived>(jdata)), args));
        }

        ReturnType operator()(const JointModelVoid)
        {
          return;
        }

        JointData & jdata;
        ArgType args;
      };

      template<typename JointModel, typename JointData>
      struct InternalVisitorModelAndData<JointModel, JointData, NoArg>
      : public boost::static_visitor<ReturnType>
      {
        InternalVisitorModelAndData(JointData & jdata)
        : jdata(jdata)
        {
        }

        template<typename JointModelDerived>
        ReturnType operator()(const JointModelBase<JointModelDerived> & jmodel) const
        {
          typedef typename helper::add_const_if_const<
            JointData, typename JointModelBase<JointModelDerived>::JointDataDerived>::type
            JointDataDerived;
          return bf::invoke(
            &JointVisitorDerived::template algo<JointModelDerived>,
            bf::make_vector(
              boost::ref(jmodel.derived()), boost::ref(boost::get<JointDataDerived>(jdata))));
        }

        JointData & jdata;
      };

      template<typename ArgType, typename Dummy = void>
      struct InternalVisitorModel : public boost::static_visitor<ReturnType>
      {
        InternalVisitorModel(ArgType args)
        : args(args)
        {
        }

        template<typename JointModelDerived>
        ReturnType operator()(const JointModelBase<JointModelDerived> & jmodel) const
        {
          return bf::invoke(
            &JointVisitorDerived::template algo<JointModelDerived>,
            bf::append(boost::ref(jmodel.derived()), args));
        }

        template<typename JointDataDerived>
        ReturnType operator()(const JointDataBase<JointDataDerived> & jdata) const
        {
          return bf::invoke(
            &JointVisitorDerived::template algo<JointDataDerived>,
            bf::append(boost::ref(jdata.derived()), args));
        }

        ReturnType operator()(const JointModelVoid)
        {
          return;
        }

        ArgType args;
      };

      template<typename Dummy>
      struct InternalVisitorModel<NoArg, Dummy> : public boost::static_visitor<ReturnType>
      {
        InternalVisitorModel()
        {
        }

        template<typename JointModelDerived>
        ReturnType operator()(const JointModelBase<JointModelDerived> & jmodel) const
        {
          return JointVisitorDerived::template algo<JointModelDerived>(jmodel.derived());
        }

        template<typename JointDataDerived>
        ReturnType operator()(const JointDataBase<JointDataDerived> & jdata) const
        {
          return JointVisitorDerived::template algo<JointDataDerived>(jdata.derived());
        }
      };
    }; // struct JointUnaryVisitorBase

    ///
    /// \brief This helper structure is now deprecated and has been replaced by
    /// JointUnaryVisitorBase.
    ///
    template<typename JointVisitorDerived, typename ReturnType = void>
    struct PINOCCHIO_DEPRECATED JointVisitorBase
    : JointUnaryVisitorBase<JointVisitorDerived, ReturnType>
    {
      typedef JointUnaryVisitorBase<JointVisitorDerived, ReturnType> Base;
      using Base::run;
    };

  } // namespace fusion
} // namespace pinocchio

#endif // ifndef __pinocchio_multibody_visitor_joint_unary_visitor_hpp__
