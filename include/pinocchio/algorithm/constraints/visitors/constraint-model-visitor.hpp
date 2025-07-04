//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_model_visitor_hpp__
#define __pinocchio_algorithm_constraints_constraint_model_visitor_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/multibody/visitor/fusion.hpp"

namespace pinocchio
{

  namespace visitors
  {

    namespace bf = boost::fusion;
    using fusion::NoArg;

    namespace internal
    {
      template<typename T>
      struct NoRun
      {
        static T run()
        {
          assert(false && "Should never happened.");
          // Hacky way to not have to return something real. The system should throw before.
          const typename std::remove_reference<T>::type * null_ptr = NULL;
          return *null_ptr;
        }
      };

      // Specialization for reference types
      template<typename T>
      struct NoRun<T &>
      {
        static T & run()
        {
          assert(false && "Should never happen.");
          // Hacky way to not have to return something real. The system should throw before.
          T * null_ptr = nullptr;
          return *null_ptr;
        }
      };

      template<>
      struct NoRun<void>
      {
        static void run()
        {
          return;
        }
      };
    } // namespace internal

    /**
     * @brief      ConstraintModelShortnameVisitor visitor
     */
    struct ConstraintModelShortnameVisitor : boost::static_visitor<std::string>
    {
      template<typename ConstraintModelDerived>
      std::string operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
      {
        return cmodel.shortname();
      }
      std::string operator()(const boost::blank &) const
      {
        PINOCCHIO_THROW_PRETTY(
          std::invalid_argument, "The constraint model is of type boost::blank.");
        return internal::NoRun<std::string>::run();
      }

      template<
        typename Scalar,
        int Options,
        template<typename S, int O> class ConstraintCollectionTpl>
      static std::string
      run(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
      {
        return boost::apply_visitor(ConstraintModelShortnameVisitor(), cmodel);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    inline std::string
    shortname(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      return ConstraintModelShortnameVisitor::run(cmodel);
    }

    /**
     * @brief      ConstraintDataShortnameVisitor visitor
     */
    struct ConstraintDataShortnameVisitor : boost::static_visitor<std::string>
    {
      template<typename ConstraintDataDerived>
      std::string operator()(const ConstraintDataBase<ConstraintDataDerived> & cdata) const
      {
        return cdata.shortname();
      }
      std::string operator()(const boost::blank &) const
      {
        PINOCCHIO_THROW_PRETTY(
          std::invalid_argument, "The constraint data is of type boost::blank.");
        return internal::NoRun<std::string>::run();
      }

      template<
        typename Scalar,
        int Options,
        template<typename S, int O> class ConstraintCollectionTpl>
      static std::string
      run(const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata)
      {
        return boost::apply_visitor(ConstraintDataShortnameVisitor(), cdata);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    inline std::string
    shortname(const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata)
    {
      return ConstraintDataShortnameVisitor::run(cdata);
    }

    ///
    /// \brief Base structure for \b Unary visitation of a ConstraintModel.
    ///        This structure provides runners to call the right visitor according to the number of
    ///        arguments.
    ///

    template<typename ConstraintModelVisitorDerived, typename ReturnType = void>
    struct ConstraintUnaryVisitorBase
    {
      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      static ReturnType run(
        ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
        ArgsTmp args)
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        ModelAndDataVisitor<ConstraintModel, ConstraintData, ArgsTmp> visitor(cdata, args);

        return boost::apply_visitor(visitor, cmodel);
      }

      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      static ReturnType run(
        const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
        ArgsTmp args)
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        ModelAndDataVisitor<ConstraintModel, ConstraintData, ArgsTmp> visitor(cdata, args);

        return boost::apply_visitor(visitor, cmodel);
      }

      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      static ReturnType run(
        const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
        ArgsTmp args)
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        ModelAndDataVisitor<ConstraintModel, const ConstraintData, ArgsTmp> visitor(cdata, args);

        return boost::apply_visitor(visitor, cmodel);
      }

      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      static ReturnType
      run(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel, ArgsTmp args)
      {
        ModelVisitor<Scalar, Options, ConstraintCollectionTpl, ArgsTmp> visitor(args);
        return boost::apply_visitor(visitor, cmodel);
      }

      template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
      static ReturnType
      run(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
      {
        ModelVisitor<Scalar, Options, ConstraintCollectionTpl, NoArg> visitor;
        return boost::apply_visitor(visitor, cmodel);
      }

      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      static ReturnType
      run(ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel, ArgsTmp args)
      {
        ModelVisitor<Scalar, Options, ConstraintCollectionTpl, ArgsTmp> visitor(args);
        return boost::apply_visitor(visitor, cmodel);
      }

      template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
      static ReturnType run(ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
      {
        ModelVisitor<Scalar, Options, ConstraintCollectionTpl, NoArg> visitor;
        return boost::apply_visitor(visitor, cmodel);
      }

    private:
      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      struct ModelVisitor : public boost::static_visitor<ReturnType>
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        ModelVisitor(ArgsTmp args)
        : args(args)
        {
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>,
            bf::append(boost::ref(cmodel.derived()), args));
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>,
            bf::append(boost::ref(cmodel.derived()), args));
        }

        template<typename ConstraintDataDerived>
        ReturnType operator()(const ConstraintDataBase<ConstraintDataDerived> & cdata) const
        {
          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintDataDerived>,
            bf::append(boost::ref(cdata.derived()), args));
        }

        ReturnType operator()(const boost::blank &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint model is of type boost::blank.");
          return internal::NoRun<ReturnType>::run();
        }

        template<typename S, int O>
        ReturnType operator()(const FictiousConstraintModelTpl<S, O> &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint model is of type FictiousConstraintModelTpl.");
          return internal::NoRun<ReturnType>::run();
        }

        template<typename S, int O>
        ReturnType operator()(const FictiousConstraintDataTpl<S, O> &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint data is of type FictiousConstraintDataTpl.");
          return internal::NoRun<ReturnType>::run();
        }

        ArgsTmp args;
      };

      template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
      struct ModelVisitor<Scalar, Options, ConstraintCollectionTpl, NoArg>
      : public boost::static_visitor<ReturnType>
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        ModelVisitor()
        {
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          return ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>(
            cmodel.derived());
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          return ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>(
            cmodel.derived());
        }

        template<typename ConstraintDataDerived>
        ReturnType operator()(const ConstraintDataBase<ConstraintDataDerived> & cdata) const
        {
          return ConstraintModelVisitorDerived::template algo<ConstraintDataDerived>(
            cdata.derived());
        }

        ReturnType operator()(const boost::blank &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint model is of type boost::blank.");
          return internal::NoRun<ReturnType>::run();
        }

        template<typename S, int O>
        ReturnType operator()(const FictiousConstraintModelTpl<S, O> &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint model is of type FictiousConstraintModelTpl.");
          return internal::NoRun<ReturnType>::run();
        }

        template<typename S, int O>
        ReturnType operator()(const FictiousConstraintDataTpl<S, O> &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint data is of type FictiousConstraintDataTpl.");
          return internal::NoRun<ReturnType>::run();
        }
      }; // struct ModelVisitor

      template<typename ConstraintModel, typename ConstraintData, typename ArgsTmp>
      struct ModelAndDataVisitor : public boost::static_visitor<ReturnType>
      {

        ModelAndDataVisitor(ConstraintData & cdata, ArgsTmp args)
        : cdata(cdata)
        , args(args)
        {
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          typedef typename ConstraintModelBase<ConstraintModelDerived>::ConstraintData
            ConstraintDataDerived;
          using ConstraintDataGet = typename std::conditional<
            std::is_const<ConstraintData>::value, const ConstraintDataDerived,
            ConstraintDataDerived>::type;

          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>,
            bf::append(
              boost::ref(cmodel.derived()), boost::ref(boost::get<ConstraintDataGet>(cdata)),
              args));
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          typedef typename ConstraintModelBase<ConstraintModelDerived>::ConstraintData
            ConstraintDataDerived;
          using ConstraintDataGet = typename std::conditional<
            std::is_const<ConstraintData>::value, const ConstraintDataDerived,
            ConstraintDataDerived>::type;

          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>,
            bf::append(
              boost::ref(cmodel.derived()), boost::ref(boost::get<ConstraintDataGet>(cdata)),
              args));
        }

        template<typename S, int O>
        ReturnType operator()(const FictiousConstraintModelTpl<S, O> &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint model is of type FictiousConstraintModelTpl.");
          return ReturnType();
        }

        ReturnType operator()(const boost::blank &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint model is of type boost::blank.");
          return internal::NoRun<ReturnType>::run();
        }

        ConstraintData & cdata;
        ArgsTmp args;
      }; // struct ModelAndDataVisitor
    };

    /**
     * @brief      ConstraintModelCalcVisitor fusion visitor
     */
    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    struct ConstraintModelCalcVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelCalcVisitor<Scalar, Options, JointCollectionTpl>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::vector<const Model &, const Data &> ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data)
      {
        cmodel.calc(model, data, cdata.derived());
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl>
    void calc(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata)
    {
      typedef ConstraintModelCalcVisitor<Scalar, Options, JointCollectionTpl> Algo;
      Algo::run(cmodel, cdata, typename Algo::ArgsType(model, data));
    }

    /**
     * @brief      ConstraintModelResizeVisitor fusion visitor
     */
    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    struct ConstraintModelResizeVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelResizeVisitor<Scalar, Options, JointCollectionTpl>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::vector<const Model &, const Data &> ArgsType;

      template<typename ConstraintModel>
      static void algo(
        pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data)
      {
        cmodel.resize(model, data, cdata.derived());
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl>
    void resize(
      ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata)
    {
      typedef ConstraintModelResizeVisitor<Scalar, Options, JointCollectionTpl> Algo;
      Algo::run(cmodel, cdata, typename Algo::ArgsType(model, data));
    }

    /**
     * @brief      ConstraintModelJacobianVisitor fusion visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      typename JacobianMatrix>
    struct ConstraintModelJacobianVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelJacobianVisitor<Scalar, Options, JointCollectionTpl, JacobianMatrix>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::vector<const Model &, const Data &, JacobianMatrix &> ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const Eigen::MatrixBase<JacobianMatrix> & jacobian_matrix)
      {
        cmodel.jacobian(model, data, cdata.derived(), jacobian_matrix.const_cast_derived());
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename JacobianMatrix>
    void jacobian(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & jacobian_matrix)
    {
      typedef ConstraintModelJacobianVisitor<Scalar, Options, JointCollectionTpl, JacobianMatrix>
        Algo;
      Algo::run(
        cmodel, cdata, typename Algo::ArgsType(model, data, jacobian_matrix.const_cast_derived()));
    }

    /**
     * @brief      ConstraintModelCreateDataVisitor fusion visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    struct ConstraintModelCreateDataVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelCreateDataVisitor<Scalar, Options, ConstraintCollectionTpl>,
        typename ConstraintCollectionTpl<Scalar, Options>::ConstraintDataVariant>
    {
      typedef NoArg ArgsType;
      typedef ConstraintCollectionTpl<Scalar, Options> ConstraintCollection;
      typedef typename ConstraintCollection::ConstraintModelVariant ConstraintModelVariant;
      typedef typename ConstraintCollection::ConstraintDataVariant ConstraintDataVariant;

      template<typename ConstraintModel>
      static ConstraintDataVariant
      algo(const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel)
      {
        return cmodel.createData();
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl>
    createData(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef ConstraintModelCreateDataVisitor<Scalar, Options, ConstraintCollectionTpl> Algo;
      return Algo::run(cmodel);
    }

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename ConstraintDataDerived>
    struct ConstraintDataComparisonOperatorVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintDataComparisonOperatorVisitor<
          Scalar,
          Options,
          ConstraintCollectionTpl,
          ConstraintDataDerived>,
        bool>
    {
      typedef boost::fusion::vector<const ConstraintDataDerived &> ArgsType;

      template<typename ConstraintData>
      static bool algo(
        const ConstraintDataBase<ConstraintData> & cdata_lhs,
        const ConstraintDataDerived & cdata_rhs)
      {
        return cdata_lhs.derived() == cdata_rhs;
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename ConstraintDataDerived>
    bool isEqual(
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata_generic,
      const ConstraintDataBase<ConstraintDataDerived> & cdata)
    {
      typedef ConstraintDataComparisonOperatorVisitor<
        Scalar, Options, ConstraintCollectionTpl, ConstraintDataDerived>
        Algo;
      return Algo::run(cdata_generic, typename Algo::ArgsType(boost::ref(cdata.derived())));
    }

    /**
     * @brief      ConstraintModelSizeVisitor visitor
     */
    template<typename Scalar, int Options>
    struct ConstraintModelSizeVisitor
    : visitors::ConstraintUnaryVisitorBase<ConstraintModelSizeVisitor<Scalar, Options>, int>
    {
      typedef NoArg ArgsType;

      template<typename ConstraintModel>
      static int algo(const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel)
      {
        return cmodel.size();
      }
    };

    template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
    int size(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef ConstraintModelSizeVisitor<Scalar, Options> Algo;
      return Algo::run(cmodel);
    }

    /**
     * @brief      ConstraintModelActiveSizeVisitor visitor
     */
    template<typename Scalar, int Options>
    struct ConstraintModelActiveSizeVisitor
    : visitors::ConstraintUnaryVisitorBase<ConstraintModelActiveSizeVisitor<Scalar, Options>, int>
    {

      typedef NoArg ArgsType;

      template<typename ConstraintModel>
      static int algo(const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel)
      {
        return cmodel.activeSize();
      }
    };

    template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
    int activeSize(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef ConstraintModelActiveSizeVisitor<Scalar, Options> Algo;
      return Algo::run(cmodel);
    }

    /**
     * @brief      ConstraintModelGetRowActivableIndexesVisitor visitor
     */
    template<typename Scalar, int Options>
    struct ConstraintModelGetRowActivableIndexesVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelGetRowActivableIndexesVisitor<Scalar, Options>,
        const std::vector<Eigen::DenseIndex> &>
    {
      typedef const std::vector<Eigen::DenseIndex> & ReturnType;

      typedef boost::fusion::vector<const Eigen::DenseIndex> ArgsType;

      template<typename ConstraintModel>
      static ReturnType algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::DenseIndex row_id)
      {
        return cmodel.getRowActivableIndexes(row_id);
      }
    };

    template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
    const std::vector<Eigen::DenseIndex> & getRowActivableIndexes(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const Eigen::DenseIndex row_id)
    {
      typedef ConstraintModelGetRowActivableIndexesVisitor<Scalar, Options> Algo;
      return Algo::run(cmodel, typename Algo::ArgsType(row_id));
    }

    /**
     * @brief      ConstraintModelGetRowActiveIndexesVisitor visitor
     */
    template<typename Scalar, int Options>
    struct ConstraintModelGetRowActiveIndexesVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelGetRowActiveIndexesVisitor<Scalar, Options>,
        const std::vector<Eigen::DenseIndex> &>
    {
      typedef const std::vector<Eigen::DenseIndex> & ReturnType;

      typedef boost::fusion::vector<const Eigen::DenseIndex> ArgsType;

      template<typename ConstraintModel>
      static ReturnType algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::DenseIndex row_id)
      {
        return cmodel.getRowActiveIndexes(row_id);
      }
    };

    template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
    const std::vector<Eigen::DenseIndex> & getRowActiveIndexes(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const Eigen::DenseIndex row_id)
    {
      typedef ConstraintModelGetRowActiveIndexesVisitor<Scalar, Options> Algo;
      return Algo::run(cmodel, typename Algo::ArgsType(row_id));
    }

    /**
     * @brief      ConstraintModelGetRowActivableSparsityPatternVisitor visitor
     */
    template<typename Scalar, int Options>
    struct ConstraintModelGetRowActivableSparsityPatternVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelGetRowActivableSparsityPatternVisitor<Scalar, Options>,
        const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> &>
    {
      typedef const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> & ReturnType;

      typedef boost::fusion::vector<const Eigen::DenseIndex> ArgsType;

      template<typename ConstraintModel>
      static ReturnType algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::DenseIndex row_id)
      {
        return cmodel.getRowActivableSparsityPattern(row_id);
      }
    };

    template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
    const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> & getRowActivableSparsityPattern(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const Eigen::DenseIndex row_id)
    {
      typedef ConstraintModelGetRowActivableSparsityPatternVisitor<Scalar, Options> Algo;
      return Algo::run(cmodel, typename Algo::ArgsType(row_id));
    }

    /**
     * @brief      ConstraintModelGetRowActiveSparsityPatternVisitor visitor
     */
    template<typename Scalar, int Options>
    struct ConstraintModelGetRowActiveSparsityPatternVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelGetRowActiveSparsityPatternVisitor<Scalar, Options>,
        const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> &>
    {
      typedef const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> & ReturnType;

      typedef boost::fusion::vector<const Eigen::DenseIndex> ArgsType;

      template<typename ConstraintModel>
      static ReturnType algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::DenseIndex row_id)
      {
        return cmodel.getRowActiveSparsityPattern(row_id);
      }
    };

    template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
    const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> & getRowActiveSparsityPattern(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const Eigen::DenseIndex row_id)
    {
      typedef ConstraintModelGetRowActiveSparsityPatternVisitor<Scalar, Options> Algo;
      return Algo::run(cmodel, typename Algo::ArgsType(row_id));
    }

    /**
     * @brief      ConstraintModelJacobianMatrixProductVisitor visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op>
    struct ConstraintModelJacobianMatrixProductVisitor
    : visitors::ConstraintUnaryVisitorBase<ConstraintModelJacobianMatrixProductVisitor<
        Scalar,
        Options,
        JointCollectionTpl,
        InputMatrix,
        OutputMatrix,
        op>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::vector<
        const Model &,
        const Data &,
        const InputMatrix &,
        OutputMatrix &,
        AssignmentOperatorTag<op>>
        ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const Eigen::MatrixBase<InputMatrix> & input_matrix,
        const Eigen::MatrixBase<OutputMatrix> & result_matrix,
        AssignmentOperatorTag<op> aot)
      {
        cmodel.jacobianMatrixProduct(
          model, data, cdata.derived(), input_matrix.derived(), result_matrix.const_cast_derived(),
          aot);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op = SETTO>
    void jacobianMatrixProduct(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix,
      const Eigen::MatrixBase<OutputMatrix> & result_matrix,
      AssignmentOperatorTag<op> aot = SetTo())
    {
      typedef ConstraintModelJacobianMatrixProductVisitor<
        Scalar, Options, JointCollectionTpl, InputMatrix, OutputMatrix, op>
        Algo;

      typename Algo::ArgsType args(
        model, data, input_matrix.derived(), result_matrix.const_cast_derived(), aot);
      Algo::run(cmodel, cdata, args);
    }

    /**
     * @brief      ConstraintModelJacobianTransposeMatrixProductVisitor visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op>
    struct ConstraintModelJacobianTransposeMatrixProductVisitor
    : visitors::ConstraintUnaryVisitorBase<ConstraintModelJacobianTransposeMatrixProductVisitor<
        Scalar,
        Options,
        JointCollectionTpl,
        InputMatrix,
        OutputMatrix,
        op>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::vector<
        const Model &,
        const Data &,
        const InputMatrix &,
        OutputMatrix &,
        AssignmentOperatorTag<op>>
        ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const Eigen::MatrixBase<InputMatrix> & input_matrix,
        const Eigen::MatrixBase<OutputMatrix> & result_matrix,
        AssignmentOperatorTag<op> aot)
      {
        cmodel.jacobianTransposeMatrixProduct(
          model, data, cdata.derived(), input_matrix.derived(), result_matrix.const_cast_derived(),
          aot);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op = SETTO>
    void jacobianTransposeMatrixProduct(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix,
      const Eigen::MatrixBase<OutputMatrix> & result_matrix,
      AssignmentOperatorTag<op> aot = SetTo())
    {
      typedef ConstraintModelJacobianTransposeMatrixProductVisitor<
        Scalar, Options, JointCollectionTpl, InputMatrix, OutputMatrix, op>
        Algo;

      typename Algo::ArgsType args(
        model, data, input_matrix.derived(), result_matrix.const_cast_derived(), aot);
      Algo::run(cmodel, cdata, args);
    }

    /**
     * @brief      AppendCouplingConstraintInertiasVisitor visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      typename VectorNLike,
      ReferenceFrame rf>
    struct AppendCouplingConstraintInertiasVisitor
    : visitors::ConstraintUnaryVisitorBase<AppendCouplingConstraintInertiasVisitor<
        Scalar,
        Options,
        JointCollectionTpl,
        VectorNLike,
        rf>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::
        vector<const Model &, Data &, const VectorNLike &, ReferenceFrameTag<rf>>
          ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        Data & data,
        const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
        const ReferenceFrameTag<rf> reference_frame)
      {
        cmodel.appendCouplingConstraintInertias(
          model, data, cdata.derived(), diagonal_constraint_inertia.derived(), reference_frame);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename VectorNLike,
      ReferenceFrame rf>
    void appendCouplingConstraintInertias(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame)
    {
      typedef AppendCouplingConstraintInertiasVisitor<
        Scalar, Options, JointCollectionTpl, VectorNLike, rf>
        Algo;

      typename Algo::ArgsType args(
        model, data, diagonal_constraint_inertia.derived(), reference_frame);
      Algo::run(cmodel, cdata, args);
    }

    /**
     * @brief      mapConstraintForceToJointSpace visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForceLike,
      class ForceAllocator,
      typename JointTorquesLike,
      ReferenceFrame rf>
    struct MapConstraintForceToJointSpaceVisitor
    : visitors::ConstraintUnaryVisitorBase<MapConstraintForceToJointSpaceVisitor<
        Scalar,
        Options,
        JointCollectionTpl,
        ConstraintForceLike,
        ForceAllocator,
        JointTorquesLike,
        rf>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef std::vector<ForceTpl<Scalar, Options>, ForceAllocator> ForceVector;
      typedef boost::fusion::vector<
        const Model &,
        const Data &,
        const ConstraintForceLike &,
        ForceVector &,
        JointTorquesLike &,
        const ReferenceFrameTag<rf>>
        ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const Eigen::MatrixBase<ConstraintForceLike> & constraint_forces,
        std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
        const Eigen::MatrixBase<JointTorquesLike> & joint_torques,
        const ReferenceFrameTag<rf> reference_frame)
      {
        cmodel.mapConstraintForceToJointSpace(
          model, data, cdata, constraint_forces, joint_forces, joint_torques.const_cast_derived(),
          reference_frame);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename ConstraintForceLike,
      class ForceAllocator,
      typename JointTorquesLike,
      ReferenceFrame rf>
    void mapConstraintForceToJointSpace(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques,
      const ReferenceFrameTag<rf> reference_frame)
    {
      typedef MapConstraintForceToJointSpaceVisitor<
        Scalar, Options, JointCollectionTpl, ConstraintForceLike, ForceAllocator, JointTorquesLike,
        rf>
        Algo;

      typename Algo::ArgsType args(
        model, data, constraint_forces.derived(), joint_forces, joint_torques.const_cast_derived(),
        reference_frame);
      Algo::run(cmodel, cdata, args);
    }

    /**
     * @brief      mapConstraintForceToJointSpace visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      class MotionAllocator,
      typename GeneralizedVelocityLike,
      typename ConstraintMotionLike,
      ReferenceFrame rf>
    struct MapJointSpaceToConstraintMotionVisitor
    : visitors::ConstraintUnaryVisitorBase<MapJointSpaceToConstraintMotionVisitor<
        Scalar,
        Options,
        JointCollectionTpl,
        MotionAllocator,
        GeneralizedVelocityLike,
        ConstraintMotionLike,
        rf>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef std::vector<MotionTpl<Scalar, Options>, MotionAllocator> MotionVector;
      typedef boost::fusion::vector<
        const Model &,
        const Data &,
        const MotionVector &,
        const GeneralizedVelocityLike &,
        ConstraintMotionLike &,
        const ReferenceFrameTag<rf>>
        ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
        const Eigen::MatrixBase<GeneralizedVelocityLike> & generalized_velocity,
        const Eigen::MatrixBase<ConstraintMotionLike> & constraint_motions,
        const ReferenceFrameTag<rf> reference_frame)
      {
        cmodel.mapJointSpaceToConstraintMotion(
          model, data, cdata, joint_motions, generalized_velocity.derived(),
          constraint_motions.const_cast_derived(), reference_frame);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      class MotionAllocator,
      typename GeneralizedVelocityLike,
      typename ConstraintMotionLike,
      ReferenceFrame rf>
    void mapJointSpaceToConstraintMotion(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
      const Eigen::MatrixBase<GeneralizedVelocityLike> & generalized_velocity,
      const Eigen::MatrixBase<ConstraintMotionLike> & constraint_motions,
      const ReferenceFrameTag<rf> reference_frame)
    {
      typedef MapJointSpaceToConstraintMotionVisitor<
        Scalar, Options, JointCollectionTpl, MotionAllocator, GeneralizedVelocityLike,
        ConstraintMotionLike, rf>
        Algo;

      typename Algo::ArgsType args(
        model, data, joint_motions, generalized_velocity.derived(),
        constraint_motions.const_cast_derived(), reference_frame);
      Algo::run(cmodel, cdata, args);
    }

    /**
     * @brief      ConstraintModelComplianceVisitor visitor
     */
    template<typename ReturnType>
    struct ConstraintModelComplianceVisitor
    : ConstraintUnaryVisitorBase<ConstraintModelComplianceVisitor<ReturnType>, ReturnType>
    {
      typedef NoArg ArgsType;

      template<typename ConstraintModelDerived>
      static ReturnType algo(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.compliance();
      }
      template<typename ConstraintModelDerived>
      static ReturnType algo(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.compliance();
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    typename traits<
      ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::ComplianceVectorTypeConstRef
    compliance(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
        ComplianceVectorTypeConstRef ReturnType;
      typedef ConstraintModelComplianceVisitor<ReturnType> Algo;
      return Algo::run(cmodel);
    }

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    typename traits<
      ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::ComplianceVectorTypeRef
    compliance(ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
        ComplianceVectorTypeRef ReturnType;
      typedef ConstraintModelComplianceVisitor<ReturnType> Algo;
      return Algo::run(cmodel);
    }

    /**
     * @brief      ConstraintModelActiveComplianceVisitor visitor
     */
    template<typename ReturnType>
    struct ConstraintModelActiveComplianceVisitor
    : ConstraintUnaryVisitorBase<ConstraintModelActiveComplianceVisitor<ReturnType>, ReturnType>
    {
      typedef NoArg ArgsType;

      template<typename ConstraintModelDerived>
      static ReturnType algo(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.getActiveCompliance();
      }
      template<typename ConstraintModelDerived>
      static ReturnType algo(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.getActiveCompliance();
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    typename traits<
      ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::ComplianceVectorTypeConstRef
    getActiveCompliance(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
        ComplianceVectorTypeConstRef ReturnType;
      typedef ConstraintModelActiveComplianceVisitor<ReturnType> Algo;
      return Algo::run(cmodel);
    }

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    typename traits<
      ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::ComplianceVectorTypeRef
    getActiveCompliance(ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
        ComplianceVectorTypeRef ReturnType;
      typedef ConstraintModelActiveComplianceVisitor<ReturnType> Algo;
      return Algo::run(cmodel);
    }

    /// \brief BaumgarteCorrectorVectorParametersGetter - default behavior for false for
    /// HasBaumgarteCorrectorVector
    template<
      bool HasBaumgarteCorrectorVector,
      typename BaumgarteVector,
      typename BaumgarteVectorReturnType>
    struct BaumgarteCorrectorVectorParametersGetter
    {
      template<typename ConstraintModelDerived>
      static BaumgarteVectorReturnType
      run(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        std::stringstream ss;
        ss << cmodel.shortname() << " does not have baumgarte vector corrector parameters.\n";
        PINOCCHIO_THROW(std::invalid_argument, ss.str());
        return internal::NoRun<BaumgarteVectorReturnType>::run();
      }
      template<typename ConstraintModelDerived>
      static BaumgarteVectorReturnType run(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        std::stringstream ss;
        ss << cmodel.shortname() << " does not have baumgarte vector corrector parameters.\n";
        PINOCCHIO_THROW(std::invalid_argument, ss.str());
        return internal::NoRun<BaumgarteVectorReturnType>::run();
      }
    };

    /// \brief BaumgarteCorrectorVectorParametersGetter - partial specialization for true for
    /// HasBaumgarteCorrectorVector
    template<typename BaumgarteVector, typename BaumgarteVectorReturnType>
    struct BaumgarteCorrectorVectorParametersGetter<
      true,
      BaumgarteVector,
      BaumgarteVectorReturnType>
    {
      template<typename ConstraintModelDerived>
      static BaumgarteVectorReturnType
      run(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.baumgarte_corrector_vector_parameters().template ref<BaumgarteVector>();
      }
      template<typename ConstraintModelDerived>
      static BaumgarteVectorReturnType run(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.baumgarte_corrector_vector_parameters().template ref<BaumgarteVector>();
      }
    };

    /**
     * @brief      BaumgarteCorrectorVectorParametersVisitor visitor
     */
    template<typename BaumgarteVectorType, typename BaumgarteVectorReturnType>
    struct BaumgarteCorrectorVectorParametersVisitor
    : ConstraintUnaryVisitorBase<
        BaumgarteCorrectorVectorParametersVisitor<BaumgarteVectorType, BaumgarteVectorReturnType>,
        BaumgarteVectorReturnType>
    {
      typedef NoArg ArgsType;

      template<typename ConstraintModelDerived>
      static BaumgarteVectorReturnType
      algo(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        static constexpr bool has_baumgarte_corrector_vector =
          traits<ConstraintModelDerived>::has_baumgarte_corrector_vector;
        return BaumgarteCorrectorVectorParametersGetter<
          has_baumgarte_corrector_vector, BaumgarteVectorType,
          BaumgarteVectorReturnType>::run(cmodel);
      }
      template<typename ConstraintModelDerived>
      static BaumgarteVectorReturnType algo(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        static constexpr bool has_baumgarte_corrector_vector =
          traits<ConstraintModelDerived>::has_baumgarte_corrector_vector;
        return BaumgarteCorrectorVectorParametersGetter<
          has_baumgarte_corrector_vector, BaumgarteVectorType,
          BaumgarteVectorReturnType>::run(cmodel);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
      BaumgarteCorrectorVectorParametersConstRef
      getBaumgarteCorrectorVectorParameters(
        const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
        BaumgarteVectorType BaumgarteVectorType;
      typedef typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
        BaumgarteCorrectorVectorParametersConstRef BaumgarteCorrectorVectorParametersConstRef;
      return BaumgarteCorrectorVectorParametersVisitor<
        BaumgarteVectorType, BaumgarteCorrectorVectorParametersConstRef>::run(cmodel);
    }

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
      BaumgarteCorrectorVectorParametersRef
      getBaumgarteCorrectorVectorParameters(
        ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
        BaumgarteVectorType BaumgarteVectorType;
      typedef typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
        BaumgarteCorrectorVectorParametersRef BaumgarteCorrectorVectorParametersRef;
      return BaumgarteCorrectorVectorParametersVisitor<
        BaumgarteVectorType, BaumgarteCorrectorVectorParametersRef>::run(cmodel);
    }

    /// \brief BaumgarteCorrectorParametersGetter - default behavior for false for
    /// HasBaumgarteCorrector
    template<bool HasBaumgarteCorrector, typename BaumgarteReturnType>
    struct BaumgarteCorrectorParametersGetter
    {
      template<typename ConstraintModelDerived>
      static BaumgarteReturnType run(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        std::stringstream ss;
        ss << cmodel.shortname() << " does not have baumgarte corrector parameters.\n";
        PINOCCHIO_THROW(std::invalid_argument, ss.str());
        return internal::NoRun<BaumgarteReturnType>::run();
      }
      template<typename ConstraintModelDerived>
      static BaumgarteReturnType run(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        std::stringstream ss;
        ss << cmodel.shortname() << " does not have baumgarte corrector parameters.\n";
        PINOCCHIO_THROW(std::invalid_argument, ss.str());
        return internal::NoRun<BaumgarteReturnType>::run();
      }
    };

    /// \brief BaumgarteCorrectorParametersGetter - partial specialization for true for
    /// HasBaumgarteCorrector
    template<typename BaumgarteReturnType>
    struct BaumgarteCorrectorParametersGetter<true, BaumgarteReturnType>
    {
      template<typename ConstraintModelDerived>
      static BaumgarteReturnType run(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.baumgarte_corrector_parameters();
      }
      template<typename ConstraintModelDerived>
      static BaumgarteReturnType run(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.baumgarte_corrector_parameters();
      }
    };

    /**
     * @brief      BaumgarteCorrectorParametersVisitor visitor
     */
    template<typename BaumgarteReturnType>
    struct BaumgarteCorrectorParametersVisitor
    : ConstraintUnaryVisitorBase<
        BaumgarteCorrectorParametersVisitor<BaumgarteReturnType>,
        BaumgarteReturnType>
    {
      typedef NoArg ArgsType;

      template<typename ConstraintModelDerived>
      static BaumgarteReturnType algo(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        static constexpr bool has_baumgarte_corrector =
          traits<ConstraintModelDerived>::has_baumgarte_corrector;
        return BaumgarteCorrectorParametersGetter<
          has_baumgarte_corrector, BaumgarteReturnType>::run(cmodel);
      }

      template<typename ConstraintModelDerived>
      static BaumgarteReturnType algo(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        static constexpr bool has_baumgarte_corrector =
          traits<ConstraintModelDerived>::has_baumgarte_corrector;
        return BaumgarteCorrectorParametersGetter<
          has_baumgarte_corrector, BaumgarteReturnType>::run(cmodel);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    const BaumgarteCorrectorParametersTpl<Scalar> & getBaumgarteCorrectorParameters(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;
      typedef const BaumgarteCorrectorParameters & ReturnType;
      typedef BaumgarteCorrectorParametersVisitor<ReturnType> Algo;
      return Algo::run(cmodel);
    }

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    BaumgarteCorrectorParametersTpl<Scalar> & getBaumgarteCorrectorParameters(
      ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;
      typedef BaumgarteCorrectorParameters & ReturnType;
      typedef BaumgarteCorrectorParametersVisitor<ReturnType> Algo;
      return Algo::run(cmodel);
    }

  } // namespace visitors

} // namespace pinocchio

#endif // ifdef __pinocchio_algorithm_constraints_constraint_model_visitor_hpp__
