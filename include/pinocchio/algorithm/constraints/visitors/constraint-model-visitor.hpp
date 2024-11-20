//
// Copyright (c) 2023-2024 INRIA
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
          // Hacky way to not have to return something real. The system should throw before.
          const typename std::remove_reference<T>::type * null_ptr = NULL;
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
        const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
        ArgsTmp args)
      {
        InternalVisitorModelAndData<Scalar, Options, ConstraintCollectionTpl, ArgsTmp> visitor(
          cdata, args);
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
        InternalVisitorModel<Scalar, Options, ConstraintCollectionTpl, ArgsTmp> visitor(args);
        return boost::apply_visitor(visitor, cmodel);
      }

      template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
      static ReturnType
      run(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
      {
        InternalVisitorModel<Scalar, Options, ConstraintCollectionTpl, NoArg> visitor;
        return boost::apply_visitor(visitor, cmodel);
      }

    private:
      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      struct InternalVisitorModel : public boost::static_visitor<ReturnType>
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        InternalVisitorModel(ArgsTmp args)
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
      struct InternalVisitorModel<Scalar, Options, ConstraintCollectionTpl, NoArg>
      : public boost::static_visitor<ReturnType>
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        InternalVisitorModel()
        {
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
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
      };

      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        InternalVisitorModelAndData(ConstraintData & cdata, ArgsTmp args)
        : cdata(cdata)
        , args(args)
        {
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>,
            bf::append(
              boost::ref(cmodel.derived()),
              boost::ref(
                boost::get<typename ConstraintModelBase<ConstraintModelDerived>::ConstraintData>(
                  cdata)),
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
      };
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
      typedef visitors::NoArg ArgsType;
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
      return ConstraintModelCreateDataVisitor<Scalar, Options, ConstraintCollectionTpl>::run(
        cmodel);
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
     * @brief      ConstraintModelGetRowSparsityPatternVisitor visitor
     */
    template<typename Scalar, int Options>
    struct ConstraintModelGetRowSparsityPatternVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelGetRowSparsityPatternVisitor<Scalar, Options>,
        const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> &>
    {
      typedef const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> & ReturnType;

      typedef boost::fusion::vector<const Eigen::DenseIndex> ArgsType;

      template<typename ConstraintModel>
      static ReturnType algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::DenseIndex row_id)
      {
        return cmodel.getRowSparsityPattern(row_id);
      }
    };

    template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
    const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> & getRowSparsityPattern(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const Eigen::DenseIndex row_id)
    {
      typedef ConstraintModelGetRowSparsityPatternVisitor<Scalar, Options> Algo;
      return Algo::run(cmodel, typename Algo::ArgsType(row_id));
    }

  } // namespace visitors

} // namespace pinocchio

#endif // ifdef __pinocchio_algorithm_constraints_constraint_model_visitor_hpp__
