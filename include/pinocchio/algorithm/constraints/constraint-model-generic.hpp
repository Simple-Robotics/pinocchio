//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraint_model_generic_hpp__
#define __pinocchio_algorithm_constraint_model_generic_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-generic.hpp"
#include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hpp"

namespace pinocchio
{

  template<
    typename _Scalar,
    int _Options,
    template<typename S, int O> class ConstraintCollectionTpl>
  struct traits<ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;
    typedef ConstraintData Data;
    typedef boost::blank ConstraintSet;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef VectorXs VectorConstraintSize;

    typedef VectorXs ComplianceVectorType;
    typedef Eigen::Ref<ComplianceVectorType> ComplianceVectorTypeRef;
    typedef Eigen::Ref<const ComplianceVectorType> ComplianceVectorTypeConstRef;

    template<typename InputMatrix>
    struct JacobianMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::
        Matrix<Scalar, Eigen::Dynamic, InputMatrix::ColsAtCompileTime, InputMatrixPlain::Options>
          type;
    };

    template<typename InputMatrix>
    struct JacobianTransposeMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::Matrix<
        Scalar,
        Eigen::Dynamic,
        InputMatrixPlain::ColsAtCompileTime,
        InputMatrixPlain::Options>
        type;
    };
  };

  template<
    typename _Scalar,
    int _Options,
    template<typename S, int O> class ConstraintCollectionTpl>
  struct ConstraintModelTpl
  : ConstraintModelBase<ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl>>
  , ConstraintCollectionTpl<_Scalar, _Options>::ConstraintModelVariant
  , serialization::Serializable<ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> Self;
    typedef ConstraintModelBase<Self> Base;
    typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;
    typedef ConstraintCollectionTpl<Scalar, Options> ConstraintCollection;
    typedef typename ConstraintCollection::ConstraintDataVariant ConstraintDataVariant;
    typedef typename ConstraintCollection::ConstraintModelVariant ConstraintModelVariant;
    typedef typename traits<Self>::ComplianceVectorType ComplianceVectorType;
    typedef typename traits<Self>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Self>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;

    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    ConstraintModelTpl()
    : ConstraintModelVariant()
    {
    }

    ConstraintModelTpl(const ConstraintModelVariant & cmodel_variant)
    : ConstraintModelVariant(cmodel_variant)
    {
    }

    template<typename ContraintModelDerived>
    ConstraintModelTpl(const ConstraintModelBase<ContraintModelDerived> & cmodel)
    : ConstraintModelVariant((ConstraintModelVariant)cmodel.derived())
    {
      BOOST_MPL_ASSERT(
        (boost::mpl::contains<typename ConstraintModelVariant::types, ContraintModelDerived>));
    }

    ConstraintData createData() const
    {
      return ::pinocchio::visitors::createData<Scalar, Options, ConstraintCollectionTpl>(*this);
    }

    template<template<typename, int> class JointCollectionTpl>
    void calc(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      ::pinocchio::visitors::calc(*this, model, data, cdata);
    }

    template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobian(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & jacobian_matrix) const
    {
      ::pinocchio::visitors::jacobian(
        *this, model, data, cdata, jacobian_matrix.const_cast_derived());
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    const EigenIndexVector & getRowActiveIndexes(const Eigen::DenseIndex row_id) const
    {
      return ::pinocchio::visitors::getRowActiveIndexes(*this, row_id);
    }

    /// \brief Returns the sparsity pattern associated with a given row
    const BooleanVector & getRowSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      return ::pinocchio::visitors::getRowSparsityPattern(*this, row_id);
    }

    /// \brief Runs the underlying jacobian multiplication with a matrix.
    template<
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix>
    typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix) const
    {
      typedef typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
        ReturnType;
      ReturnType res(size(), input_matrix.cols());
      jacobianMatrixProduct(model, data, cdata, input_matrix.derived(), res);
      return res;
    }

    /// \brief Runs the underlying jacobian multiplication with a matrix.
    template<
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op = SETTO>
    void jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix,
      const Eigen::MatrixBase<OutputMatrix> & result_matrix,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      ::pinocchio::visitors::jacobianMatrixProduct(
        *this, model, data, cdata, input_matrix.derived(), result_matrix.const_cast_derived(), aot);
    }

    /// \brief Runs the underlying jacobian transpose multiplication with a matrix.
    template<
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix>
    typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix) const
    {
      typedef
        typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
          ReturnType;
      ReturnType res(model.nv, input_matrix.cols());
      jacobianTransposeMatrixProduct(*this, model, data, cdata, input_matrix.derived(), res);
      return res;
    }

    template<
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op = SETTO>
    void jacobianTransposeMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix,
      const Eigen::MatrixBase<OutputMatrix> & result_matrix,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      ::pinocchio::visitors::jacobianTransposeMatrixProduct(
        *this, model, data, cdata, input_matrix.derived(), result_matrix.const_cast_derived(), aot);
    }

    static std::string classname()
    {
      return "ConstraintModel";
    }
    std::string shortname() const
    {
      return ::pinocchio::visitors::shortname(*this);
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ComplianceVectorTypeConstRef compliance() const
    {
      return ::pinocchio::visitors::compliance<ComplianceVectorTypeConstRef>(*this);
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ComplianceVectorTypeRef compliance()
    {
      return ::pinocchio::visitors::compliance<ComplianceVectorTypeRef>(*this);
    }

    /// \brief Returns the size of the constraint
    int size() const
    {
      return ::pinocchio::visitors::size(*this);
    }

    boost::blank & set()
    {
      static boost::blank val;
      PINOCCHIO_THROW_PRETTY(
        std::runtime_error, "Set method is not accessible for ConstraintModelTpl.");
      return val;
    }

    const boost::blank & set() const
    {
      static boost::blank val;
      PINOCCHIO_THROW_PRETTY(
        std::runtime_error, "Set method is not accessible for ConstraintModelTpl.");
      return val;
    }

    ConstraintModelVariant & toVariant()
    {
      return static_cast<ConstraintModelVariant &>(*this);
    }

    const ConstraintModelVariant & toVariant() const
    {
      return static_cast<const ConstraintModelVariant &>(*this);
    }

    template<typename ConstraintModelDerived>
    bool isEqual(const ConstraintModelBase<ConstraintModelDerived> & other) const
    {
      return ::pinocchio::isEqual(*this, other.derived());
    }

    bool isEqual(const ConstraintModelTpl & other) const
    {
      return toVariant() == other.toVariant();
    }

    bool operator==(const ConstraintModelTpl & other) const
    {
      return isEqual(other);
    }

    bool operator!=(const ConstraintModelTpl & other) const
    {
      return !(*this == other);
    }
  }; // struct ConstraintModelTpl
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraint_model_generic_hpp__
