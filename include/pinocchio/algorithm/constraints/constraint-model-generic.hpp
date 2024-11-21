//
// Copyright (c) 2023-2024 INRIA
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
    typedef boost::blank ConstraintSet;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef VectorXs VectorConstraintSize;

    typedef VectorXs ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;
  };

  template<
    typename _Scalar,
    int _Options,
    template<typename S, int O> class ConstraintCollectionTpl>
  struct ConstraintModelTpl
  : ConstraintModelBase<ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl>>
  , ConstraintCollectionTpl<_Scalar, _Options>::ConstraintModelVariant
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef ConstraintModelBase<ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl>>
      Base;
    typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;
    typedef ConstraintCollectionTpl<Scalar, Options> ConstraintCollection;
    typedef typename ConstraintCollection::ConstraintDataVariant ConstraintDataVariant;
    typedef typename ConstraintCollection::ConstraintModelVariant ConstraintModelVariant;

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
  }; // struct ConstraintModelTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraint_model_generic_hpp__
