//
// Copyright (c) 2023 INRIA
//

#ifndef __pinocchio_algorithm_constraint_model_base_hpp__
#define __pinocchio_algorithm_constraint_model_base_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/algorithm/fwd.hpp"

namespace pinocchio
{

  template<class Derived>
  struct ConstraintModelBase : NumericalBase<Derived>
  {
    typedef typename traits<Derived>::Scalar Scalar;
    enum
    {
      Options = traits<Derived>::Options
    };
    typedef typename traits<Derived>::ConstraintData ConstraintData;
    typedef typename traits<Derived>::ConstraintSet ConstraintSet;

    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> BooleanVector;
    //    typedef Eigen::Matrix<Eigen::DenseIndex,Eigen::Dynamic,1,Options> IndexVector;
    typedef std::vector<Eigen::DenseIndex> EigenIndexVector;

    Derived & derived()
    {
      return static_cast<Derived &>(*this);
    }
    const Derived & derived() const
    {
      return static_cast<const Derived &>(*this);
    }

    template<typename NewScalar>
    typename CastType<NewScalar, Derived>::type cast() const
    {
      return derived().template cast<NewScalar>();
    }

    template<typename OtherDerived>
    void cast(ConstraintModelBase<OtherDerived> & other) const
    {
      other.name = name;
    }

    /// \brief Evaluate the constraint values at the current state given by data and store the
    /// results in cdata.
    template<int Options, template<typename, int> class JointCollectionTpl>
    void calc(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      derived().calc(model, data, cdata);
    }

    template<int Options, template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobian(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & jacobian_matrix) const
    {
      derived().jacobian(model, data, cdata, jacobian_matrix.const_cast_derived());
    }

    // Attributes common to all constraints

    /// \brief Name of the constraint
    std::string name;

    template<typename OtherDerived>
    bool operator==(const ConstraintModelBase<OtherDerived> & other) const
    {
      return name == other.name;
    }

    template<typename OtherDerived>
    ConstraintModelBase & operator=(const ConstraintModelBase<OtherDerived> & other)
    {
      name = other.name;

      return *this;
    }

    ConstraintData createData() const
    {
      return derived().createData();
    }

    /// \brief Returns the colwise sparsity associated with a given row
    const BooleanVector & getRowSparsityPattern(const Eigen::Index row_id) const
    {
      return derived().getRowSparsityPattern(row_id);
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    const EigenIndexVector & getRowActiveIndexes(const Eigen::DenseIndex row_id) const
    {
      return derived().getRowActiveIndexes(row_id);
    }

    int size() const
    {
      return derived().size();
    }

    ConstraintSet & set()
    {
      return derived().set();
    }
    const ConstraintSet & set() const
    {
      return derived().set();
    }

  protected:
    template<int Options, template<typename, int> class JointCollectionTpl>
    explicit ConstraintModelBase(const ModelTpl<Scalar, Options, JointCollectionTpl> & /*model*/)
    {
    }

    /// \brief Default constructor
    ConstraintModelBase()
    {
    }

    ConstraintModelBase & base()
    {
      return *this;
    }

    const ConstraintModelBase & base() const
    {
      return *this;
    }
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraint_model_base_hpp__
