//
// Copyright (c) 2015-2025 CNRS INRIA
//

#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/constraints-variant.hpp"
#include "pinocchio/bindings/python/utils/std-aligned-vector.hpp"

namespace pinocchio
{
  namespace python
  {
    void exposeConstraints()
    {
      typedef context::ConstraintCollectionDefault::ConstraintModelVariant ConstraintModelVariant;
      boost::mpl::for_each<ConstraintModelVariant::types>(ConstraintModelExposer());
      bp::to_python_converter<ConstraintModelVariant, ConstraintVariantVisitor<ConstraintModelVariant>>();
      ConstraintModelPythonVisitor<context::ConstraintModel>::expose();
      StdAlignedVectorPythonVisitor<context::ConstraintModel>::expose("StdVec_ConstraintModelVector");

      typedef context::ConstraintCollectionDefault::ConstraintDataVariant ConstraintDataVariant;
      boost::mpl::for_each<ConstraintDataVariant::types>(ConstraintDataExposer());
      bp::to_python_converter<ConstraintDataVariant, ConstraintVariantVisitor<ConstraintDataVariant>>();
      ConstraintDataPythonVisitor<context::ConstraintData>::expose();
      StdAlignedVectorPythonVisitor<context::ConstraintData>::expose("StdVec_ConstraintDataVector");

      // Special vector exposition for usage in Simple
      // TODO: A foreach ?
      StdAlignedVectorPythonVisitor<context::BilateralPointConstraintModel>::expose("StdVec_BilateralPointConstraintModel");
      StdAlignedVectorPythonVisitor<context::BilateralPointConstraintData>::expose("StdVec_BilateralPointConstraintData");
      StdAlignedVectorPythonVisitor<context::WeldConstraintModel>::expose("StdVec_WeldConstraintModel");
      StdAlignedVectorPythonVisitor<context::WeldConstraintData>::expose("StdVec_WeldConstraintData");
    }
  } // namespace python
} // namespace pinocchio
