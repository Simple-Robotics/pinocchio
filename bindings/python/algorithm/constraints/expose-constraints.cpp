//
// Copyright (c) 2025 INRIA
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
      boost::mpl::for_each<ConstraintModelVariant::types>(ConstraintStdVectorExposer());
      bp::to_python_converter<
        ConstraintModelVariant, ConstraintVariantVisitor<ConstraintModelVariant>>();
      ConstraintModelPythonVisitor<context::ConstraintModel>::expose();
      StdAlignedVectorPythonVisitor<context::ConstraintModel>::expose("StdVec_ConstraintModel");

      typedef context::ConstraintCollectionDefault::ConstraintDataVariant ConstraintDataVariant;
      boost::mpl::for_each<ConstraintDataVariant::types>(ConstraintDataExposer());
      boost::mpl::for_each<ConstraintDataVariant::types>(ConstraintStdVectorExposer());
      bp::to_python_converter<
        ConstraintDataVariant, ConstraintVariantVisitor<ConstraintDataVariant>>();
      ConstraintDataPythonVisitor<context::ConstraintData>::expose();
      StdAlignedVectorPythonVisitor<context::ConstraintData>::expose("StdVec_ConstraintData");
    }
  } // namespace python
} // namespace pinocchio
