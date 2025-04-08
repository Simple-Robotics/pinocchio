//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_utils_hpp__
#define __pinocchio_algorithm_constraints_utils_hpp__

#include "pinocchio/utils/std-vector.hpp"

namespace pinocchio
{

  template<typename ConstraintModel, class ConstraintAllocator>
  typename internal::template std_vector_with_same_allocator<
    std::vector<ConstraintModel, ConstraintAllocator>>::
    template type<typename ConstraintModel::ConstraintData>
    createData(const std::vector<ConstraintModel, ConstraintAllocator> & constraint_models)
  {
    typedef typename internal::template std_vector_with_same_allocator<
      std::vector<ConstraintModel, ConstraintAllocator>>::
      template type<typename ConstraintModel::ConstraintData>
        ReturnType;

    ReturnType constraint_datas(constraint_models.size());

    for (const auto & cm : constraint_models)
      constraint_datas.push_back(cm.createData());

    return constraint_datas;
  }

} // namespace pinocchio

#endif // __pinocchio_algorithm_constraints_utils_hpp__
