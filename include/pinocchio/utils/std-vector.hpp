//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_utils_std_vector_hpp__
#define __pinocchio_utils_std_vector_hpp__

#include <vector>

#include "pinocchio/utils/template-template-parameter.hpp"

namespace pinocchio
{
  namespace internal
  {
    template<typename V>
    struct std_vector_extract_allocator_type
    {
      template<typename T>
      using type =
        typename extract_template_template_parameter<typename V::allocator_type>::template type<T>;
    };

    template<typename V>
    struct std_vector_with_same_allocator
    {
      template<typename T>
      using allocator_type = typename std_vector_extract_allocator_type<V>::template type<T>;

      template<typename T>
      using type = std::vector<T, allocator_type<T>>;
    };

  } // namespace internal
} // namespace pinocchio

#endif // __pinocchio_utils_std_vector_hpp__
