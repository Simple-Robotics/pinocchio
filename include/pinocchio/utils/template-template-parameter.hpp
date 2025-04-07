//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_utils_template_template_parameter_hpp__
#define __pinocchio_utils_template_template_parameter_hpp__

namespace pinocchio
{
  namespace internal
  {
    template<typename C>
    struct extract_template_template_parameter;

    template<template<class...> class C, class... Parameters>
    struct extract_template_template_parameter<C<Parameters...>>
    {
      template<class... Other>
      using type = C<Other...>;
    };

  } // namespace internal
} // namespace pinocchio

#endif // __pinocchio_utils_template_template_parameter_hpp__
