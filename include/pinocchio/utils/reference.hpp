//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_utils_reference_hpp__
#define __pinocchio_utils_reference_hpp__

#include <functional>
#include <memory>

namespace pinocchio
{
  
  template<typename T>
  T & get_ref(const std::reference_wrapper<T> & ref) { return ref.get(); }
  template<typename T>
  const T & get_ref(const std::reference_wrapper<const T> & ref) { return ref.get(); }
  
  template<typename T>
  T & get_ref(const std::shared_ptr<T> & ptr) { return *ptr; }
  template<typename T>
  const T & get_ref(const std::shared_ptr<const T> & ptr) { return *ptr; }

}

#endif // ifndef __pinocchio_utils_reference_hpp__
