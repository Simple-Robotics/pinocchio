//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_utils_reference_hpp__
#define __pinocchio_utils_reference_hpp__

#include <functional>
#include <memory>

namespace pinocchio
{
  namespace helper
  {
    // std::reference_wrapper
    template<typename T>
    T & get_ref(const std::reference_wrapper<T> & ref)
    {
      return ref.get();
    }
    template<typename T>
    const T & get_ref(const std::reference_wrapper<const T> & ref)
    {
      return ref.get();
    }

    template<typename T>
    T * get_pointer(const std::reference_wrapper<T> & ref)
    {
      return &ref.get();
    }
    template<typename T>
    const T * get_pointer(const std::reference_wrapper<const T> & ref)
    {
      return &ref.get();
    }

    // std::shared_ptr
    template<typename T>
    T & get_ref(const std::shared_ptr<T> & ptr)
    {
      return *ptr;
    }
    template<typename T>
    const T & get_ref(const std::shared_ptr<const T> & ptr)
    {
      return *ptr;
    }

    template<typename T>
    T * get_pointer(const std::shared_ptr<T> & ptr)
    {
      return ptr.get();
    }
    template<typename T>
    const T * get_pointer(const std::shared_ptr<const T> & ptr)
    {
      return ptr.get();
    }

  } // namespace helper
} // namespace pinocchio

#endif // ifndef __pinocchio_utils_reference_hpp__
