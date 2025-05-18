//
// Copyright (c) 2024-2025 INRIA
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
    T * get_pointer(const std::shared_ptr<T> & ptr)
    {
      return ptr.get();
    }
    template<typename T>
    const T * get_pointer(const std::shared_ptr<const T> & ptr)
    {
      return ptr.get();
    }

    template<typename T>
    std::reference_wrapper<T> make_ref(T & value)
    {
      return std::reference_wrapper<T>(value);
    }

    template<typename T>
    std::reference_wrapper<const T> make_ref(const T & value)
    {
      return std::reference_wrapper<const T>(value);
    }

    template<typename T>
    struct remove_ref
    {
      typedef T type;
      static T & get_ref(T & v)
      {
        return v;
      }
    };

    //    template<typename T>
    //    struct remove_ref<const T>
    //    {
    //      typedef const T type;
    //      static const T & get_ref(const T & v)
    //      {
    //        return v;
    //      }
    //    };

    template<typename T>
    struct remove_ref<std::reference_wrapper<T>>
    {
      typedef typename remove_ref<T>::type type;

      static T & get_ref(const std::reference_wrapper<T> & ref)
      {
        return ref.get();
      }
    };

    template<typename T>
    struct remove_ref<const std::reference_wrapper<const T>>
    {
      typedef typename remove_ref<const T>::type type;

      static const T & get_ref(const std::reference_wrapper<const T> & ref)
      {
        return ref.get();
      }
    };

    template<typename T>
    struct remove_ref<std::shared_ptr<T>>
    {
      typedef typename remove_ref<T>::type type;

      static T & get_ref(const std::shared_ptr<T> & ptr)
      {
        return *ptr;
      }
    };

    template<typename T>
    struct remove_ref<const std::shared_ptr<T>> : remove_ref<std::shared_ptr<T>>
    {
    };

    template<typename T>
    struct remove_ref<std::shared_ptr<const T>>
    {
      typedef typename remove_ref<const T>::type type;

      static const T & get_ref(const std::shared_ptr<const T> & ptr)
      {
        return *ptr;
      }
    };

    template<typename T>
    struct remove_ref<const std::shared_ptr<const T>> : remove_ref<std::shared_ptr<const T>>
    {
    };

    template<typename T>
    struct remove_ref<std::unique_ptr<T>>
    {
      typedef typename remove_ref<T>::type type;

      static T & get_ref(const std::unique_ptr<T> & ptr)
      {
        return *ptr;
      }
    };

    template<typename T>
    struct remove_ref<const std::unique_ptr<T>> : remove_ref<std::unique_ptr<T>>
    {
    };

    template<typename T>
    struct remove_ref<std::unique_ptr<const T>>
    {
      typedef typename remove_ref<const T>::type type;

      static const T & get_ref(const std::unique_ptr<const T> & ptr)
      {
        return *ptr;
      }
    };

    template<typename T>
    struct remove_ref<const std::unique_ptr<const T>> : remove_ref<std::unique_ptr<const T>>
    {
    };

    template<typename T>
    typename remove_ref<T>::type & get_ref(T & v)
    {
      return remove_ref<T>::get_ref(v);
    }

    template<typename T>
    const typename remove_ref<const T>::type & get_ref(const T & v)
    {
      return remove_ref<const T>::get_ref(v);
    }

    template<typename T>
    struct is_type_holder
    {
      static constexpr bool value = false;
    };

    template<typename T>
    struct is_type_holder<std::reference_wrapper<T>>
    {
      static constexpr bool value = true;
    };

    template<typename T>
    struct is_type_holder<std::shared_ptr<T>>
    {
      static constexpr bool value = true;
    };

  } // namespace helper
} // namespace pinocchio

#endif // ifndef __pinocchio_utils_reference_hpp__
