//
// Copyright (c) 2018-2025 INRIA
// Copyright (c) 2017-2018 CNRS
//

#ifndef __pinocchio_macros_hpp__
#define __pinocchio_macros_hpp__

#include "pinocchio/deprecated.hpp"
#include "pinocchio/warning.hpp"

// On Windows, __cplusplus is not necessarily set to the C++ version being used.
// See https://docs.microsoft.com/fr-fr/cpp/build/reference/zc-cplusplus?view=vs-2019 for further
// information.

#if (__cplusplus >= 202302L)
  #define PINOCCHIO_WITH_CXX23_SUPPORT
#endif

#if (__cplusplus >= 202002L)
  #define PINOCCHIO_WITH_CXX20_SUPPORT
#endif

#if (__cplusplus >= 201703L || (defined(_MSVC_LANG) && _MSVC_LANG >= 201703))
  #define PINOCCHIO_WITH_CXX17_SUPPORT
#endif

#if (__cplusplus >= 201402L || (defined(_MSVC_LANG) && _MSVC_LANG >= 201403))
  #define PINOCCHIO_WITH_CXX14_SUPPORT
#endif

#if (__cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1600))
  #define PINOCCHIO_WITH_CXX11_SUPPORT
#endif

#if defined(__APPLE__) && defined(__aarch64__)
  #define PINOCCHIO_MAC_ARM64
#endif

#define PINOCCHIO_NOEXCEPT noexcept

/// \brief Function attribute to forbid inlining.
/// This is a compiler hint that can be not respected.
#define PINOCCHIO_DONT_INLINE EIGEN_DONT_INLINE

// This macro can be used to prevent from macro expansion, similarly to EIGEN_NOT_A_MACRO
#define PINOCCHIO_NOT_A_MACRO

namespace pinocchio
{
  namespace helper
  {
    template<typename T>
    struct argument_type;
    template<typename T, typename U>
    struct argument_type<T(U)>
    {
      typedef U type;
    };
  } // namespace helper
} // namespace pinocchio

/// \brief Empty macro argument
#define PINOCCHIO_MACRO_EMPTY_ARG

/// \brief Helper to declare that a parameter is unused
#define PINOCCHIO_UNUSED_VARIABLE(var) (void)(var)
#ifndef NDEBUG
  #define PINOCCHIO_ONLY_USED_FOR_DEBUG(var)
#else
  #define PINOCCHIO_ONLY_USED_FOR_DEBUG(var) PINOCCHIO_UNUSED_VARIABLE(var)
#endif

#ifdef PINOCCHIO_WITH_CXX17_SUPPORT
  #define PINOCCHIO_MAYBE_UNUSED [[maybe_unused]]
#else
  #define PINOCCHIO_MAYBE_UNUSED
#endif

#define PINOCCHIO_DEPRECATED_MOVED_HEADER(old_header, new_header)                                  \
  PINOCCHIO_DEPRECATED_HEADER("#old_header has been replaced by #new_header.")

/// Ensure that a matrix (or vector) is of correct size (compile-time and run-time assertion)
#define PINOCCHIO_ASSERT_MATRIX_SPECIFIC_SIZE(type, M, nrows, ncols)                               \
  EIGEN_STATIC_ASSERT(                                                                             \
    (type::RowsAtCompileTime == Eigen::Dynamic || type::RowsAtCompileTime == nrows)                \
      && (type::ColsAtCompileTime == Eigen::Dynamic || type::ColsAtCompileTime == ncols),          \
    THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);                                          \
  assert(M.rows() == nrows && M.cols() == ncols);

/// Static assertion.
/// \param condition a boolean convertible expression
/// \param msg a valid C++ variable name.
#define PINOCCHIO_STATIC_ASSERT(condition, msg)                                                    \
  {                                                                                                \
    int msg[(condition) ? 1 : -1]; /*avoid unused-variable warning*/                               \
    (void)msg;                                                                                     \
  }

namespace pinocchio
{
  namespace helper
  {
    template<typename D, template<typename> class TypeAccess>
    struct handle_return_type_without_typename
    {
      typedef typename TypeAccess<typename argument_type<void(D)>::type>::type type;
    };
  } // namespace helper
} // namespace pinocchio

/// \brief macros for pragma push/pop/ignore deprecated warnings
#if defined(__GNUC__) || defined(__clang__)
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH _Pragma("GCC diagnostic push")
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_POP _Pragma("GCC diagnostic pop")
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS                          \
    _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_VARIADIC_MACROS                                    \
    _Pragma("GCC diagnostic ignored \"-Wvariadic-macros\"")
  #if defined(__clang__)
    #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_MAYBE_UNINITIALIZED
    #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_SELF_ASSIGN_OVERLOADED                           \
      _Pragma("GCC diagnostic ignored \"-Wself-assign-overloaded\"")
  #else
    #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_MAYBE_UNINITIALIZED                              \
      _Pragma("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
    #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_SELF_ASSIGN_OVERLOADED
  #endif
#elif defined(WIN32)
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH _Pragma("warning(push)")
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_POP _Pragma("warning(pop)")
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS                          \
    _Pragma("warning(disable : 4996)")
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_VARIADIC_MACROS
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_SELF_ASSIGN_OVERLOADED
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_MAYBE_UNINITIALIZED
#else
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_POP
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_VARIADIC_MACROS
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_SELF_ASSIGN_OVERLOADED
  #define PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_MAYBE_UNINITIALIZED
#endif // __GNUC__ || __clang__

// Handle explicitely the GCC boring warning: 'anonymous variadic macros were introduced in C++11'
#include <exception>
#include <stdexcept>
#include <sstream>

#if defined(__GNUC__)
  #pragma GCC system_header
#endif

PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_VARIADIC_MACROS

#if WIN32
  #define PINOCCHIO_PRETTY_FUNCTION __FUNCSIG__
#else
  #define PINOCCHIO_PRETTY_FUNCTION __PRETTY_FUNCTION__
#endif

/// \brief Generic macro to throw an exception in Pinocchio if the condition is not met with a given
/// input message.
#if !defined(PINOCCHIO_NO_THROW)
  #define PINOCCHIO_THROW(exception_type, message)                                                 \
    {                                                                                              \
      throw exception_type(message);                                                               \
    }
  #define PINOCCHIO_THROW_IF(condition, exception_type, message)                                   \
    if (condition)                                                                                 \
    {                                                                                              \
      PINOCCHIO_THROW(exception_type, message);                                                    \
    }

  #define PINOCCHIO_THROW_PRETTY(exception, message)                                               \
    {                                                                                              \
      std::stringstream ss;                                                                        \
      ss << "From file: " << __FILE__ << "\n";                                                     \
      ss << "in function: " << PINOCCHIO_PRETTY_FUNCTION << "\n";                                  \
      ss << "at line: " << __LINE__ << "\n";                                                       \
      ss << "message:\n" << message << "\n";                                                       \
      throw exception(ss.str());                                                                   \
    }
  #define PINOCCHIO_THROW_PRETTY_IF(condition, exception_type, message)                            \
    if (condition)                                                                                 \
    {                                                                                              \
      PINOCCHIO_THROW_PRETTY(exception_type, message);                                             \
    }
#else
  #define PINOCCHIO_THROW(exception_type, message)
  #define PINOCCHIO_THROW_IF(condition, exception_type, message)
  #define PINOCCHIO_THROW_PRETTY(exception, message)
  #define PINOCCHIO_THROW_PRETTY_IF(condition, exception, message)
#endif

#define _PINOCCHIO_EXPAND(x) x
#define _PINOCCHIO_GET_OVERRIDE_FOR_CHECK_INPUT_ARGUMENT(_1, _2, MACRO_NAME, ...) MACRO_NAME

#define _PINOCCHIO_CHECK_INPUT_ARGUMENT_2(condition, message)                                      \
  {PINOCCHIO_THROW_IF(!(condition), std::invalid_argument, message)}

#define _PINOCCHIO_CHECK_INPUT_ARGUMENT_1(condition)                                               \
  _PINOCCHIO_CHECK_INPUT_ARGUMENT_2(                                                               \
    condition, "The following check on the input argument has failed: " #condition)

#define _PINOCCHIO_CHECK_INPUT_ARGUMENT_0

/// \brief Macro to check an assert-like condition and throw a std::invalid_argument exception (with
/// a message) if violated.
#define PINOCCHIO_CHECK_INPUT_ARGUMENT(...)                                                        \
  _PINOCCHIO_EXPAND(_PINOCCHIO_EXPAND(_PINOCCHIO_GET_OVERRIDE_FOR_CHECK_INPUT_ARGUMENT(            \
    __VA_ARGS__, _PINOCCHIO_CHECK_INPUT_ARGUMENT_2, _PINOCCHIO_CHECK_INPUT_ARGUMENT_1,             \
    _PINOCCHIO_CHECK_INPUT_ARGUMENT_0))(__VA_ARGS__))

#define _PINOCCHIO_GET_OVERRIDE_FOR_CHECK_ARGUMENT_SIZE(_1, _2, _3, MACRO_NAME, ...) MACRO_NAME

#define _PINOCCHIO_CHECK_ARGUMENT_SIZE_3(size, expected_size, message)                             \
  if (size != expected_size)                                                                       \
  {                                                                                                \
    std::ostringstream oss;                                                                        \
    oss << "wrong argument size: expected " << expected_size << ", got " << size << std::endl;     \
    oss << "hint: " << message << std::endl;                                                       \
    PINOCCHIO_THROW(std::invalid_argument, oss.str());                                             \
  }

#define _PINOCCHIO_CHECK_ARGUMENT_SIZE_2(size, expected_size)                                      \
  _PINOCCHIO_CHECK_ARGUMENT_SIZE_3(                                                                \
    size, expected_size,                                                                           \
    PINOCCHIO_WARN_STRINGISE(size) " is different from " PINOCCHIO_WARN_STRINGISE(expected_size))

#define _PINOCCHIO_CHECK_ARGUMENT_SIZE_1

/// \brief Macro to check if the size of an element is equal to the expected size.
#define PINOCCHIO_CHECK_ARGUMENT_SIZE(...)                                                         \
  _PINOCCHIO_EXPAND(_PINOCCHIO_EXPAND(_PINOCCHIO_GET_OVERRIDE_FOR_CHECK_ARGUMENT_SIZE(             \
    __VA_ARGS__, _PINOCCHIO_CHECK_ARGUMENT_SIZE_3, _PINOCCHIO_CHECK_ARGUMENT_SIZE_2,               \
    _PINOCCHIO_CHECK_ARGUMENT_SIZE_1))(__VA_ARGS__))

/// \brief Macro to check whether two matrices have the same size.
#define PINOCCHIO_CHECK_SAME_MATRIX_SIZE(mat1, mat2)                                               \
  if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())                                    \
  {                                                                                                \
    std::ostringstream oss;                                                                        \
    oss << "wrong matrix size: expected (" << mat1.rows() << ", " << mat1.cols() << "), got ("     \
        << mat2.rows() << ", " << mat2.cols() << ")" << std::endl;                                 \
    PINOCCHIO_THROW_PRETTY(std::invalid_argument, oss.str());                                      \
  }

/// \brief Macro to check whether a given matrix is square
#define PINOCCHIO_CHECK_SQUARE_MATRIX(mat)                                                         \
  if (mat.rows() != mat.cols())                                                                    \
  {                                                                                                \
    std::ostringstream oss;                                                                        \
    oss << "the matrix is not square: expected (" << mat.rows() << " == " << mat.cols()            \
        << "), got (" << mat.rows() << " != " << mat.cols() << ")" << std::endl;                   \
    PINOCCHIO_THROW_PRETTY(std::invalid_argument, oss.str());                                      \
  }

PINOCCHIO_COMPILER_DIAGNOSTIC_POP

#endif // ifndef __pinocchio_macros_hpp__
