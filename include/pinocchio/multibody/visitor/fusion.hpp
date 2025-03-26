//
// Copyright (c) 2015-2023 CNRS INRIA
//

#ifndef __pinocchio_multibody_visitor_fusion_hpp__
#define __pinocchio_multibody_visitor_fusion_hpp__

#define BOOST_FUSION_INVOKE_MAX_ARITY 12

#include <boost/variant/static_visitor.hpp>
#include <boost/fusion/include/invoke.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>

namespace pinocchio
{
  namespace fusion
  {

    namespace bf = boost::fusion;
    typedef boost::blank NoArg;

  } // namespace fusion
} // namespace pinocchio

namespace boost
{
  namespace fusion
  {

    // Declare struct that takes at least 1 type
    template<typename T, typename... Ts>
    struct AppendReturnType;

    // Specialization for one type V: V is the return type
    template<typename V>
    struct AppendReturnType<V>
    {
      typedef V type;
    };

    // Specialization for two or more types: return as pushfront of the first type on the return
    // type of the remainings
    template<typename T1, typename T, typename... Ts>
    struct AppendReturnType<T1, T, Ts...>
    {
      typedef
        typename result_of::push_front<typename AppendReturnType<T, Ts...>::type const, T1>::type
          type;
    };

    // Append of one type
    template<typename V>
    V append(V const & v)
    {
      return v;
    }

    // Append of two or more types
    template<typename T1, typename T, typename... Ts>
    typename AppendReturnType<T1, T, Ts...>::type
    append(T1 const & t1, T const & t, Ts const &... ts)
    {
      return push_front(append(t, ts...), t1);
    }
  } // namespace fusion
} // namespace boost

#endif // ifndef __pinocchio_multibody_visitor_fusion_hpp__
