//
// Copyright (c) 2017-2018 CNRS
// Copyright (c) 2018-2025 INRIA
//

#ifndef __pinocchio_serialization_eigen_matrix_hpp__
#define __pinocchio_serialization_eigen_matrix_hpp__

#include "pinocchio/math/tensor.hpp"

#include <boost/serialization/split_free.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>

// If hpp-fcl < 3.0.0 The GCC Eigen/Boost.Serialization workaround
// is already defined.
// If we don't link against hpp-fcl or hpp-fcl >= 3.0.0 then we must define
// the workaround.
#if defined PINOCCHIO_WITH_HPP_FCL
  #include <hpp/fcl/config.hh>
  #if !HPP_FCL_VERSION_AT_LEAST(3, 0, 0) // hpp-fcl < 3.0.0
    #define HPP_FCL_SKIP_EIGEN_BOOST_SERIALIZATION
    #include <hpp/fcl/serialization/eigen.h>
  #else // hpp-fcl >= 3.0.0
    // Workaround a bug in GCC >= 7 and C++17.
    // ref. https://gitlab.com/libeigen/eigen/-/issues/1676
    #ifdef __GNUC__
      #if __GNUC__ >= 7 && __cplusplus >= 201703L
namespace boost
{
  namespace serialization
  {
    struct U;
  }
} // namespace boost
namespace Eigen
{
  namespace internal
  {
    template<>
    struct traits<boost::serialization::U>
    {
      enum
      {
        Flags = 0
      };
    };
  } // namespace internal
} // namespace Eigen
      #endif
    #endif
  #endif
#else // !PINOCCHIO_WITH_HPP_FCL
  // Workaround a bug in GCC >= 7 and C++17.
  // ref. https://gitlab.com/libeigen/eigen/-/issues/1676
  #ifdef __GNUC__
    #if __GNUC__ >= 7 && __cplusplus >= 201703L
namespace boost
{
  namespace serialization
  {
    struct U;
  }
} // namespace boost
namespace Eigen
{
  namespace internal
  {
    template<>
    struct traits<boost::serialization::U>
    {
      enum
      {
        Flags = 0
      };
    };
  } // namespace internal
} // namespace Eigen
    #endif
  #endif
#endif

// Similar workaround but for MSVC when C++17 is enabled.
// TODO Find _MSC_VER range.
#if (defined(_MSVC_LANG) && _MSVC_LANG >= 201703)
namespace boost
{
  namespace archive
  {
    class binary_iarchive;
    class xml_iarchive;
    class text_iarchive;
  } // namespace archive
} // namespace boost
namespace Eigen
{
  namespace internal
  {
    template<>
    struct traits<boost::archive::binary_iarchive>
    {
      enum
      {
        Flags = 0
      };
    };
    template<>
    struct traits<boost::archive::xml_iarchive>
    {
      enum
      {
        Flags = 0
      };
    };
    template<>
    struct traits<boost::archive::text_iarchive>
    {
      enum
      {
        Flags = 0
      };
    };
  } // namespace internal
} // namespace Eigen
#endif // MSVC with C++17

namespace boost
{
  namespace serialization
  {
    namespace internal
    {
      namespace Eigen
      {
        template<class Archive, typename EigenPlainObjectBase>
        void serialize_eigen_plain_object(
          Archive & ar, EigenPlainObjectBase & m, const unsigned int /*version*/)
        {

          ::Eigen::DenseIndex rows(m.rows()), cols(m.cols());
          if (EigenPlainObjectBase::RowsAtCompileTime == ::Eigen::Dynamic)
            ar & BOOST_SERIALIZATION_NVP(rows);
          if (EigenPlainObjectBase::ColsAtCompileTime == ::Eigen::Dynamic)
            ar & BOOST_SERIALIZATION_NVP(cols);

          if (Archive::is_loading::value)
            m.resize(rows, cols);
          ar & make_nvp("data", make_array(m.data(), (size_t)m.size()));
        }
      } // namespace Eigen
    } // namespace internal

    template<class Archive, typename Derived>
    void serialize(Archive & ar, ::Eigen::PlainObjectBase<Derived> & m, const unsigned int version)
    {
      internal::Eigen::serialize_eigen_plain_object(ar, m, version);
    }

    template<
      class Archive,
      typename Scalar,
      int Rows,
      int Cols,
      int Options,
      int MaxRows,
      int MaxCols>
    void serialize(
      Archive & ar,
      ::Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> & m,
      const unsigned int version)
    {
      typedef ::Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> Self;
      typedef typename Self::Base Base;
      serialize(ar, static_cast<Base &>(m), version);
    }

    template<class Archive, typename PlainObjectBase, int MapOptions, typename StrideType>
    void serialize(
      Archive & ar,
      ::Eigen::Map<PlainObjectBase, MapOptions, StrideType> & m,
      const unsigned int version)
    {
      internal::Eigen::serialize_eigen_plain_object(ar, m, version);
    }

    template<
      class Archive,
      typename Scalar,
      int Rows,
      int Cols,
      int Options,
      int MaxRows,
      int MaxCols>
    void serialize(
      Archive & ar,
      ::Eigen::Array<Scalar, Rows, Cols, Options, MaxRows, MaxCols> & m,
      const unsigned int version)
    {
      typedef ::Eigen::Array<Scalar, Rows, Cols, Options, MaxRows, MaxCols> Self;
      typedef typename Self::Base Base;
      serialize(ar, static_cast<Base &>(m), version);
    }

#if !defined(PINOCCHIO_WITH_EIGEN_TENSOR_MODULE)                                                   \
  && ((__cplusplus <= 199711L && EIGEN_COMP_MSVC < 1900) || defined(__CUDACC__) || defined(EIGEN_AVOID_STL_ARRAY))
    template<class Archive, typename _IndexType, std::size_t _NumIndices>
    void save(
      Archive & ar, const Eigen::array<_IndexType, _NumIndices> & a, const unsigned int /*version*/)
    {
      ar & make_nvp("array", make_array(&a.front(), _NumIndices));
    }

    template<class Archive, typename _IndexType, std::size_t _NumIndices>
    void
    load(Archive & ar, Eigen::array<_IndexType, _NumIndices> & a, const unsigned int /*version*/)
    {
      ar >> make_nvp("array", make_array(&a.front(), _NumIndices));
    }

    template<class Archive, typename _IndexType, std::size_t _NumIndices>
    void
    serialize(Archive & ar, Eigen::array<_IndexType, _NumIndices> & a, const unsigned int version)
    {
      split_free(ar, a, version);
    }
#else
    template<class Archive, class T, std::size_t N>
    void save(Archive & ar, const std::array<T, N> & a, const unsigned int version)
    {
      typedef std::array<T, N> Array;
      serialize(ar, const_cast<Array &>(a), version);
    }

    template<class Archive, class T, std::size_t N>
    void load(Archive & ar, std::array<T, N> & a, const unsigned int version)
    {
      serialize(ar, a, version);
    }
#endif

#ifdef PINOCCHIO_WITH_EIGEN_TENSOR_MODULE

    template<class Archive, typename _IndexType, int _NumIndices>
    void save(
      Archive & ar, const Eigen::DSizes<_IndexType, _NumIndices> & ds, const unsigned int version)
    {
      save(ar, static_cast<const Eigen::array<_IndexType, _NumIndices> &>(ds), version);
    }

    template<class Archive, typename _IndexType, int _NumIndices>
    void load(Archive & ar, Eigen::DSizes<_IndexType, _NumIndices> & ds, const unsigned int version)
    {
      load(ar, static_cast<Eigen::array<_IndexType, _NumIndices> &>(ds), version);
    }

    template<class Archive, typename _IndexType, int _NumIndices>
    void
    serialize(Archive & ar, Eigen::DSizes<_IndexType, _NumIndices> & ds, const unsigned int version)
    {
      split_free(ar, static_cast<Eigen::array<_IndexType, _NumIndices> &>(ds), version);
    }

#endif

    template<class Archive, typename _Scalar, int _NumIndices, int _Options, typename _IndexType>
    void save(
      Archive & ar,
      const ::pinocchio::Tensor<_Scalar, _NumIndices, _Options, _IndexType> & t,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::Tensor<_Scalar, _NumIndices, _Options, _IndexType> Tensor;
      const typename Tensor::Dimensions & dimensions = t.dimensions();

      ar & BOOST_SERIALIZATION_NVP(dimensions);
      ar & make_nvp("data", make_array(t.data(), (size_t)t.size()));
    }

    template<class Archive, typename _Scalar, int _NumIndices, int _Options, typename _IndexType>
    void load(
      Archive & ar,
      ::pinocchio::Tensor<_Scalar, _NumIndices, _Options, _IndexType> & t,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::Tensor<_Scalar, _NumIndices, _Options, _IndexType> Tensor;
      typename Tensor::Dimensions dimensions;

      ar >> BOOST_SERIALIZATION_NVP(dimensions);
      t.resize(dimensions);

      ar >> make_nvp("data", make_array(t.data(), (size_t)t.size()));
    }

    template<class Archive, typename _Scalar, int _NumIndices, int _Options, typename _IndexType>
    void serialize(
      Archive & ar,
      ::pinocchio::Tensor<_Scalar, _NumIndices, _Options, _IndexType> & t,
      const unsigned int version)
    {
      split_free(ar, t, version);
    }

  } // namespace serialization
} // namespace boost

#endif // ifndef __pinocchio_serialization_eigen_matrix_hpp__
