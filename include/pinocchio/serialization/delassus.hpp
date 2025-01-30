//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_delassus_serialization_hpp__
#define __pinocchio_algorithm_delassus_serialization_hpp__

#include "pinocchio/serialization/eigen.hpp"
#include "pinocchio/algorithm/delassus-operator.hpp"

namespace boost
{
  namespace serialization
  {

    template<class Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::DelassusOperatorBase<Derived> & delassus,
      const unsigned int version)
    {
      // Nothing to do yet
      PINOCCHIO_UNUSED_VARIABLE(ar);
      PINOCCHIO_UNUSED_VARIABLE(delassus);
      PINOCCHIO_UNUSED_VARIABLE(version);
    }

    namespace internal
    {

      template<typename Scalar, int Options>
      struct DelassusOperatorDenseAccessor
      : public ::pinocchio::DelassusOperatorDenseTpl<Scalar, Options>
      {
        typedef ::pinocchio::DelassusOperatorDenseTpl<Scalar, Options> Base;
        using Base::damping;
        using Base::delassus_matrix;
        using Base::llt;
        using Base::mat_tmp;
      };

    } // namespace internal

    template<class Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::DelassusOperatorDenseTpl<Scalar, Options> & delassus,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::DelassusOperatorDenseTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      typedef typename Self::CholeskyDecomposition CholeskyDecomposition;
      typedef typename Self::Matrix Matrix;
      typedef ::pinocchio::DelassusOperatorDenseTpl<Scalar, Options> Self;

      ar & make_nvp("base", boost::serialization::base_object<Base>(delassus));

      typedef internal::DelassusOperatorDenseAccessor<Scalar, Options> Accessor;
      auto & delassus_ = reinterpret_cast<Accessor &>(delassus);
      ar & make_nvp("delassus_matrix", delassus_.delassus_matrix);
      ar & make_nvp("damping", delassus_.damping);

      if (Archive::is_loading::value)
      {
        delassus_.llt = CholeskyDecomposition(delassus_.delassus_matrix);
        delassus_.mat_tmp =
          Matrix(delassus_.delassus_matrix.rows(), delassus_.delassus_matrix.cols());
      }
    }

  } // namespace serialization
} // namespace boost

#endif // __pinocchio_algorithm_delassus_serialization_hpp__
