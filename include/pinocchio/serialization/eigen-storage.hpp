//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_serialization_eigen_storage_hpp__
#define __pinocchio_serialization_eigen_storage_hpp__

#include "pinocchio/serialization/eigen.hpp"

#include "pinocchio/container/storage.hpp"

namespace boost
{
  namespace serialization
  {

    namespace internal
    {
      template<typename MatrixLike>
      struct EigenStorageAccessor : public ::pinocchio::EigenStorageTpl<MatrixLike>
      {
        typedef ::pinocchio::EigenStorageTpl<MatrixLike> Base;
        using Base::m_map;
        using Base::m_storage;
      };
    } // namespace internal

    template<typename Archive, typename MatrixLike>
    void serialize(
      Archive & ar,
      ::pinocchio::EigenStorageTpl<MatrixLike> & storage,
      const unsigned int /*version*/)
    {
      Eigen::Index rows = storage.rows();
      Eigen::Index cols = storage.cols();
      ar & make_nvp("rows", rows);
      ar & make_nvp("cols", cols);

      typedef internal::EigenStorageAccessor<MatrixLike> Accessor;
      Accessor & storage_ = static_cast<Accessor &>(storage);
      ar & make_nvp("m_storage", storage_.m_storage);

      if (Archive::is_loading::value)
      {
        storage.resize(rows, cols); // reset internal map to point to storage
      }
    }

  } // namespace serialization
} // namespace boost

#endif // __pinocchio_serialization_eigen_storage_hpp__
