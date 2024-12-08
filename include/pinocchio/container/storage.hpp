//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_container_storage_hpp__
#define __pinocchio_container_storage_hpp__

#include "pinocchio/fwd.hpp"

namespace pinocchio
{

  template<typename Scalar, int MaxSizeAtCompileTime = Eigen::Dynamic>
  struct EigenStorageTpl;

  template<typename _Scalar, int _MaxSizeAtCompileTime>
  struct EigenStorageTpl
  {
    typedef _Scalar Scalar;
    enum
    {
      MaxSizeAtCompileTime = _MaxSizeAtCompileTime
    };

    typedef Eigen::Matrix<Scalar, MaxSizeAtCompileTime, 1> StorageVector;

    /// \brief Default constructor from a given maximum size.
    ///
    /// \param[in] max_size Size of the allocated memory chunk (max_size * sizeof(Scalar)).
    explicit EigenStorageTpl(const Eigen::DenseIndex max_size)
    : m_storage(max_size)
    {
    }

    /// \brief Resize the current capacity of the internal storage.
    ///
    /// \remarks The resizing only happens when the new_size is greater than the current capacity
    void resize(const Eigen::DenseIndex new_size)
    {
      if (new_size > capacity())
        m_storage.resize(2 * new_size); // Double the size of the storage
    }

    /// \brief Conservative resize of the current capacity of the internal storage. The data are
    /// kepts in memory.
    ///
    /// \remarks The resizing only happens when the new_size is greater than the current capacity
    void conservativeResize(const Eigen::DenseIndex new_size)
    {
      if (new_size > capacity())
        m_storage.conservativeResize(2 * new_size); // Double the size of the storage
    }

    template<typename MatrixLike>
    Eigen::Map<typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike)>
    as(const Eigen::DenseIndex rows, const Eigen::DenseIndex cols)
    {
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike) PlainType;
      typedef Eigen::Map<PlainType> ReturnType;

      const Eigen::DenseIndex size = rows * cols;
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        size <= capacity(), "You ask for a matrix that exceeds the storage capacity");

      return ReturnType(m_storage.data(), rows, cols);
    }

    template<typename MatrixLike>
    const Eigen::Map<const typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike)>
    as(const Eigen::DenseIndex rows, const Eigen::DenseIndex cols) const
    {
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike) PlainType;
      typedef const Eigen::Map<const PlainType> ReturnType;

      const Eigen::DenseIndex size = rows * cols;
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        size <= capacity(), "You ask for a matrix that exceeds the storage capacity");

      return ReturnType(m_storage.data(), rows, cols);
    }

    template<typename VectorLike>
    Eigen::Map<typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike)> as(const Eigen::DenseIndex size)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike);
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike) PlainType;
      typedef Eigen::Map<PlainType> ReturnType;

      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        size <= capacity(), "You ask for a vector that exceeds the storage capacity");

      return ReturnType(m_storage.data(), size);
    }

    template<typename VectorLike>
    const Eigen::Map<const typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike)>
    as(const Eigen::DenseIndex size) const
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike);
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike) PlainType;
      typedef const Eigen::Map<const PlainType> ReturnType;

      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        size <= capacity(), "You ask for a vector that exceeds the storage capacity");

      return ReturnType(m_storage.data(), size);
    }

    ///  \brief Returns the size of the storage space currently allocated.
    Eigen::DenseIndex capacity() const
    {
      return m_storage.size();
    }

    /// \brief Returns a const reference of the internal storage.
    const StorageVector & storage() const
    {
      return m_storage;
    }

    /// \brief Returns the internal pointer to the data.
    const Scalar * data() const
    {
      return m_storage.data();
    }

  protected:
    /// \brief Internal vector containing the stored quantities
    StorageVector m_storage;
  }; // struct EigenStorageTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_container_storage_hpp__
