//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_container_storage_hpp__
#define __pinocchio_container_storage_hpp__

#include "pinocchio/fwd.hpp"

namespace pinocchio
{

  template<typename MatrixLike>
  struct EigenStorageTpl
  {
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike) PlainMatrixType;
    typedef typename MatrixLike::Scalar Scalar;

    typedef Eigen::Map<PlainMatrixType> MapType;
    typedef Eigen::Map<PlainMatrixType> & RefMapType;
    typedef const Eigen::Map<PlainMatrixType> & ConstRefMapType;
    typedef const Eigen::Map<const PlainMatrixType> ConstMapType;
    typedef Eigen::DenseIndex Index;

    enum
    {
      MaxSizeAtCompileTime =
        ((PlainMatrixType::MaxRowsAtCompileTime != Eigen::Dynamic)
         && (PlainMatrixType::MaxRowsAtCompileTime != Eigen::Dynamic))
          ? PlainMatrixType::MaxRowsAtCompileTime * PlainMatrixType::MaxColsAtCompileTime
          : Eigen::Dynamic,
      IsVectorAtCompileTime = MatrixLike::IsVectorAtCompileTime,
      Options = PlainMatrixType::Options
    };

    typedef Eigen::Matrix<Scalar, MaxSizeAtCompileTime, 1, Options> StorageVector;

    /// \brief Default constructor from given matrix dimension (rows, cols) and maximum rows and
    /// columns
    ///
    /// \param[in] rows Number of rows
    /// \param[in] cols Number of columns
    /// \param[in] max_rows Maximum number of rows
    /// \param[in] max_cols Maximum number of columns
    ///
    EigenStorageTpl(const Index rows, const Index cols, const Index max_rows, const Index max_cols)
    : m_storage(max_rows * max_cols)
    , m_map(MapType(m_storage.data(), rows, cols))
    {
    }

    /// \brief Default constructor from given matrix dimension (rows, cols).
    ///
    /// \param[in] rows Number of rows.
    /// \param[in] cols Number of columns.
    ///
    EigenStorageTpl(const Index rows, const Index cols)
    : m_storage(Eigen::DenseIndex(rows * cols))
    , m_map(MapType(m_storage.data(), rows, cols))
    {
    }

    /// \brief Resize the current capacity of the internal storage.
    ///
    /// \remarks The resizing only happens when the new_size is greater than the current capacity
    void resize(const Index rows, const Index cols)
    {
      const Index new_size = rows * cols;
      if (new_size > capacity())
        m_storage.resize(2 * new_size); // Double the size of the storage
      new (&m_map) MapType(m_storage.data(), rows, cols);
    }

    /// \brief Conservative resize of the current capacity of the internal storage. The data are
    /// kepts in memory.
    ///
    /// \remarks The resizing only happens when the new_size is greater than the current capacity
    void conservativeResize(const Index rows, const Index cols)
    {
      const Index old_rows = this->rows(), old_cols = this->cols();
      const PlainMatrixType copy(map()); // save current value in the storage
      this->resize(rows, cols);

      // copy back values
      const Index min_rows = (std::min)(rows, old_rows), min_cols = (std::min)(cols, old_cols);
      map().topLeftCorner(min_rows, min_cols) = copy.topLeftCorner(min_rows, min_cols);
    }

    ///  \brief Returns the size of the storage space currently allocated.
    Index capacity() const
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

    /// \brief Returns a map toward the internal matrix.
    ConstRefMapType map() const
    {
      return m_map;
    }

    /// \brief Returns a map toward the internal matrix.
    RefMapType map()
    {
      return m_map;
    }

    /// \brief Returns the number of rows
    Index rows() const
    {
      return map().rows();
    }

    /// \brief Returns the number of columns
    Index cols() const
    {
      return map().cols();
    }

    ///  \brief Returns the size of the underlying matrix or vector.
    Index size() const
    {
      return map().size();
    }

  protected:
    /// \brief Internal vector containing the stored quantities
    StorageVector m_storage;

    /// \brief Map
    MapType m_map;
  }; // struct EigenStorageTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_container_storage_hpp__
