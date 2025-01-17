//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_container_double_entry_container_hpp__
#define __pinocchio_container_double_entry_container_hpp__

#include "pinocchio/fwd.hpp"
#include <vector>

namespace pinocchio
{
  namespace container
  {
    template<typename T, class Allocator = std::allocator<T>>
    struct DoubleEntryContainer;

    template<typename T, class Allocator>
    struct DoubleEntryContainer
    {
      typedef Eigen::Array<long, Eigen::Dynamic, Eigen::Dynamic> Array;
      typedef std::vector<T, Allocator> Vector;
      typedef Eigen::Index Index;
      typedef std::pair<Index, Index> IndexPair;
      typedef typename Vector::iterator iterator;
      typedef typename Vector::const_iterator const_iterator;

      /// \brief Empty constructor
      DoubleEntryContainer()
      {
      }

      /// \brief Default contructor from two array dimension
      DoubleEntryContainer(const Index nrows, const Index ncols)
      : m_keys(Array::Constant(nrows, ncols, -1))
      {
      }

      /// \brief Equality comparison operator
      bool operator==(const DoubleEntryContainer & other) const
      {
        return (m_keys == other.m_keys).all() && m_values == other.m_values;
      }

      bool operator!=(const DoubleEntryContainer & other) const
      {
        return !(*this == other);
      }

      /// \brief Returns the number of rows of the double entry table.
      Eigen::Index rows() const
      {
        return m_keys.rows();
      }

      /// \brief Returns the number of columns of the double entry table.
      Eigen::Index cols() const
      {
        return m_keys.cols();
      }

      /// \brief Returns the current size of the double entry table.
      size_t size() const
      {
        return m_values.size();
      }

      /// \brief Clear the content of the double entry table.
      void clear()
      {
        m_keys.fill(-1);
        m_values.clear();
      }

      /// \brief Fill with the same input value all the elements of the double entry table.
      void fill(const T & new_value)
      {
        std::fill(begin(), end(), new_value);
      }

      /// \brief Returns a const reference to the array of keys.
      const Array & keys() const
      {
        return m_keys;
      }

      /// \brief Returns a reference to the array of keys.
      Array & keys()
      {
        return m_keys;
      }

      /// \brief Returns the vector of values contained in the double entry table.
      const Vector & values() const
      {
        return m_values;
      }
      /// \brief Returns the vector of values contained in the double entry table.
      Vector & values()
      {
        return m_values;
      }

      ///  \brief Returns true if the key (entry1,entry2) has been succesfully added.
      bool insert(const Index entry1, const Index entry2, const T & value = T())
      {
        if (!(entry1 >= 0 && entry1 < rows()) || !(entry2 >= 0 && entry2 < cols()))
          return false;
        if (m_keys(entry1, entry2) >= 0)
          return false;

        m_keys(entry1, entry2) = long(m_values.size());
        m_values.push_back(value);
        return true;
      }

      bool insert(const IndexPair & key, const T & value = T())
      {
        return this->insert(key.first, key.second, value);
      }

      ///  \brief Returns true if the key (entry1,entry2) has been succesfully removed.
      bool remove(const Index entry1, const Index entry2)
      {
        if (!(entry1 >= 0 && entry1 < rows()) || !(entry2 >= 0 && entry2 < cols()))
          return false;
        const long index = m_keys(entry1, entry2);
        if (index < 0)
          return false;

        m_values.erase(std::next(m_values.begin(), index));
        m_keys.coeffRef(entry1, entry2) = -1;

        typedef Eigen::Array<long, Eigen::Dynamic, 1> OneDArray;
        typedef Eigen::Map<OneDArray> OneDArrayMap;
        OneDArrayMap keys_map = OneDArrayMap(m_keys.data(), m_keys.size(), 1);

        for (Eigen::Index elt_id = 0; elt_id < keys_map.size(); ++elt_id)
        {
          if (keys_map[elt_id] > index)
          {
            keys_map[elt_id]--;
          }
        }

        return true;
      }

      bool remove(const IndexPair & key)
      {
        return this->remove(key.first, key.second);
      }

      /// \brief Finds the element associated with the given input key (entry1,entry2).
      /// \returns an iterator to the element associated with the input key if it exists.
      iterator find(const Index entry1, const Index entry2)
      {
        if (!(entry1 >= 0 && entry1 < rows()) || !(entry2 >= 0 && entry2 < cols()))
          return m_values.end();

        const long index = m_keys(entry1, entry2);
        if (index < 0)
          return m_values.end();

        return std::next(m_values.begin(), index);
      }

      iterator find(const IndexPair & key)
      {
        return this->find(key.first, key.second);
      }

      /// \brief Finds the element associated with the given input key (entry1,entry2).
      /// \returns an iterator to the element associated with the input key if it exists.
      const_iterator find(const Index entry1, const Index entry2) const
      {
        if (!(entry1 >= 0 && entry1 < rows()) || !(entry2 >= 0 && entry2 < cols()))
          return m_values.end();

        const long index = m_keys(entry1, entry2);
        if (index < 0)
          return m_values.end();

        return std::next(m_values.begin(), index);
      }

      const_iterator find(const IndexPair & key) const
      {
        return this->find(key.first, key.second);
      }

      /// \brief Check whether the key (entry1,entry2) exists.
      bool exist(const Index entry1, const Index entry2) const
      {
        if (!(entry1 >= 0 && entry1 < rows()) || !(entry2 >= 0 && entry2 < cols()))
          return false;
        if (m_keys(entry1, entry2) < 0)
          return false;

        return true;
      }

      bool exist(const IndexPair & key) const
      {
        return exist(key.first, key.second);
      }

      T & operator[](const IndexPair & key)
      {
        const Index entry1 = key.first;
        const Index entry2 = key.second;

        if (!this->exist(entry1, entry2))
          this->insert(entry1, entry2);

        const long index = m_keys(entry1, entry2);

        return m_values[size_t(index)];
      }

      /// \brief Getter to access to a given value referenced by the input key without prior check.
      T & get(const IndexPair & key)
      {
        assert(this->exist(key));
        const Index entry1 = key.first;
        const Index entry2 = key.second;
        const long index = m_keys(entry1, entry2);

        return m_values[size_t(index)];
      }

      /// \brief Getter to access to a given value referenced by the input key without prior check.
      const T & get(const IndexPair & key) const
      {
        assert(this->exist(key));
        const Index entry1 = key.first;
        const Index entry2 = key.second;
        const long index = m_keys(entry1, entry2);

        return m_values[size_t(index)];
      }

#ifdef PINOCCHIO_WITH_CXX23_SUPPORT
      T & operator[](const Index entry1, const Index entry2)
      {
        return this->operator[]({entry1, entry2});
      }
#endif

      iterator begin()
      {
        return m_values.begin();
      }

      iterator end()
      {
        return m_values.end();
      }

      iterator rbegin()
      {
        return m_values.cbegin();
      }

      iterator rend()
      {
        return m_values.cend();
      }

      /// \brief Increase the capacity of the vector (the total number of elements that the vector
      /// can hold without requiring reallocation) to a value that's greater or equal to new_cap
      void reserve(const size_t new_cap)
      {
        m_values.reserve(new_cap);
      }

    protected:
      Array m_keys;
      Vector m_values;
    };

  } // namespace container

} // namespace pinocchio

#endif // ifndef __pinocchio_container_double_entry_container_hpp__
