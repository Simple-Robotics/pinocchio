//
// Copyright (c) 2017-2021 CNRS INRIA
//

#ifndef __pinocchio_serialization_serializable_hpp__
#define __pinocchio_serialization_serializable_hpp__

#include "pinocchio/serialization/archive.hpp"

namespace pinocchio
{
  namespace serialization
  {

    template<class Derived>
    struct Serializable
    {
    private:
      Derived & serializableDerived()
      {
        return *static_cast<Derived *>(this);
      }
      const Derived & serializableDerived() const
      {
        return *static_cast<const Derived *>(this);
      }

    public:
      /// \brief Loads a Derived object from a text file.
      void loadFromText(const std::string & filename)
      {
        pinocchio::serialization::loadFromText(serializableDerived(), filename);
      }

      /// \brief Saves a Derived object as a text file.
      void saveToText(const std::string & filename) const
      {
        pinocchio::serialization::saveToText(serializableDerived(), filename);
      }

      /// \brief Loads a Derived object from a stream string.
      void loadFromStringStream(std::istringstream & is)
      {
        pinocchio::serialization::loadFromStringStream(serializableDerived(), is);
      }

      /// \brief Saves a Derived object to a string stream.
      void saveToStringStream(std::stringstream & ss) const
      {
        pinocchio::serialization::saveToStringStream(serializableDerived(), ss);
      }

      /// \brief Loads a Derived object from a  string.
      void loadFromString(const std::string & str)
      {
        pinocchio::serialization::loadFromString(serializableDerived(), str);
      }

      /// \brief Saves a Derived object to a string.
      std::string saveToString() const
      {
        return pinocchio::serialization::saveToString(serializableDerived());
      }

      /// \brief Loads a Derived object from an XML file.
      void loadFromXML(const std::string & filename, const std::string & tag_name)
      {
        pinocchio::serialization::loadFromXML(serializableDerived(), filename, tag_name);
      }

      /// \brief Saves a Derived object as an XML file.
      void saveToXML(const std::string & filename, const std::string & tag_name) const
      {
        pinocchio::serialization::saveToXML(serializableDerived(), filename, tag_name);
      }

      /// \brief Loads a Derived object from an binary file.
      void loadFromBinary(const std::string & filename)
      {
        pinocchio::serialization::loadFromBinary(serializableDerived(), filename);
      }

      /// \brief Saves a Derived object as an binary file.
      void saveToBinary(const std::string & filename) const
      {
        pinocchio::serialization::saveToBinary(serializableDerived(), filename);
      }

      /// \brief Loads a Derived object from a binary container.
      void loadFromBinary(boost::asio::streambuf & container)
      {
        pinocchio::serialization::loadFromBinary(serializableDerived(), container);
      }

      /// \brief Saves a Derived object as a binary container.
      void saveToBinary(boost::asio::streambuf & container) const
      {
        pinocchio::serialization::saveToBinary(serializableDerived(), container);
      }

      /// \brief Loads a Derived object from a static binary container.
      void loadFromBinary(StaticBuffer & container)
      {
        pinocchio::serialization::loadFromBinary(serializableDerived(), container);
      }

      /// \brief Saves a Derived object as a static binary container.
      void saveToBinary(StaticBuffer & container) const
      {
        pinocchio::serialization::saveToBinary(serializableDerived(), container);
      }
    };

  } // namespace serialization
} // namespace pinocchio

#endif // ifndef __pinocchio_serialization_serializable_hpp__
