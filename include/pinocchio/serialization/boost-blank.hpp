//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_serialization_boost_blank_hpp__
#define __pinocchio_serialization_boost_blank_hpp__

namespace boost
{
  namespace serialization
  {

    template<typename Archive>
    void serialize(Archive &, ::boost::blank &, const unsigned int)
    {
      // do nothing
    }
  } // namespace serialization
} // namespace boost

#endif // __pinocchio_serialization_boost_blank_hpp__
