//
// Copyright (c) 2015-2022 CNRS INRIA
//

#ifndef __pinocchio_multibody_geometry_hxx__
#define __pinocchio_multibody_geometry_hxx__

#include <algorithm>

#include "pinocchio/multibody/model.hpp"

#if BOOST_VERSION / 100 % 1000 >= 60
  #include <boost/bind/bind.hpp>
  #include <boost/utility.hpp>
#else
  #include <boost/bind.hpp>
#endif
#include <boost/foreach.hpp>

#include <sstream>
/// @cond DEV

namespace pinocchio
{

  inline CollisionPair::CollisionPair()
  : Base((std::numeric_limits<GeomIndex>::max)(), (std::numeric_limits<GeomIndex>::max)())
  {
  }

  ///
  /// \brief Default constructor of a collision pair from two collision object indexes.
  /// \remarks The two indexes must be different, otherwise the constructor throws.
  ///
  /// \param[in] co1 Index of the first collision object.
  /// \param[in] co2 Index of the second collision object.
  ///
  inline CollisionPair::CollisionPair(const GeomIndex co1, const GeomIndex co2)
  : Base(co1, co2)
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(co1 != co2, "The index of collision objects must not be equal.");
  }

  inline bool CollisionPair::operator==(const CollisionPair & rhs) const
  {
    return (first == rhs.first && second == rhs.second)
           || (first == rhs.second && second == rhs.first);
  }

  inline bool CollisionPair::operator!=(const CollisionPair & other) const
  {
    return !(*this == other);
  }

  inline void CollisionPair::disp(std::ostream & os) const
  {
    os << "collision pair (" << first << "," << second << ")\n";
  }

  inline std::ostream & operator<<(std::ostream & os, const CollisionPair & X)
  {
    X.disp(os);
    return os;
  }

  inline GeometryModel GeometryModel::clone() const
  {
    GeometryModel res;
    res.ngeoms = ngeoms;
    res.collisionPairs = collisionPairs;
    res.collisionPairMapping = collisionPairMapping;

    res.geometryObjects.reserve(geometryObjects.size());
    for (const GeometryObject & geometry_object : geometryObjects)
    {
      res.geometryObjects.push_back(geometry_object.clone());
    }

    return res;
  }

  inline GeometryData::GeometryData(const GeometryModel & geom_model)
  : oMg(geom_model.ngeoms)
  , activeCollisionPairs(geom_model.collisionPairs.size(), true)
#ifdef PINOCCHIO_WITH_HPP_FCL
  , distanceRequests(geom_model.collisionPairs.size(), hpp::fcl::DistanceRequest(true))
  , distanceResults(geom_model.collisionPairs.size())
  , collisionRequests(
      geom_model.collisionPairs.size(), hpp::fcl::CollisionRequest(::hpp::fcl::NO_REQUEST, 1))
  , collisionResults(geom_model.collisionPairs.size())
  , contactPatchRequests(geom_model.collisionPairs.size()) // use default constructor
  , contactPatchResults(geom_model.collisionPairs.size())  // use default constructor
  , radius()
  , collisionPairIndex(0)
#endif // PINOCCHIO_WITH_HPP_FCL
  , innerObjects()
  , outerObjects()
  {
#ifdef PINOCCHIO_WITH_HPP_FCL
    BOOST_FOREACH (hpp::fcl::CollisionRequest & creq, collisionRequests)
    {
      creq.enable_cached_gjk_guess = true;
    }
  #if HPP_FCL_VERSION_AT_LEAST(1, 4, 5)
    BOOST_FOREACH (hpp::fcl::DistanceRequest & dreq, distanceRequests)
    {
      dreq.enable_cached_gjk_guess = true;
    }
  #endif
    collision_functors.reserve(geom_model.collisionPairs.size());
    contact_patch_functors.reserve(geom_model.collisionPairs.size());
    distance_functors.reserve(geom_model.collisionPairs.size());

    for (size_t cp_index = 0; cp_index < geom_model.collisionPairs.size(); ++cp_index)
    {
      const CollisionPair & cp = geom_model.collisionPairs[cp_index];
      const GeometryObject & obj_1 = geom_model.geometryObjects[cp.first];
      const GeometryObject & obj_2 = geom_model.geometryObjects[cp.second];

      collision_functors.push_back(ComputeCollision(obj_1, obj_2));
      contact_patch_functors.push_back(ComputeContactPatch(obj_1, obj_2));
      distance_functors.push_back(ComputeDistance(obj_1, obj_2));
    }
#endif
    fillInnerOuterObjectMaps(geom_model);
  }

  inline GeometryData::GeometryData(const GeometryData & other)
  : oMg(other.oMg)
  , activeCollisionPairs(other.activeCollisionPairs)
#ifdef PINOCCHIO_WITH_HPP_FCL
  , distanceRequests(other.distanceRequests)
  , distanceResults(other.distanceResults)
  , collisionRequests(other.collisionRequests)
  , collisionResults(other.collisionResults)
  , contactPatchRequests(other.contactPatchRequests)
  , contactPatchResults(other.contactPatchResults)
  , radius(other.radius)
  , collisionPairIndex(other.collisionPairIndex)
  , collision_functors(other.collision_functors)
  , contact_patch_functors(other.contact_patch_functors)
  , distance_functors(other.distance_functors)
#endif // PINOCCHIO_WITH_HPP_FCL
  , innerObjects(other.innerObjects)
  , outerObjects(other.outerObjects)
  {
  }

  inline GeometryData & GeometryData::operator=(const GeometryData & other)
  {
    if (this != &other)
    {
      oMg = other.oMg;
      activeCollisionPairs = other.activeCollisionPairs;
#ifdef PINOCCHIO_WITH_HPP_FCL
      distanceRequests = other.distanceRequests;
      distanceResults = other.distanceResults;
      collisionRequests = other.collisionRequests;
      collisionResults = other.collisionResults;
      contactPatchRequests = other.contactPatchRequests;
      contactPatchResults = other.contactPatchResults;
      radius = other.radius;
      collisionPairIndex = other.collisionPairIndex;
      collision_functors = other.collision_functors;
      contact_patch_functors = other.contact_patch_functors;
      distance_functors = other.distance_functors;
#endif // PINOCCHIO_WITH_HPP_FCL
      innerObjects = other.innerObjects;
      outerObjects = other.outerObjects;
    }
    return *this;
  }

  inline GeometryData::~GeometryData()
  {
  }

  template<typename S2, int O2, template<typename, int> class JointCollectionTpl>
  GeomIndex GeometryModel::addGeometryObject(
    const GeometryObject & object, const ModelTpl<S2, O2, JointCollectionTpl> & model)
  {
    if (object.parentFrame < (FrameIndex)model.nframes)
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        model.frames[object.parentFrame].parentJoint == object.parentJoint,
        "The object joint parent and its frame joint parent do not match.");

    Eigen::DenseIndex idx = (Eigen::DenseIndex)(ngeoms++);
    geometryObjects.push_back(object);
    geometryObjects.back().parentJoint = model.frames[object.parentFrame].parentJoint;

    collisionPairMapping.conservativeResize(idx + 1, idx + 1);
    collisionPairMapping.col(idx).fill(-1);
    collisionPairMapping.row(idx).fill(-1);
    assert(collisionPairMapping.cols() == (Eigen::DenseIndex)ngeoms);

    return (GeomIndex)idx;
  }

  inline GeomIndex GeometryModel::addGeometryObject(const GeometryObject & object)
  {
    Eigen::DenseIndex idx = (Eigen::DenseIndex)(ngeoms++);
    geometryObjects.push_back(object);

    collisionPairMapping.conservativeResize(idx + 1, idx + 1);
    collisionPairMapping.col(idx).fill(-1);
    collisionPairMapping.row(idx).fill(-1);
    assert(collisionPairMapping.cols() == (Eigen::DenseIndex)ngeoms);

    return (GeomIndex)idx;
  }

  inline void GeometryModel::removeGeometryObject(const std::string & name)
  {
    GeomIndex i = 0;
    GeometryObjectVector::iterator it;
    for (it = geometryObjects.begin(); it != geometryObjects.end(); ++it, ++i)
    {
      if (it->name == name)
      {
        break;
      }
    }
    PINOCCHIO_THROW_IF(
      (it == geometryObjects.end()), std::invalid_argument,
      (std::string("Object ") + name + std::string(" does not belong to model")).c_str());
    // Remove all collision pairs that contain i as first or second index,
    for (CollisionPairVector::iterator itCol = collisionPairs.begin();
         itCol != collisionPairs.end(); ++itCol)
    {
      if ((itCol->first == i) || (itCol->second == i))
      {
        itCol = collisionPairs.erase(itCol);
        itCol--;
      }
      else
      {
        // Indices of objects after the one that is removed should be decreased by one.
        if (itCol->first > i)
          itCol->first--;
        if (itCol->second > i)
          itCol->second--;
      }
    }
    geometryObjects.erase(it);
    ngeoms--;
  }

  inline GeomIndex GeometryModel::getGeometryId(const std::string & name) const
  {
#if BOOST_VERSION / 100 % 1000 >= 60
    using namespace boost::placeholders;
#endif
    GeometryObjectVector::const_iterator it = std::find_if(
      geometryObjects.begin(), geometryObjects.end(),
      boost::bind(&GeometryObject::name, _1) == name);
    return GeomIndex(it - geometryObjects.begin());
  }

  inline bool GeometryModel::existGeometryName(const std::string & name) const
  {
#if BOOST_VERSION / 100 % 1000 >= 60
    using namespace boost::placeholders;
#endif
    return std::find_if(
             geometryObjects.begin(), geometryObjects.end(),
             boost::bind(&GeometryObject::name, _1) == name)
           != geometryObjects.end();
  }

  inline void GeometryData::fillInnerOuterObjectMaps(const GeometryModel & geomModel)
  {
    innerObjects.clear();
    outerObjects.clear();

    for (GeomIndex gid = 0; gid < geomModel.geometryObjects.size(); gid++)
      innerObjects[geomModel.geometryObjects[gid].parentJoint].push_back(gid);

    BOOST_FOREACH (const CollisionPair & pair, geomModel.collisionPairs)
    {
      outerObjects[geomModel.geometryObjects[pair.first].parentJoint].push_back(pair.second);
    }
  }

  inline std::ostream & operator<<(std::ostream & os, const GeometryModel & geomModel)
  {
    os << "Nb geometry objects = " << geomModel.ngeoms << std::endl;

    for (GeomIndex i = 0; i < (GeomIndex)(geomModel.ngeoms); ++i)
    {
      os << geomModel.geometryObjects[i] << std::endl;
    }

    return os;
  }

  inline std::ostream & operator<<(std::ostream & os, const GeometryData & geomData)
  {
#ifdef PINOCCHIO_WITH_HPP_FCL
    os << "Number of collision pairs = " << geomData.activeCollisionPairs.size() << std::endl;

    for (PairIndex i = 0; i < (PairIndex)(geomData.activeCollisionPairs.size()); ++i)
    {
      os << "Pairs " << i << (geomData.activeCollisionPairs[i] ? " active" : " inactive")
         << std::endl;
    }
#else
    os << "WARNING** Without fcl library, no collision checking or distance computations are "
          "possible. Only geometry placements can be computed."
       << std::endl;
#endif
    os << "Number of geometry objects = " << geomData.oMg.size() << std::endl;

    return os;
  }

  inline void GeometryModel::addCollisionPair(const CollisionPair & pair)
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      pair.first < ngeoms, "The input pair.first is larger than the number of geometries contained "
                           "in the GeometryModel");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      pair.second < ngeoms, "The input pair.second is larger than the number of geometries "
                            "contained in the GeometryModel");
    if (!existCollisionPair(pair))
    {
      collisionPairs.push_back(pair);

      collisionPairMapping((Eigen::DenseIndex)pair.second, (Eigen::DenseIndex)pair.first) =
        (int)(collisionPairs.size() - 1);
      collisionPairMapping((Eigen::DenseIndex)pair.first, (Eigen::DenseIndex)pair.second) =
        collisionPairMapping(
          (Eigen::DenseIndex)pair.second,
          (Eigen::DenseIndex)pair.first); // make symmetric
    }
  }

  inline void GeometryModel::setCollisionPairs(const MatrixXb & map, const bool upper)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      map.rows(), (Eigen::DenseIndex)ngeoms, "Input map does not have the correct number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      map.cols(), (Eigen::DenseIndex)ngeoms,
      "Input map does not have the correct number of columns.");
    removeAllCollisionPairs();

    for (Eigen::DenseIndex i = 0; i < (Eigen::DenseIndex)ngeoms; ++i)
    {
      for (Eigen::DenseIndex j = i + 1; j < (Eigen::DenseIndex)ngeoms; ++j)
      {
        if (upper)
        {
          if (map(i, j))
            addCollisionPair(CollisionPair((std::size_t)i, (std::size_t)j));
        }
        else
        {
          if (map(j, i))
            addCollisionPair(CollisionPair((std::size_t)i, (std::size_t)j));
        }
      }
    }
  }

  inline void GeometryModel::addAllCollisionPairs()
  {
    removeAllCollisionPairs();
    for (GeomIndex i = 0; i < ngeoms; ++i)
    {
      const JointIndex joint_i = geometryObjects[i].parentJoint;
      for (GeomIndex j = i + 1; j < ngeoms; ++j)
      {
        const JointIndex joint_j = geometryObjects[j].parentJoint;
        if (joint_i != joint_j)
          addCollisionPair(CollisionPair(i, j));
      }
    }
  }

  inline void GeometryModel::removeCollisionPair(const CollisionPair & pair)
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      pair.first < ngeoms, "The input pair.first is larger than the number of geometries contained "
                           "in the GeometryModel");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      pair.second < ngeoms, "The input pair.second is larger than the number of geometries "
                            "contained in the GeometryModel");

    const long index = (long)findCollisionPair(pair);

    if (index != (long)collisionPairs.size())
    {
      collisionPairMapping((Eigen::DenseIndex)pair.first, (Eigen::DenseIndex)pair.second) =
        collisionPairMapping((Eigen::DenseIndex)pair.second, (Eigen::DenseIndex)pair.first) = -1;
      collisionPairs.erase(collisionPairs.begin() + index);

      for (Eigen::DenseIndex i = 0; i < collisionPairMapping.cols(); ++i)
      {
        for (Eigen::DenseIndex j = i + 1; j < collisionPairMapping.cols(); ++j)
        {
          if (collisionPairMapping(i, j) > index)
          {
            collisionPairMapping(i, j)--;
            collisionPairMapping(j, i) = collisionPairMapping(i, j);
          }
        }
      }
    }
  }

  inline void GeometryModel::removeAllCollisionPairs()
  {
    collisionPairs.clear();
    collisionPairMapping.fill(-1);
  }

  inline bool GeometryModel::existCollisionPair(const CollisionPair & pair) const
  {
    return collisionPairMapping((Eigen::DenseIndex)pair.first, (Eigen::DenseIndex)pair.second)
           != -1;
  }

  inline PairIndex GeometryModel::findCollisionPair(const CollisionPair & pair) const
  {
    int res = collisionPairMapping((Eigen::DenseIndex)pair.first, (Eigen::DenseIndex)pair.second);
    if (res == -1)
      return collisionPairs.size();
    else
      return (PairIndex)res;
  }

  inline void GeometryData::activateCollisionPair(const PairIndex pair_id)
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      pair_id < activeCollisionPairs.size(),
      "The input argument pair_id is larger than the number of collision pairs contained in "
      "activeCollisionPairs.");
    activeCollisionPairs[pair_id] = true;
  }

  inline void GeometryData::activateAllCollisionPairs()
  {
    std::fill(activeCollisionPairs.begin(), activeCollisionPairs.end(), true);
  }

  inline void GeometryData::setActiveCollisionPairs(
    const GeometryModel & geom_model, const MatrixXb & map, const bool upper)
  {
    const Eigen::DenseIndex ngeoms = (Eigen::DenseIndex)geom_model.ngeoms;
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      map.rows(), ngeoms, "Input map does not have the correct number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      map.cols(), ngeoms, "Input map does not have the correct number of columns.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      geom_model.collisionPairs.size(), activeCollisionPairs.size(),
      "Current geometry data and the input geometry model are not conistent.");

    for (size_t k = 0; k < geom_model.collisionPairs.size(); ++k)
    {
      const CollisionPair & cp = geom_model.collisionPairs[k];

      Eigen::DenseIndex i, j;
      if (upper)
      {
        j = (Eigen::DenseIndex)std::max(cp.first, cp.second);
        i = (Eigen::DenseIndex)std::min(cp.first, cp.second);
      }
      else
      {
        i = (Eigen::DenseIndex)std::max(cp.first, cp.second);
        j = (Eigen::DenseIndex)std::min(cp.first, cp.second);
      }

      activeCollisionPairs[k] = map(i, j);
    }
  }

  inline void GeometryData::setGeometryCollisionStatus(
    const GeometryModel & geom_model, const GeomIndex geom_id, const bool enable_collision)
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      geom_id < geom_model.ngeoms, "The index of the geometry is not valid");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      geom_model.collisionPairs.size(), activeCollisionPairs.size(),
      "Current geometry data and the input geometry model are not conistent.");

    for (size_t k = 0; k < geom_model.collisionPairs.size(); ++k)
    {
      const CollisionPair & cp = geom_model.collisionPairs[k];
      if (cp.first == geom_id || cp.second == geom_id)
      {
        activeCollisionPairs[k] = enable_collision;
      }
    }
  }

#ifdef PINOCCHIO_WITH_HPP_FCL
  inline void GeometryData::setSecurityMargins(
    const GeometryModel & geom_model,
    const MatrixXs & security_margin_map,
    const bool upper,
    const bool sync_distance_upper_bound)
  {
    const Eigen::DenseIndex ngeoms = (Eigen::DenseIndex)geom_model.ngeoms;
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      security_margin_map.rows(), ngeoms, "Input map does not have the correct number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      security_margin_map.cols(), ngeoms, "Input map does not have the correct number of columns.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      geom_model.collisionPairs.size(), collisionRequests.size(),
      "Current geometry data and the input geometry model are not consistent.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      geom_model.collisionPairs.size(), contactPatchRequests.size(),
      "Current geometry data and the input geometry model are not consistent.");

    for (size_t k = 0; k < geom_model.collisionPairs.size(); ++k)
    {
      const CollisionPair & cp = geom_model.collisionPairs[k];

      Eigen::DenseIndex i, j;
      if (upper)
      {
        j = (Eigen::DenseIndex)std::max(cp.first, cp.second);
        i = (Eigen::DenseIndex)std::min(cp.first, cp.second);
      }
      else
      {
        i = (Eigen::DenseIndex)std::max(cp.first, cp.second);
        j = (Eigen::DenseIndex)std::min(cp.first, cp.second);
      }

      collisionRequests[k].security_margin = security_margin_map(i, j);
      if (sync_distance_upper_bound)
        collisionRequests[k].distance_upper_bound = collisionRequests[k].security_margin;
    }
  }
#endif // ifdef PINOCCHIO_WITH_HPP_FCL

  inline void GeometryData::deactivateCollisionPair(const PairIndex pair_id)
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      pair_id < activeCollisionPairs.size(),
      "The input argument pair_id is larger than the number of collision pairs contained in "
      "activeCollisionPairs.");
    activeCollisionPairs[pair_id] = false;
  }

  inline void GeometryData::deactivateAllCollisionPairs()
  {
    std::fill(activeCollisionPairs.begin(), activeCollisionPairs.end(), false);
  }

} // namespace pinocchio

/// @endcond

#endif // ifndef __pinocchio_multibody_geometry_hxx__
