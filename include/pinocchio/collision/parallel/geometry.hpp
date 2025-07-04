//
// Copyright (c) 2021-2024 INRIA
//

#ifndef __pinocchio_collision_parallel_geometry_hpp__
#define __pinocchio_collision_parallel_geometry_hpp__

#include "pinocchio/multibody/pool/geometry.hpp"
#include "pinocchio/algorithm/geometry.hpp"
#include "pinocchio/collision/collision.hpp"
#include "pinocchio/utils/openmp.hpp"

#include <cstdint>

namespace pinocchio
{

  inline bool computeCollisionsInParallel(
    const size_t num_threads,
    const GeometryModel & geom_model,
    GeometryData & geom_data,
    const bool stopAtFirstCollision = false)
  {
    bool is_colliding = false;

    setDefaultOpenMPSettings(num_threads);
    std::size_t cp_index = 0;

    OpenMPException openmp_exception;

#pragma omp parallel for schedule(dynamic)
    for (cp_index = 0; cp_index < geom_model.collisionPairs.size(); ++cp_index)
    {
      if (stopAtFirstCollision && is_colliding)
        continue;

      openmp_exception.run([=, &is_colliding, &geom_model, &geom_data] {
        // lambda start corpus
        const CollisionPair & collision_pair = geom_model.collisionPairs[cp_index];

        if (
          geom_data.activeCollisionPairs[cp_index]
          && !(
            geom_model.geometryObjects[collision_pair.first].disableCollision
            || geom_model.geometryObjects[collision_pair.second].disableCollision))
        {
          bool res = computeCollision(geom_model, geom_data, cp_index);
          if (!is_colliding && res)
          {
            is_colliding = true;
            geom_data.collisionPairIndex = cp_index; // first pair to be in collision
          }
        }
        // lambda end corpus
      });
    }

    openmp_exception.rethrowException();

    return is_colliding;
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType>
  bool computeCollisionsInParallel(
    const size_t num_threads,
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const GeometryModel & geom_model,
    GeometryData & geom_data,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const bool stopAtFirstCollision = false)
  {
    updateGeometryPlacements(model, data, geom_model, geom_data, q);
    return computeCollisionsInParallel(num_threads, geom_model, geom_data, stopAtFirstCollision);
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorPool,
    typename CollisionVectorResult>
  void computeCollisionsInParallel(
    const size_t num_threads,
    GeometryPoolTpl<Scalar, Options, JointCollectionTpl> & pool,
    const Eigen::MatrixBase<ConfigVectorPool> & q,
    const Eigen::MatrixBase<CollisionVectorResult> & res,
    const bool stopAtFirstCollisionInConfiguration = false,
    const bool stopAtFirstCollisionInBatch = false)
  {
    typedef GeometryPoolTpl<Scalar, Options, JointCollectionTpl> Pool;
    typedef typename Pool::Model Model;
    typedef typename Pool::Data Data;
    typedef typename Pool::GeometryModel GeometryModel;
    typedef typename Pool::GeometryData GeometryData;
    typedef typename Pool::ModelVector ModelVector;
    typedef typename Pool::DataVector DataVector;
    typedef typename Pool::GeometryModelVector GeometryModelVector;
    typedef typename Pool::GeometryDataVector GeometryDataVector;

    PINOCCHIO_CHECK_INPUT_ARGUMENT(pool.size() > 0, "The pool should have at least one element");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(num_threads <= pool.size(), "The pool is too small");

    const ModelVector & models = pool.getModels();
    const Model & model_check = models[0];
    const GeometryModelVector & geometry_models = pool.getGeometryModels();
    DataVector & datas = pool.getDatas();
    GeometryDataVector & geometry_datas = pool.getGeometryDatas();
    CollisionVectorResult & res_ = res.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(q.rows(), model_check.nq);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(q.cols(), res.size());
    res_.fill(false);

    setDefaultOpenMPSettings(num_threads);
    const Eigen::DenseIndex batch_size = res.size();

    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(ConfigVectorPool) ConfigVectorPoolPlain;
    PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(ConfigVectorPoolPlain) q_thread(num_threads, q);

    OpenMPException openmp_exception;

    // TODO(jcarpent): set one res_ per thread to enhance efficiency
    if (stopAtFirstCollisionInBatch)
    {
      bool is_colliding = false;
      Eigen::DenseIndex i = 0;
#pragma omp parallel for schedule(static)
      for (i = 0; i < batch_size; i++)
      {
        if (is_colliding)
          continue;

        const int thread_id = omp_get_thread_num();

        const Model & model = models[size_t(thread_id)];
        Data & data = datas[size_t(thread_id)];
        const GeometryModel & geometry_model = geometry_models[size_t(thread_id)];
        GeometryData & geometry_data = geometry_datas[size_t(thread_id)];
        const ConfigVectorPoolPlain & q = q_thread[size_t(thread_id)];

        openmp_exception.run(
          [=, &is_colliding, &model, &data, &geometry_model, &geometry_data, &q, &res_] {
            // lambda start corpus

            res_[i] = computeCollisions(
              model, data, geometry_model, geometry_data, q.col(i),
              stopAtFirstCollisionInConfiguration);

            if (res_[i])
            {
              is_colliding = true;
            }
            // lambda end corpus
          });
      }
    }
    else
    {
      Eigen::DenseIndex i = 0;
#pragma omp parallel for schedule(static)
      for (i = 0; i < batch_size; i++)
      {
        const int thread_id = omp_get_thread_num();

        const Model & model = models[(size_t)thread_id];
        Data & data = datas[(size_t)thread_id];
        const GeometryModel & geometry_model = geometry_models[(size_t)thread_id];
        GeometryData & geometry_data = geometry_datas[(size_t)thread_id];
        const ConfigVectorPoolPlain & q = q_thread[size_t(thread_id)];

        openmp_exception.run([=, &model, &data, &geometry_model, &geometry_data, &q, &res_] {
          // lambda start corpus
          res_[i] = computeCollisions(
            model, data, geometry_model, geometry_data, q.col(i),
            stopAtFirstCollisionInConfiguration);
          // lambda end corpus
        });
      }
    }

    openmp_exception.rethrowException();
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_collision_parallel_geometry_hpp__
