//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_collision_parallel_broadphase_hpp__
#define __pinocchio_collision_parallel_broadphase_hpp__

#include "pinocchio/collision/pool/broadphase-manager.hpp"
#include "pinocchio/collision/broadphase.hpp"
#include "pinocchio/algorithm/geometry.hpp"
#include "pinocchio/utils/openmp.hpp"

#include <cstdint>

namespace pinocchio
{

  template<
    typename BroadPhaseManagerDerived,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorPool,
    typename CollisionVectorResult>
  void computeCollisionsInParallel(
    const size_t num_threads,
    BroadPhaseManagerPoolBase<BroadPhaseManagerDerived, Scalar, Options, JointCollectionTpl> & pool,
    const Eigen::MatrixBase<ConfigVectorPool> & q,
    const Eigen::MatrixBase<CollisionVectorResult> & res,
    const bool stopAtFirstCollisionInConfiguration = false,
    const bool stopAtFirstCollisionInBatch = false)
  {
    typedef BroadPhaseManagerPoolBase<BroadPhaseManagerDerived, Scalar, Options, JointCollectionTpl>
      Pool;
    typedef typename Pool::Model Model;
    typedef typename Pool::Data Data;
    typedef typename Pool::ModelVector ModelVector;
    typedef typename Pool::DataVector DataVector;
    typedef typename Pool::BroadPhaseManager BroadPhaseManager;
    typedef typename Pool::BroadPhaseManagerVector BroadPhaseManagerVector;

    const ModelVector & models = pool.getModels();
    const Model & model_check = models[0];
    DataVector & datas = pool.getDatas();
    BroadPhaseManagerVector & broadphase_managers = pool.getBroadPhaseManagers();
    CollisionVectorResult & res_ = res.const_cast_derived();

    PINOCCHIO_CHECK_INPUT_ARGUMENT(num_threads <= pool.size(), "The pool is too small");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(q.rows(), model_check.nq);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(q.cols(), res.size());
    res_.fill(false);

    setDefaultOpenMPSettings(num_threads);
    const Eigen::DenseIndex batch_size = res.size();

    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(ConfigVectorPool) ConfigVectorPoolPlain;
    PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(ConfigVectorPoolPlain) q_thread(num_threads, q);

    OpenMPException openmp_exception;

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

        const Model & model = models[(size_t)thread_id];
        Data & data = datas[(size_t)thread_id];
        BroadPhaseManager & manager = broadphase_managers[(size_t)thread_id];
        const ConfigVectorPoolPlain & q = q_thread[size_t(thread_id)];

        openmp_exception.run([=, &is_colliding, &model, &data, &manager, &q, &res_] {
          // lambda start corpus

          res_[i] =
            computeCollisions(model, data, manager, q.col(i), stopAtFirstCollisionInConfiguration);

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
        BroadPhaseManager & manager = broadphase_managers[(size_t)thread_id];
        const ConfigVectorPoolPlain & q = q_thread[size_t(thread_id)];

        openmp_exception.run([=, &model, &data, &manager, &q, &res_] {
          // lambda start corpus
          res_[i] =
            computeCollisions(model, data, manager, q.col(i), stopAtFirstCollisionInConfiguration);
          // lambda end corpus
        });
      }
    }

    openmp_exception.rethrowException();
  }

  ///
  /// \brief Evaluate the collision over a set of trajectories and return whether a trajectory
  /// contains a collision
  ///
  template<
    typename BroadPhaseManagerDerived,
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl>
  void computeCollisionsInParallel(
    const size_t num_threads,
    BroadPhaseManagerPoolBase<BroadPhaseManagerDerived, Scalar, Options, JointCollectionTpl> & pool,
    const std::vector<Eigen::MatrixXd> & trajectories,
    std::vector<VectorXb> & res,
    const bool stopAtFirstCollisionInTrajectory = false)
  {
    typedef BroadPhaseManagerPoolBase<BroadPhaseManagerDerived, Scalar, Options, JointCollectionTpl>
      Pool;
    typedef typename Pool::Model Model;
    typedef typename Pool::Data Data;
    typedef typename Pool::ModelVector ModelVector;
    typedef typename Pool::DataVector DataVector;
    typedef typename Pool::BroadPhaseManager BroadPhaseManager;
    typedef typename Pool::BroadPhaseManagerVector BroadPhaseManagerVector;

    const ModelVector & models = pool.getModels();
    const Model & model_check = models[0];
    DataVector & datas = pool.getDatas();
    BroadPhaseManagerVector & broadphase_managers = pool.getBroadPhaseManagers();

    PINOCCHIO_CHECK_INPUT_ARGUMENT(num_threads <= pool.size(), "The pool is too small");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(trajectories.size(), res.size());

    for (size_t k = 0; k < trajectories.size(); ++k)
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(trajectories[k].cols(), res[k].size());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(trajectories[k].rows(), model_check.nq);
    }

    setDefaultOpenMPSettings(num_threads);
    const size_t batch_size = trajectories.size();

    OpenMPException openmp_exception;

    size_t i = 0;
#pragma omp parallel for schedule(static)
    for (i = 0; i < batch_size; i++)
    {
      const int thread_id = omp_get_thread_num();

      const Model & model = models[size_t(thread_id)];
      Data & data = datas[(size_t)thread_id];
      const Eigen::MatrixXd & current_traj = trajectories[i];
      VectorXb & res_current_traj = res[i];
      res_current_traj.fill(false);
      BroadPhaseManager & manager = broadphase_managers[size_t(thread_id)];

      openmp_exception.run([=, &model, &data, &manager, &current_traj, &res_current_traj] {
        // lambda start corpus
        for (Eigen::DenseIndex col_id = 0; col_id < current_traj.cols(); ++col_id)
        {
          res_current_traj[col_id] =
            computeCollisions(model, data, manager, current_traj.col(col_id), true);

          if (res_current_traj[col_id] && stopAtFirstCollisionInTrajectory)
            break;
        }
        // lambda end corpus
      });
    }

    openmp_exception.rethrowException();
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_collision_parallel_broadphase_hpp__
