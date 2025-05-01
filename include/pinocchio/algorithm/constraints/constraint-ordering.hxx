//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_ordering_hxx__
#define __pinocchio_algorithm_constraints_constraint_ordering_hxx__

#include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hpp"

/// @cond DEV

namespace pinocchio
{

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  struct MinimalOrderingConstraintStepVisitor
  : visitors::ConstraintUnaryVisitorBase<
      MinimalOrderingConstraintStepVisitor<Scalar, Options, JointCollectionTpl>>
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;

    typedef boost::fusion::vector<const Model &, Data &> ArgsType;

    typedef visitors::ConstraintUnaryVisitorBase<
      MinimalOrderingConstraintStepVisitor<Scalar, Options, JointCollectionTpl>>
      Base;

    template<typename ConstraintModel>
    static void
    algo(const ConstraintModelBase<ConstraintModel> & cmodel, const Model & model, Data & data)
    {
      algo_step(cmodel.derived(), model, data);
    }

    template<typename ConstraintModel>
    static void algo_step(
      const BinaryConstraintModelBase<ConstraintModel> & cmodel, const Model & model, Data & data)
    {
      typedef std::pair<JointIndex, JointIndex> JointPair;
      typedef typename Data::Matrix6 Matrix6;

      PINOCCHIO_UNUSED_VARIABLE(model);

      const JointIndex joint1_id = cmodel.joint1_id;
      const JointIndex joint2_id = cmodel.joint2_id;
      auto & neighbours = data.neighbour_links;

      const auto constraint_size = cmodel.size();
      data.constraints_supported_dim[joint1_id] += constraint_size;
      data.constraints_supported_dim[joint2_id] += constraint_size;

      if (joint2_id > 0 && joint1_id > 0)
      {
        const JointPair joint_pair =
          joint1_id > joint2_id ? JointPair{joint2_id, joint1_id} : JointPair{joint1_id, joint2_id};

        if (!data.joint_cross_coupling.exist(joint_pair))
          data.joint_cross_coupling[joint_pair] = Matrix6::Zero();

        auto & joint1_neighbours = neighbours[joint1_id];
        auto & joint2_neighbours = neighbours[joint2_id];

        if (
          std::find(joint1_neighbours.begin(), joint1_neighbours.end(), joint2_id)
          == joint1_neighbours.end())
          joint1_neighbours.push_back(joint2_id);
        if (
          std::find(joint2_neighbours.begin(), joint2_neighbours.end(), joint1_id)
          == joint2_neighbours.end())
          joint2_neighbours.push_back(joint1_id);
      }
    }
    template<typename _Scalar, int _Options>
    static void algo_step(
      const FrictionalJointConstraintModelTpl<_Scalar, _Options> & cmodel,
      const Model & model,
      Data & data)
    {
      for (const JointIndex joint_id : cmodel.getActiveJoints())
      {
        data.constraints_supported_dim[joint_id] += model.nvs[joint_id];
      }
    }

    using Base::run;

    template<typename ConstraintModel>
    static void run(
      const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
      const Model & model,
      Data & data)
    {
      algo(cmodel.derived(), model, data);
    }

    template<
      typename _Scalar,
      int _Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    static void run(
      const pinocchio::ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl> & cmodel,
      const Model & model,
      Data & data)
    {
      ArgsType args(model, data);
      run(cmodel.derived(), args);
    }
  }; // struct MinimalOrderingConstraintStepVisitor

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator>
  inline void computeJointMinimalOrdering(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models)
  {
    typedef typename Model::JointIndex JointIndex;
    typedef std::pair<JointIndex, JointIndex> JointPair;
    typedef Data::Matrix6 Matrix6;

    // First step: for each joint, collect their neighbourds
    auto & neighbours = data.neighbour_links;
    for (auto & neighbour_elt : neighbours)
      neighbour_elt.clear();
    data.joint_cross_coupling.clear();

    // Get links supporting constraints
    std::fill(data.constraints_supported_dim.begin(), data.constraints_supported_dim.end(), 0);
    typedef MinimalOrderingConstraintStepVisitor<Scalar, Options, JointCollectionTpl> Step;
    for (std::size_t i = 0; i < constraint_models.size(); ++i)
    {
      const ConstraintModel & cmodel = constraint_models[i];
      Step::run(cmodel, model, data);
    }

    // Second step: order the joints according to the minimum degree heuristic
    auto & elimination_order = data.elimination_order;
    elimination_order.clear(); // clearing in case inited once more

    std::vector<size_t> num_children(size_t(model.njoints), 0);
    std::list<JointIndex> leaf_vertices;

    for (JointIndex joint_id = JointIndex(model.njoints - 1); joint_id > 0; --joint_id)
    {
      num_children[joint_id] = model.children[joint_id].size();
      if (num_children[joint_id] == 0)
        leaf_vertices.push_back(joint_id);
    }

    while (leaf_vertices.size() > 0)
    {
      JointIndex joint_id_with_least_neighbors = std::numeric_limits<JointIndex>::max();
      ;
      size_t least_neighbours = std::numeric_limits<size_t>::max();

      for (const auto joint_id : leaf_vertices)
      {
        assert(joint_id != 0);
        const auto leaf_num_neighours = neighbours[joint_id].size();
        if (leaf_num_neighours < least_neighbours)
        {
          least_neighbours = leaf_num_neighours;
          joint_id_with_least_neighbors = joint_id;
        }
      }

      const JointIndex joint_id = joint_id_with_least_neighbors;
      assert(joint_id != 0);
      leaf_vertices.remove(joint_id);
      elimination_order.push_back(joint_id);

      const JointIndex parent_id = model.parents[joint_id];
      num_children[parent_id]--;
      // If the number of children joints of parent is reaching zero, this means that parent is now
      // a leaf node.
      if (num_children[parent_id] == 0 && parent_id != 0)
        leaf_vertices.push_front(parent_id);

      data.constraints_supported_dim[parent_id] += data.constraints_supported_dim[joint_id];
      const auto & joint_neighbours = neighbours[joint_id];
      auto & parent_neighbours = neighbours[parent_id];
      for (size_t j = 0; j < joint_neighbours.size(); j++)
      {
        const JointIndex neighbour_j = joint_neighbours[j];
        auto & neighbour_j_neighbours = neighbours[neighbour_j];
        if (neighbour_j != parent_id)
        {
          const JointPair jp_pair = neighbour_j < parent_id ? JointPair(neighbour_j, parent_id)
                                                            : JointPair(parent_id, neighbour_j);

          if (!data.joint_cross_coupling.exist(jp_pair))
          {
            data.joint_cross_coupling[jp_pair] = Matrix6::Zero(); // add joint_cross_coupling

            if (
              std::find(parent_neighbours.begin(), parent_neighbours.end(), neighbour_j)
              == parent_neighbours.end())
            {
              parent_neighbours.push_back(neighbour_j);
              neighbour_j_neighbours.push_back(parent_id);
            }
          }
        }

        // Remove joint_id from the list of neighbours for neighbour_j_neighbours
        neighbour_j_neighbours.erase(
          std::remove(neighbour_j_neighbours.begin(), neighbour_j_neighbours.end(), joint_id),
          neighbour_j_neighbours.end());

        for (size_t k = j + 1; k < joint_neighbours.size(); ++k)
        {
          const JointIndex neighbour_k = joint_neighbours[k];
          auto & neighbour_k_neighbours = neighbours[neighbour_k];
          assert(neighbour_k != neighbour_j && "Must never happen!");
          const JointPair cross_coupling = neighbour_j < neighbour_k
                                             ? JointPair{neighbour_j, neighbour_k}
                                             : JointPair{neighbour_k, neighbour_j};

          if (!data.joint_cross_coupling.exist(cross_coupling))
          {
            data.joint_cross_coupling[cross_coupling] = Matrix6::Zero(); // add edge
            neighbour_j_neighbours.push_back(neighbour_k);
            neighbour_k_neighbours.push_back(neighbour_j);
          }
        }
      }
    }
  }
} // namespace pinocchio

/// @endcond

#endif // ifndef __pinocchio_algorithm_constraints_constraint_ordering_hxx__
