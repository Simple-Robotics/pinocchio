//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_visitors_hxx__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_visitors_hxx__

namespace pinocchio
{

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodySystemsTplSolveInPlaceBackwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodySystemsTplSolveInPlaceBackwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
    typedef typename DelassusOperator::CustomData CustomData;

    typedef boost::fusion::vector<const Model &, const Data &, CustomData &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      const pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      const Data & data,
      CustomData & custom_data)
    {
      typedef typename Data::Force Force;
      typedef typename Data::Motion Motion;
      typedef typename Motion::Vector6 Vector6;
      typedef typename Data::Matrix6 Matrix6;
      typedef std::pair<JointIndex, JointIndex> JointPair;

      const auto & neighbours = data.neighbour_links;
      auto & joint_cross_coupling = data.joint_cross_coupling;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      const auto Jcols = jmodel.jointCols(data.J);

      Force & ofi = custom_data.of_augmented[i];

      jmodel.jointVelocitySelector(custom_data.u).noalias() -= Jcols.transpose() * ofi.toVector();

      const auto res = (jdata.Dinv() * jmodel.jointVelocitySelector(custom_data.u))
                         .eval(); // Abuse of notation to reuse existing unused data variable

      if (neighbours[i].size())
      {
        const Vector6 a_tmp = Jcols * res;

        for (JointIndex vertex_j : neighbours[i])
        {
          const Matrix6 & crosscoupling_ij =
            (i > vertex_j) ? joint_cross_coupling.get(JointPair(vertex_j, i))
                           : joint_cross_coupling.get(JointPair(i, vertex_j)).transpose();
          custom_data.of_augmented[vertex_j].toVector().noalias() += crosscoupling_ij * a_tmp;
        }
      }

      if (parent > 0)
      {
        ofi.toVector().noalias() += jdata.UDinv() * jmodel.jointVelocitySelector(custom_data.u);
        custom_data.of_augmented[parent] += ofi;
      }
    }
  };

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodySystemsTplSolveInPlaceForwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodySystemsTplSolveInPlaceForwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
    typedef typename DelassusOperator::CustomData CustomData;
    typedef std::pair<JointIndex, JointIndex> JointPair;

    typedef boost::fusion::vector<const Model &, const Data &, CustomData &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      const pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      const Data & data,
      CustomData & custom_data)
    {
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Matrix6x Matrix6x;
      typedef typename Data::Matrix6 Matrix6;

      const auto J_cols = jmodel.jointCols(data.J);

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];
      const auto & joint_neighbours = data.neighbour_links[i];

      auto & oai = custom_data.oa_augmented[i];
      oai = custom_data.oa_augmented[parent];

      if (joint_neighbours.size())
      {
        Force coupling_forces = Force::Zero();
        for (const JointIndex vertex_j : joint_neighbours)
        {
          const Matrix6 & edge_ij =
            (i > vertex_j) ? data.joint_cross_coupling.get(JointPair(vertex_j, i)).transpose()
                           : data.joint_cross_coupling.get(JointPair(i, vertex_j));

          coupling_forces.toVector().noalias() +=
            edge_ij * custom_data.oa_augmented[vertex_j].toVector();
        }

        jmodel.jointVelocitySelector(custom_data.u).noalias() -=
          J_cols.transpose() * coupling_forces.toVector();
      }

      // Abuse of notation using custom_data.ddq for storing delta ddq
      jmodel.jointVelocitySelector(custom_data.ddq).noalias() =
        jdata.Dinv() * jmodel.jointVelocitySelector(custom_data.u)
        - jdata.UDinv().transpose() * oai.toVector();
      oai.toVector().noalias() += J_cols * jmodel.jointVelocitySelector(custom_data.ddq);
    }
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_visitors_hxx__
