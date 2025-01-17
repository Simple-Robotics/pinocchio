//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_loop_constrained_aba_hxx__
#define __pinocchio_algorithm_loop_constrained_aba_hxx__

#include <algorithm>

/// @cond DEV

namespace pinocchio
{

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class Allocator>
  inline void initLcaba(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<RigidConstraintModelTpl<Scalar, Options>, Allocator> & contact_models)
  {

    typedef typename Model::JointIndex JointIndex;
    typedef std::pair<JointIndex, JointIndex> JointPair;
    typedef Data::Matrix6 Matrix6;

    // Ensure only LOCAL_WORLD_ALIGNED constraints are accepted
    for (std::size_t i = 0; i < contact_models.size(); ++i)
    {
      const RigidConstraintModelTpl<Scalar, Options> & contact_model = contact_models[i];
      switch (contact_model.reference_frame)
      {
      case LOCAL_WORLD_ALIGNED:
        break;
      default:
        assert(false && "Frames other than LOCAL_WORLD_ALIGNED not accepted");
        break;
      }
    }

    auto & neighbours = data.neighbour_links;
    neighbours.resize(static_cast<size_t>(model.njoints));

    // Get links supporting constraints
    std::fill(data.constraints_supported_dim.begin(), data.constraints_supported_dim.end(), 0);
    for (std::size_t i = 0; i < contact_models.size(); ++i)
    {
      const RigidConstraintModelTpl<Scalar, Options> & contact_model = contact_models[i];
      const JointIndex joint1_id = contact_model.joint1_id;
      const JointIndex joint2_id = contact_model.joint2_id;

      const auto constraint_size = contact_model.size();
      data.constraints_supported_dim[joint1_id] += constraint_size;
      data.constraints_supported_dim[joint2_id] += constraint_size;

      if (joint2_id > 0)
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
      const auto leaf_with_least_neighbors_it =
        ::std::min_element(leaf_vertices.begin(), leaf_vertices.end());
      const JointIndex leaf_with_least_neighbors = *leaf_with_least_neighbors_it;

      const JointIndex joint_id = leaf_with_least_neighbors;
      leaf_vertices.remove(joint_id);
      elimination_order.push_back(joint_id);

      const JointIndex parent_id = model.parents[joint_id];
      num_children[parent_id]--;
      // If the number of children joints of parent is reaching zero, this means that parent is now
      // a leaf node.
      if (num_children[parent_id] == 0 && parent_id != 0)
        leaf_vertices.push_back(parent_id);

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

        // Remove joint_id form the list of neighbours for neighbour_j_neighbours
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

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentVectorType>
  struct LCABAForwardStep1
  : public fusion::JointUnaryVisitorBase<
      LCABAForwardStep1<Scalar, Options, JointCollectionTpl, ConfigVectorType, TangentVectorType>>
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;

    typedef boost::fusion::
      vector<const Model &, Data &, const ConfigVectorType &, const TangentVectorType &>
        ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      Data & data,
      const Eigen::MatrixBase<ConfigVectorType> & q,
      const Eigen::MatrixBase<TangentVectorType> & v)
    {
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Motion Motion;

      const JointIndex i = jmodel.id();
      Motion & ov = data.ov[i];
      jmodel.calc(jdata.derived(), q.derived(), v.derived());

      const JointIndex parent = model.parents[i];
      data.liMi[i] = model.jointPlacements[i] * jdata.M();
      if (parent > 0)
        data.oMi[i] = data.oMi[parent] * data.liMi[i];
      else
        data.oMi[i] = data.liMi[i];

      jmodel.jointCols(data.J) = data.oMi[i].act(jdata.S());

      ov = data.oMi[i].act(jdata.v());
      if (parent > 0)
        ov += data.ov[parent];

      data.oa_gf[i] = data.oMi[i].act(jdata.c());
      if (parent > 0)
        data.oa_gf[i] += (data.ov[parent] ^ ov);

      data.oinertias[i] = data.oYcrb[i] = data.oMi[i].act(model.inertias[i]);
      data.oYaba[i] = data.oYcrb[i].matrix();

      data.oh[i] = data.oYcrb[i] * ov; // necessary for ABA derivatives
      data.of[i] = ov.cross(data.oh[i]);
    }
  };

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  struct LCABABackwardStep
  : public fusion::JointUnaryVisitorBase<LCABABackwardStep<Scalar, Options, JointCollectionTpl>>
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;

    typedef boost::fusion::vector<const Model &, Data &> ArgsType;

    template<typename JointModel>
    static void algo(
      const JointModelBase<JointModel> & jmodel,
      JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      Data & data)
    {

      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Inertia Inertia;
      typedef typename Data::Force Force;
      typedef typename Data::Matrix6 Matrix6;

      const auto & neighbours = data.neighbour_links;

      typedef std::pair<JointIndex, JointIndex> JointPair;

      auto & joint_cross_coupling = data.joint_cross_coupling;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];
      typename Inertia::Matrix6 & Ia = data.oYaba[i];

      typedef
        typename SizeDepType<JointModel::NV>::template ColsReturn<pinocchio::Data::Matrix6x>::Type
          ColBlock;
      ColBlock Jcols = jmodel.jointCols(data.J);

      Force & fi = data.of[i];

      jmodel.jointVelocitySelector(data.u).noalias() -= Jcols.transpose() * fi.toVector();

      jdata.U().noalias() = Ia * Jcols;
      jdata.StU().noalias() = Jcols.transpose() * jdata.U();

      // Account for the rotor inertia contribution
      jdata.StU().diagonal() += jmodel.jointVelocitySelector(model.armature);

      pinocchio::internal::PerformStYSInversion<Scalar>::run(jdata.StU(), jdata.Dinv());
      auto & JDinv = jdata.UDinv();
      JDinv.noalias() = Jcols * jdata.Dinv();

      data.oL[i].setIdentity();
      data.oL[i].noalias() -= JDinv * jdata.U().transpose();

      Motion a_tmp;
      a_tmp.toVector().noalias() = data.oL[i] * data.oa_gf[i].toVector();
      a_tmp.toVector().noalias() += JDinv * jmodel.jointVelocitySelector(data.u);

      const auto & joint_neighbours = neighbours[i];
      for (size_t j = 0; j < joint_neighbours.size(); j++)
      {
        const JointIndex vertex_j = joint_neighbours[j];
        const Matrix6 & crosscoupling_ij =
          (i > vertex_j)
            ? joint_cross_coupling[JointPair(vertex_j, i)]
            : joint_cross_coupling[JointPair(i, vertex_j)].transpose(); // avoid memalloc

        auto & crosscoupling_ix_Jcols = jdata.UDinv();
        crosscoupling_ix_Jcols.noalias() =
          crosscoupling_ij * Jcols; // Warning: UDinv() is actually edge_ij * J

        auto & crosscoupling_ij_Jcols_Dinv = jdata.U();
        crosscoupling_ij_Jcols_Dinv.noalias() = crosscoupling_ix_Jcols * jdata.Dinv();

        data.oYaba[vertex_j].noalias() -=
          crosscoupling_ij_Jcols_Dinv
          * crosscoupling_ix_Jcols.transpose(); // Warning: UDinv() is actually edge_ij * J, U() is
                                                // actually edge_ij * J_cols * Dinv
        data.of[vertex_j].toVector().noalias() += crosscoupling_ij * a_tmp.toVector();

        const Matrix6 crosscoupling_ij_oL = crosscoupling_ij * data.oL[i];
        if (vertex_j == parent)
        {
          data.oYaba[parent].noalias() += crosscoupling_ij_oL + crosscoupling_ij_oL.transpose();
        }
        else
        {
          if (vertex_j < parent)
          {
            joint_cross_coupling[{vertex_j, parent}].noalias() += crosscoupling_ij_oL;
          }
          else
          {
            joint_cross_coupling[{parent, vertex_j}].noalias() += crosscoupling_ij_oL.transpose();
          }
        }

        for (size_t k = j + 1; k < joint_neighbours.size(); ++k)
        {
          const JointIndex vertex_k = joint_neighbours[k];

          const Matrix6 & edge_ik = (i > vertex_k) ? edges[JointPair(vertex_k, i)]
                                                   : edges[JointPair(i, vertex_k)].transpose();

          crosscoupling_ix_Jcols.noalias() = edge_ik * Jcols;

          assert(vertex_j != vertex_k && "Must never happen!");
          if (vertex_j < vertex_k)
          {
            joint_cross_coupling[{vertex_j, vertex_k}].noalias() -=
              crosscoupling_ij_Jcols_Dinv
              * crosscoupling_ix_Jcols.transpose(); // Warning: UDinv() is actually edge_ik * J_col,
                                                    // U() is edge_ij * J_col * Dinv
          }
          else // if (vertex_k < vertex_j)
          {
            joint_cross_coupling[{vertex_k, vertex_j}].transpose().noalias() -=
              crosscoupling_ij_Jcols_Dinv
              * crosscoupling_ix_Jcols.transpose(); // Warning: UDinv() is actually edge_ik *
                                                    // J_col, U() is edge_ij * J_col * Dinv
          }
        }
      }

      jdata.U().noalias() = Ia * Jcols;
      jdata.UDinv().noalias() =
        jdata.U() * jdata.Dinv(); // TODO:check where its used when parent == 0
      if (parent > 0)
      {
        Ia.noalias() -= jdata.UDinv() * jdata.U().transpose();
        fi.toVector().noalias() +=
          Ia * data.oa_gf[i].toVector() + jdata.UDinv() * jmodel.jointVelocitySelector(data.u);
        data.oYaba[parent] += Ia;
        data.of[parent] += fi;
      }
    }
  };

  // A reduced backward sweep that only propagates the affine terms
  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  struct LCABABackwardStepReduced
  : public fusion::JointUnaryVisitorBase<
      LCABABackwardStepReduced<Scalar, Options, JointCollectionTpl>>
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;

    typedef boost::fusion::vector<const Model &, Data &> ArgsType;

    template<typename JointModel>
    static void algo(
      const JointModelBase<JointModel> & jmodel,
      JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      Data & data)
    {

      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Force Force;
      typedef typename Data::Motion Motion;
      typedef typename Motion::Vector6 Vector6;
      typedef typename Data::Matrix6 Matrix6;

      const auto & neighbours = data.neighbour_links;

      typedef std::pair<JointIndex, JointIndex> JointPair;

      auto & joint_cross_coupling = data.joint_cross_coupling;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      typedef
        typename SizeDepType<JointModel::NV>::template ColsReturn<pinocchio::Data::Matrix6x>::Type
          ColBlock;
      const ColBlock & Jcols = jmodel.jointCols(data.J);

      Force & fi = data.of[i];

      jmodel.jointVelocitySelector(data.u).noalias() = -Jcols.transpose() * fi.toVector();

      jmodel.jointVelocitySelector(data.g).noalias() =
        jdata.Dinv()
        * jmodel.jointVelocitySelector(
          data.u); // Abuse of notation to reuse existing unused data variable

      const Vector6 a_tmp = Jcols * jmodel.jointVelocitySelector(data.g);

      for (JointIndex vertex_j : neighbours[i])
      {
        const Matrix6 & edge_ij = (i > vertex_j) ? edges[JointPair(vertex_j, i)]
                                                 : edges[JointPair(i, vertex_j)].transpose();
        data.of[vertex_j].toVector().noalias() += edge_ij * a_tmp;
      }

      if (parent > 0)
      {
        data.of[parent].toVector().noalias() +=
          fi.toVector() + jdata.U() * jmodel.jointVelocitySelector(data.g);
      }
    }
  };

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  struct LCABAForwardStep2
  : public fusion::JointUnaryVisitorBase<LCABAForwardStep2<Scalar, Options, JointCollectionTpl>>
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;

    typedef std::pair<JointIndex, JointIndex> JointPair;
    typedef typename Data::Matrix6 Matrix6;

    typedef boost::fusion::vector<const Model &, Data &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      Data & data)
    {
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Matrix6x Matrix6x;

      typedef typename SizeDepType<JointModel::NV>::template ColsReturn<Matrix6x>::Type ColBlock;
      ColBlock J_cols = jmodel.jointCols(data.J);

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];
      const std::vector<JointIndex> & neighbours = data.neighbour_links[i];

      data.oa_gf[i] += data.oa_gf[parent]; // does take into account the gravity field

      Force coupling_forces = Force::Zero();
      for (JointIndex vertex_j : neighbours)
      {
        const Matrix6 & edge_ij = (i > vertex_j)
                                    ? data.joint_cross_coupling[JointPair(vertex_j, i)].transpose()
                                    : data.joint_cross_coupling[JointPair(i, vertex_j)];
        coupling_forces.toVector().noalias() += edge_ij * data.oa_gf[vertex_j].toVector();
      }

      jmodel.jointVelocitySelector(data.u).noalias() -=
        J_cols.transpose() * coupling_forces.toVector();

      jmodel.jointVelocitySelector(data.ddq).noalias() =
        jdata.Dinv() * jmodel.jointVelocitySelector(data.u)
        - jdata.UDinv().transpose() * data.oa_gf[i].toVector();
      data.oa_gf[i].toVector().noalias() += J_cols * jmodel.jointVelocitySelector(data.ddq);

      // Handle consistent output
      data.oa[i] = data.oa_gf[i]; // + model.gravity;
      // data.of[i] = data.oinertias[i] * data.oa_gf[i] + data.ov[i].cross(data.oh[i]);
    }
  };

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  struct LCABAReducedForwardStep
  : public fusion::JointUnaryVisitorBase<
      LCABAReducedForwardStep<Scalar, Options, JointCollectionTpl>>
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;

    typedef std::pair<JointIndex, JointIndex> JointPair;
    typedef typename Data::Matrix6 Matrix6;

    typedef boost::fusion::vector<const Model &, Data &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      Data & data)
    {
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Matrix6x Matrix6x;

      typedef typename SizeDepType<JointModel::NV>::template ColsReturn<Matrix6x>::Type ColBlock;
      ColBlock J_cols = jmodel.jointCols(data.J);

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];
      const auto & neighbours = data.neighbour_links[i];

      data.oa_gf[i] = data.oa_gf[parent]; // does take into account the gravity field

      Force & fi = data.of[i];
      for (JointIndex vertex_j : neighbours)
      {

        const Matrix6 & edge_ij = (i > vertex_j)
                                    ? data.joint_cross_coupling[JointPair(vertex_j, i)].transpose()
                                    : data.joint_cross_coupling[JointPair(i, vertex_j)];
        fi.toVector().noalias() += edge_ij * data.oa_gf[vertex_j].toVector();
      }

      jmodel.jointVelocitySelector(data.u).noalias() = -J_cols.transpose() * fi.toVector();

      // Abuse of notation using data.g for storing delta ddq
      jmodel.jointVelocitySelector(data.g).noalias() =
        jdata.Dinv()
        * (jmodel.jointVelocitySelector(data.u) - jdata.U().transpose() * data.oa_gf[i].toVector());
      data.oa_gf[i].toVector().noalias() += J_cols * jmodel.jointVelocitySelector(data.g);

      data.oa[i] += data.oa_gf[i];
    }
  };

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentVectorType1,
    typename TangentVectorType2,
    class ContactModelAllocator,
    class ContactDataAllocator>
  inline const typename DataTpl<Scalar, Options, JointCollectionTpl>::TangentVectorType & lcaba(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType1> & v,
    const Eigen::MatrixBase<TangentVectorType2> & tau,
    const std::vector<RigidConstraintModelTpl<Scalar, Options>, ContactModelAllocator> &
      contact_models,
    std::vector<RigidConstraintDataTpl<Scalar, Options>, ContactDataAllocator> & contact_datas,
    ProximalSettingsTpl<Scalar> & settings)
  {

    assert(model.check(data) && "data is not consistent with model.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The joint configuration vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      v.size(), model.nv, "The joint velocity vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      tau.size(), model.nv, "The joint torque vector is not of right size");

    typedef typename ModelTpl<Scalar, Options, JointCollectionTpl>::JointIndex JointIndex;
    typedef std::pair<JointIndex, JointIndex> JointPair;
    typedef Data::Matrix6 Matrix6;
    typedef RigidConstraintModel::Matrix36 Matrix36;

    data.u = tau;
    data.oa_gf[0] = -model.gravity;
    data.oa[0] = data.oa_gf[0];
    data.of[0].setZero();

    const Scalar mu = Scalar(1) / settings.mu;

    for (auto & coupling : data.joint_cross_coupling)
      coupling.setZero();

    typedef LCABAForwardStep1<
      Scalar, Options, JointCollectionTpl, ConfigVectorType, TangentVectorType1>
      Pass1;
    for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    {
      Pass1::run(
        model.joints[i], data.joints[i],
        typename Pass1::ArgsType(model, data, q.derived(), v.derived()));
    }

    for (std::size_t i = 0; i < contact_models.size(); ++i)
    {
      RigidConstraintData & cdata = contact_datas[i];
      cdata.contact_force.setZero();
      const RigidConstraintModelTpl<Scalar, Options> & contact_model = contact_models[i];
      typename RigidConstraintData::Motion & vc1 = contact_datas[i].contact1_velocity;
      typename RigidConstraintData::Motion & vc2 = contact_datas[i].contact2_velocity;
      const JointIndex joint_id = contact_model.joint1_id;
      const JointIndex joint2_id = contact_model.joint2_id;

      const typename RigidConstraintModel::BaumgarteCorrectorParameters & corrector =
        contact_model.corrector;
      typename RigidConstraintData::Motion & contact_vel_err = cdata.contact_velocity_error;

      SE3 & oMc1 = cdata.oMc1;
      oMc1 = data.oMi[joint_id] * contact_model.joint1_placement;

      SE3 & oMc2 = cdata.oMc2;
      oMc2 = data.oMi[joint2_id] * contact_model.joint2_placement;

      SE3 & c1Mc2 = cdata.c1Mc2;
      c1Mc2 = oMc1.actInv(oMc2);

      vc1 = oMc1.actInv(data.ov[joint_id]);
      if (joint2_id > 0)
        vc2 = oMc2.actInv(data.ov[joint2_id]);
      else
        vc2.setZero();
      const Motion vc2_in_frame1 = c1Mc2.act(vc2);

      if (contact_model.type == CONTACT_6D)
      {
        cdata.contact_placement_error = -log6(c1Mc2);
        contact_vel_err = vc1 - vc2_in_frame1;
        const Matrix6 A1 = oMc1.toActionMatrixInverse();
        const Matrix6 A1tA1 = A1.transpose() * A1;
        data.oYaba[joint_id].noalias() += mu * A1tA1;

        // Baumgarte
        if (check_expression_if_real<Scalar, false>(
              isZero(corrector.Kp, static_cast<Scalar>(0.))
              && isZero(corrector.Kd, static_cast<Scalar>(0.))))
        {
          cdata.contact_acceleration_desired.setZero();
        }
        else
        {
          cdata.contact_acceleration_desired.toVector().noalias() =
            -(corrector.Kd.asDiagonal() * contact_vel_err.toVector())
            - (corrector.Kp.asDiagonal() * cdata.contact_placement_error.toVector());
        }

        if (joint2_id > 0)
        {
          const Matrix6 A2 = -A1; // only for 6D case. also used below for computing A2tA2 and A1tA2
          data.oYaba[joint2_id].noalias() += mu * A1tA1;
          data.of[joint2_id].toVector().noalias() +=
            A2.transpose()
            * (cdata.contact_force.toVector() - mu * cdata.contact_acceleration_desired.toVector());

          const JointPair jp =
            joint_id < joint2_id ? JointPair{joint_id, joint2_id} : JointPair{joint2_id, joint_id};
          data.joint_cross_coupling[jp] -= mu * A1tA1;
        }
        else
        {
          cdata.contact_acceleration_desired.toVector().noalias() -= A1 * model.gravity.toVector();
        }

        data.of[joint_id].toVector().noalias() +=
          A1.transpose()
          * (cdata.contact_force.toVector() - mu * cdata.contact_acceleration_desired.toVector());
      }
      else if (contact_model.type == CONTACT_3D)
      {
        const Matrix36 & A1 = oMc1.toActionMatrixInverse().topRows<3>();
        data.oYaba[joint_id].noalias() += mu * A1.transpose() * A1;

        if (check_expression_if_real<Scalar, false>(
              isZero(corrector.Kp, static_cast<Scalar>(0.))
              && isZero(corrector.Kd, static_cast<Scalar>(0.))))
        {
          cdata.contact_acceleration_desired.setZero();
        }
        else
        {
          cdata.contact_acceleration_desired.linear().noalias() =
            -(corrector.Kd.asDiagonal() * contact_vel_err.linear())
            - (corrector.Kp.asDiagonal() * cdata.contact_placement_error.linear());
          cdata.contact_acceleration_desired.angular().setZero();
        }

        cdata.contact_acceleration_desired.linear().noalias() -= vc1.angular().cross(vc1.linear());
        if (joint2_id > 0)
        {
          const Matrix36 A2 =
            -c1Mc2.rotation() * (oMc2.toActionMatrixInverse().topRows<3>()); // TODO:remove memalloc

          cdata.contact_acceleration_desired.linear().noalias() +=
            c1Mc2.rotation() * vc2.angular().cross(vc2.linear());
          data.oYaba[joint2_id].noalias() += mu * A2.transpose() * A2;
          data.of[joint2_id].toVector().noalias() +=
            A2.transpose()
            * (cdata.contact_force.linear() - mu * cdata.contact_acceleration_desired.linear());

          if (joint_id < joint2_id)
          {
            data.joint_cross_coupling[{joint_id, joint2_id}].noalias() += mu * A1.transpose() * A2;
          }
          else
          {
            data.joint_cross_coupling[{joint2_id, joint_id}].noalias() += mu * A2.transpose() * A1;
          }
        }
        else
        {
          cdata.contact_acceleration_desired.linear().noalias() -= A1 * model.gravity.toVector();
        }

        data.of[joint_id].toVector().noalias() +=
          A1.transpose()
          * (cdata.contact_force.linear() - mu * cdata.contact_acceleration_desired.linear());
      }
      else
      {
        assert(false && "Must never happen");
      }
    }

    typedef LCABABackwardStep<Scalar, Options, JointCollectionTpl> Pass2;

    const std::vector<JointIndex> elim_order = data.elimination_order;

    for (JointIndex i : elim_order)
    {
      Pass2::run(model.joints[i], data.joints[i], typename Pass2::ArgsType(model, data));
    }

    assert(settings.mu > 0 && "constrainedABA requires mu > 0.");

    typedef LCABAForwardStep2<Scalar, Options, JointCollectionTpl> Pass3;

    for (int it = int(elim_order.size()) - 1; it >= 0; --it)
    {
      const JointIndex i = elim_order[size_t(it)];
      if (data.constraints_supported_dim[i] > 0)
        Pass3::run(model.joints[i], data.joints[i], typename Pass3::ArgsType(model, data));
    }

    typedef LCABABackwardStepReduced<Scalar, Options, JointCollectionTpl> ReducedPass2;
    typedef LCABAReducedForwardStep<Scalar, Options, JointCollectionTpl> ReducedPass3;
    data.g.setZero();
    int iter = 0;
    for (iter = 1; iter < settings.max_iter; iter++)
    {
      settings.absolute_residual = 0.0;
      for (JointIndex j = 1; j < (JointIndex)model.njoints; ++j)
      {
        if (data.constraints_supported_dim[j] > 0)
          data.of[j].setZero();
      }
      for (std::size_t j = 0; j < contact_models.size(); ++j)
      {
        const RigidConstraintModelTpl<Scalar, Options> & contact_model = contact_models[j];
        RigidConstraintData & cdata = contact_datas[j];
        const JointIndex joint1_id = contact_model.joint1_id;
        const JointIndex joint2_id = contact_model.joint2_id;
        typename RigidConstraintData::Motion & contact_acc_err = cdata.contact_acceleration_error;

        if (contact_model.type == CONTACT_6D)
        {
          if (joint2_id > 0)
            contact_acc_err = cdata.oMc1.actInv((data.oa[joint1_id] - data.oa[joint2_id]))
                              - cdata.contact_acceleration_desired;
          else
            contact_acc_err =
              cdata.oMc1.actInv((data.oa[joint1_id])) - cdata.contact_acceleration_desired;

          cdata.contact_force.toVector().noalias() += mu * contact_acc_err.toVector();

          data.of[joint1_id] += cdata.oMc1.act(Force(mu * contact_acc_err.toVector()));

          if (joint2_id > 0)
            data.of[joint2_id] -= cdata.oMc1.act(Force(mu * contact_acc_err.toVector()));
        }
        else
        {
          if (joint2_id > 0)
            contact_acc_err.linear() =
              cdata.oMc1.actInv(data.oa[joint1_id]).linear()
              - cdata.c1Mc2.rotation() * cdata.oMc2.actInv(data.oa[joint2_id]).linear()
              - cdata.contact_acceleration_desired.linear();
          else
            contact_acc_err.linear() = cdata.oMc1.actInv(data.oa[joint1_id]).linear()
                                       - cdata.contact_acceleration_desired.linear();

          cdata.contact_force.linear() += mu * contact_acc_err.linear();

          data.of[joint1_id] += cdata.oMc1.act(Force(mu * contact_acc_err.toVector()));

          if (joint2_id > 0)
            data.of[joint2_id].toVector().noalias() -=
              cdata.oMc2.toActionMatrixInverse().topRows<3>().transpose()
              * (cdata.c1Mc2.rotation().transpose() * (mu * contact_acc_err.linear()));
        }
        Scalar c_err_residual = contact_acc_err.toVector().template lpNorm<Eigen::Infinity>();
        if (settings.absolute_residual < c_err_residual)
          settings.absolute_residual = c_err_residual;
      }
      if (settings.absolute_residual < settings.absolute_accuracy)
        break;

      // reduced backward sweep
      for (JointIndex j : elim_order)
      {
        if (data.constraints_supported_dim[j] > 0)
          ReducedPass2::run(
            model.joints[j], data.joints[j], typename ReducedPass2::ArgsType(model, data));
      }

      // forward sweep
      data.oa_gf[0].setZero();
      for (int it = int(elim_order.size()) - 1; it >= 0; it--)
      {
        const JointIndex j = elim_order[size_t(it)];
        if (data.constraints_supported_dim[j] > 0)
        {
          ReducedPass3::run(
            model.joints[j], data.joints[j], typename ReducedPass3::ArgsType(model, data));
        }
      }
      data.ddq += data.g;
      // related to least-squares residual
      if (data.g.template lpNorm<Eigen::Infinity>() < settings.absolute_accuracy)
        break;
    }
    settings.iter = iter;

    // outward sweep after convergence/timeout for joints not supporting a constraint
    data.oa_gf[0] = -model.gravity;
    for (int it = int(elim_order.size()) - 1; it >= 0; --it)
    {
      const JointIndex j = elim_order[size_t(it)];
      if (data.constraints_supported_dim[j] == 0)
        Pass3::run(model.joints[j], data.joints[j], typename Pass3::ArgsType(model, data));
      else
        data.oa_gf[j] = data.oa[j];
    }

    return data.ddq;
  }

} // namespace pinocchio

/// @endcond

#endif // ifndef __pinocchio_algorithm_loop_constrained_aba_hxx__
