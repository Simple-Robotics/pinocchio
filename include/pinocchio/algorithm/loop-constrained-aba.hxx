//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_loop_constrained_aba_hxx__
#define __pinocchio_algorithm_loop_constrained_aba_hxx__

/// @cond DEV

namespace pinocchio
{

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
      data.oYaba_augmented[i] = data.oYcrb[i].matrix();

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
      typedef typename JointModel::JointDataDerived JointData;
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Force Force;
      typedef typename Data::Vector6 Vector6;
      typedef typename Data::Matrix6 Matrix6;

      typedef std::pair<JointIndex, JointIndex> JointPair;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];
      auto & Ia = data.oYaba_augmented[i];

      const auto Jcols = jmodel.jointCols(data.J);

      Force & fi = data.of[i];

      jmodel.jointVelocitySelector(data.u).noalias() -= Jcols.transpose() * fi.toVector();

      jdata.U().noalias() = Ia * Jcols;
      jdata.StU().noalias() = Jcols.transpose() * jdata.U();

      // Account for the rotor inertia contribution
      jdata.StU().diagonal() += jmodel.jointVelocitySelector(model.armature);

      pinocchio::internal::PerformStYSInversion<Scalar>::run(jdata.StU(), jdata.Dinv());

      jdata.UDinv().noalias() =
        jdata.U() * jdata.Dinv(); // TODO:check where its used when parent == 0
      if (parent > 0)
      {
        Ia.noalias() -= jdata.UDinv() * jdata.U().transpose();
        data.oYaba_augmented[parent] += Ia;

        fi.toVector().noalias() +=
          Ia * data.oa_gf[i].toVector() + jdata.UDinv() * jmodel.jointVelocitySelector(data.u);
        data.of[parent] += fi;
      }

      // End of the classic ABA backward pass - beginning of cross-coupling handling
      const auto & neighbours = data.neighbour_links;
      auto & joint_cross_coupling = data.joint_cross_coupling;
      const auto & joint_neighbours = neighbours[i];

      if (joint_neighbours.size() == 0)
        return; // We can return from this point as this joint has no neighbours

      using Matrix6xNV = typename std::remove_reference<typename JointData::UDTypeRef>::type;
      typedef Eigen::Map<Matrix6xNV> MapMatrix6xNV;
      MapMatrix6xNV mat1_tmp = MapMatrix6xNV(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, 6, jmodel.nv()));
      MapMatrix6xNV mat2_tmp = MapMatrix6xNV(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, 6, jmodel.nv()));

      auto & JDinv = mat1_tmp;
      JDinv.noalias() = Jcols * jdata.Dinv();

      // oL == data.oL[i]
      Matrix6 oL = -JDinv * jdata.U().transpose();
      oL += Matrix6::Identity();

      // a_tmp is a Spatial Acceleration
      Vector6 a_tmp = oL * data.oa_gf[i].toVector();
      a_tmp.noalias() += JDinv * jmodel.jointVelocitySelector(data.u);

      for (size_t j = 0; j < joint_neighbours.size(); j++)
      {
        const JointIndex vertex_j = joint_neighbours[j];
        const Matrix6 & crosscoupling_ij =
          (i > vertex_j)
            ? joint_cross_coupling.get(JointPair(vertex_j, i))
            : joint_cross_coupling.get(JointPair(i, vertex_j)).transpose(); // avoid memalloc

        auto & crosscoupling_ix_Jcols = mat1_tmp;
        crosscoupling_ix_Jcols.noalias() =
          crosscoupling_ij * Jcols; // Warning: UDinv() is actually edge_ij * J

        auto & crosscoupling_ij_Jcols_Dinv = mat2_tmp;
        crosscoupling_ij_Jcols_Dinv.noalias() = crosscoupling_ix_Jcols * jdata.Dinv();

        data.oYaba_augmented[vertex_j].noalias() -=
          crosscoupling_ij_Jcols_Dinv
          * crosscoupling_ix_Jcols.transpose(); // Warning: UDinv() is actually edge_ij * J, U() is
                                                // actually edge_ij * J_cols * Dinv
        data.of[vertex_j].toVector().noalias() += crosscoupling_ij * a_tmp;

        const Matrix6 crosscoupling_ij_oL = crosscoupling_ij * oL;
        if (vertex_j == parent)
        {
          data.oYaba_augmented[parent].noalias() +=
            crosscoupling_ij_oL + crosscoupling_ij_oL.transpose();
        }
        else
        {
          if (vertex_j < parent)
          {
            joint_cross_coupling.get({vertex_j, parent}).noalias() += crosscoupling_ij_oL;
          }
          else
          {
            joint_cross_coupling.get({parent, vertex_j}).noalias() +=
              crosscoupling_ij_oL.transpose();
          }
        }

        for (size_t k = j + 1; k < joint_neighbours.size(); ++k)
        {
          const JointIndex vertex_k = joint_neighbours[k];

          const Matrix6 & edge_ik =
            (i > vertex_k) ? joint_cross_coupling.get(JointPair(vertex_k, i))
                           : joint_cross_coupling.get(JointPair(i, vertex_k)).transpose();

          crosscoupling_ix_Jcols.noalias() = edge_ik * Jcols;

          assert(vertex_j != vertex_k && "Must never happen!");
          if (vertex_j < vertex_k)
          {
            joint_cross_coupling.get({vertex_j, vertex_k}).noalias() -=
              crosscoupling_ij_Jcols_Dinv
              * crosscoupling_ix_Jcols.transpose(); // Warning: UDinv() is actually edge_ik * J_col,
                                                    // U() is edge_ij * J_col * Dinv
          }
          else // if (vertex_k < vertex_j)
          {
            joint_cross_coupling.get({vertex_k, vertex_j}).transpose().noalias() -=
              crosscoupling_ij_Jcols_Dinv
              * crosscoupling_ix_Jcols.transpose(); // Warning: UDinv() is actually edge_ik *
                                                    // J_col, U() is edge_ij * J_col * Dinv
          }
        }
      }
    }
  };

  // A reduced backward sweep that only propagates the affine terms
  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  struct LCABAReducedBackwardStep
  : public fusion::JointUnaryVisitorBase<
      LCABAReducedBackwardStep<Scalar, Options, JointCollectionTpl>>
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
      typedef std::pair<JointIndex, JointIndex> JointPair;

      const auto & neighbours = data.neighbour_links;
      auto & joint_cross_coupling = data.joint_cross_coupling;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      typedef
        typename SizeDepType<JointModel::NV>::template ColsReturn<pinocchio::Data::Matrix6x>::Type
          ColBlock;
      const auto Jcols = jmodel.jointCols(data.J);

      Force & fi = data.of[i];

      jmodel.jointVelocitySelector(data.u).noalias() = -Jcols.transpose() * fi.toVector();

      jmodel.jointVelocitySelector(data.g).noalias() =
        jdata.Dinv()
        * jmodel.jointVelocitySelector(
          data.u); // Abuse of notation to reuse existing unused data variable

      const Vector6 a_tmp = Jcols * jmodel.jointVelocitySelector(data.g);

      for (JointIndex vertex_j : neighbours[i])
      {
        const Matrix6 & edge_ij = (i > vertex_j)
                                    ? joint_cross_coupling.get(JointPair(vertex_j, i))
                                    : joint_cross_coupling.get(JointPair(i, vertex_j)).transpose();
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
        const Matrix6 & edge_ij =
          (i > vertex_j) ? data.joint_cross_coupling.get(JointPair(vertex_j, i)).transpose()
                         : data.joint_cross_coupling.get(JointPair(i, vertex_j));
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
      typedef typename Data::Matrix6 Matrix6;

      typedef typename SizeDepType<JointModel::NV>::template ColsReturn<Matrix6x>::Type ColBlock;
      ColBlock J_cols = jmodel.jointCols(data.J);

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];
      const auto & neighbours = data.neighbour_links[i];

      data.oa_gf[i] = data.oa_gf[parent]; // does take into account the gravity field

      Force & fi = data.of[i];
      for (JointIndex vertex_j : neighbours)
      {

        const Matrix6 & edge_ij =
          (i > vertex_j) ? data.joint_cross_coupling.get(JointPair(vertex_j, i)).transpose()
                         : data.joint_cross_coupling.get(JointPair(i, vertex_j));
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

  namespace internal
  {
    template<typename ConstraintModel>
    struct LCABAConstraintCalcStep
    {
      typedef typename ConstraintModel::ConstraintData ConstraintData;
      typedef typename ConstraintModel::Scalar Scalar;

      template<typename Model, typename Data>
      static void run(
        const Model & model,
        Data & data,
        const ConstraintModel & cmodel,
        ConstraintData & cdata,
        const Scalar mu);
    };

    template<typename Scalar, int Options>
    struct LCABAConstraintCalcStep<RigidConstraintModelTpl<Scalar, Options>>
    {
      typedef RigidConstraintModelTpl<Scalar, Options> ConstraintModel;
      typedef typename ConstraintModel::ConstraintData ConstraintData;

      template<typename Model, typename Data>
      static void run(
        const Model & model,
        Data & data,
        const ConstraintModel & cmodel,
        ConstraintData & cdata,
        const Scalar mu)
      {
        typedef typename Data::SE3 SE3;
        typedef typename Model::JointIndex JointIndex;
        typedef std::pair<JointIndex, JointIndex> JointPair;
        typedef typename Data::Matrix6 Matrix6;
        typedef typename ConstraintModel::Matrix36 Matrix36;

        cdata.contact_force.setZero();

        cmodel.calc(model, data, cdata);

        SE3 & oMc1 = cdata.oMc1;
        SE3 & oMc2 = cdata.oMc2;
        SE3 & c1Mc2 = cdata.c1Mc2;
        typename ConstraintData::Motion & vc1 = cdata.contact1_velocity;
        typename ConstraintData::Motion & vc2 = cdata.contact2_velocity;
        const JointIndex joint1_id = cmodel.joint1_id;
        const JointIndex joint2_id = cmodel.joint2_id;

        const auto & corrector = cmodel.corrector;
        auto & contact_velocity_error = cdata.contact_velocity_error;

        if (joint1_id > 0)
          vc1 = oMc1.actInv(data.ov[joint1_id]);
        else
          vc1.setZero();
        if (joint2_id > 0)
          vc2 = oMc2.actInv(data.ov[joint2_id]);
        else
          vc2.setZero();
        const Motion vc2_in_frame1 = c1Mc2.act(vc2);

        if (cmodel.type == CONTACT_6D)
        {
          cdata.contact_placement_error = -log6(c1Mc2);
          contact_velocity_error = vc1 - vc2_in_frame1;
          const Matrix6 A1 = oMc1.toActionMatrixInverse();
          const Matrix6 A1tA1 = A1.transpose() * A1;
          data.oYaba_augmented[joint1_id].noalias() += mu * A1tA1;

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
              -(corrector.Kd.asDiagonal() * contact_velocity_error.toVector())
              - (corrector.Kp.asDiagonal() * cdata.contact_placement_error.toVector());
          }

          cdata.contact_acceleration_desired -= oMc1.actInv(data.oa[joint1_id]);
          cdata.contact_acceleration_desired -= cdata.contact_velocity_error.cross(vc2_in_frame1);

          if (joint2_id > 0)
          {
            cdata.contact_acceleration_desired += oMc1.actInv(data.oa[joint2_id]);

            const Matrix6 A2 =
              -A1; // only for 6D case. also used below for computing A2tA2 and A1tA2
            data.oYaba_augmented[joint2_id].noalias() += mu * A1tA1;
            data.of[joint2_id].toVector().noalias() +=
              A2.transpose()
              * (/*cdata.contact_force.toVector()*/ -mu * cdata.contact_acceleration_desired.toVector());

            const JointPair jp = joint1_id < joint2_id ? JointPair{joint1_id, joint2_id}
                                                       : JointPair{joint2_id, joint1_id};
            assert(data.joint_cross_coupling.exist(jp) && "Must never happen");
            data.joint_cross_coupling.get(jp) -= mu * A1tA1;
          }
          else
          {
            cdata.contact_acceleration_desired.toVector().noalias() -=
              A1 * model.gravity.toVector();
          }

          data.of[joint1_id].toVector().noalias() +=
            A1.transpose()
            * (/*cdata.contact_force.toVector()*/ -mu * cdata.contact_acceleration_desired.toVector());
        }
        else if (cmodel.type == CONTACT_3D)
        {
          const Matrix36 & A1 = oMc1.toActionMatrixInverse().template topRows<3>();
          data.oYaba_augmented[joint1_id].noalias() += mu * A1.transpose() * A1;

          if (check_expression_if_real<Scalar, false>(
                isZero(corrector.Kp, static_cast<Scalar>(0.))
                && isZero(corrector.Kd, static_cast<Scalar>(0.))))
          {
            cdata.contact_acceleration_desired.setZero();
          }
          else
          {
            cdata.contact_acceleration_desired.linear().noalias() =
              -(corrector.Kd.asDiagonal() * contact_velocity_error.linear())
              - (corrector.Kp.asDiagonal() * cdata.contact_placement_error.linear());
            cdata.contact_acceleration_desired.angular().setZero();
          }

          cdata.contact_acceleration_desired.linear().noalias() -=
            vc1.angular().cross(vc1.linear());
          if (joint2_id > 0)
          {
            const Matrix36 A2 =
              -c1Mc2.rotation()
              * (oMc2.toActionMatrixInverse().template topRows<3>()); // TODO:remove memalloc

            cdata.contact_acceleration_desired.linear().noalias() +=
              c1Mc2.rotation() * vc2.angular().cross(vc2.linear());
            data.oYaba_augmented[joint2_id].noalias() += mu * A2.transpose() * A2;
            data.of[joint2_id].toVector().noalias() +=
              A2.transpose()
              * (/*cdata.contact_force.toVector()*/ -mu * cdata.contact_acceleration_desired.linear());

            if (joint1_id < joint2_id)
            {
              data.joint_cross_coupling.get({joint1_id, joint2_id}).noalias() +=
                mu * A1.transpose() * A2;
            }
            else
            {
              data.joint_cross_coupling.get({joint2_id, joint1_id}).noalias() +=
                mu * A2.transpose() * A1;
            }
          }
          else
          {
            cdata.contact_acceleration_desired.linear().noalias() -= A1 * model.gravity.toVector();
          }

          data.of[joint1_id].toVector().noalias() +=
            A1.transpose()
            * (/*cdata.contact_force.toVector()*/ -mu * cdata.contact_acceleration_desired.linear());
        }
        else
        {
          assert(false && "Must never happen");
        }
      }
    };
  } // namespace internal

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentVectorType1,
    typename TangentVectorType2,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  inline const typename DataTpl<Scalar, Options, JointCollectionTpl>::TangentVectorType & lcaba(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType1> & v,
    const Eigen::MatrixBase<TangentVectorType2> & tau,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    ProximalSettingsTpl<Scalar> & settings)
  {

    assert(model.check(data) && "data is not consistent with model.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      q.size(), model.nq, "The joint configuration vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      v.size(), model.nv, "The joint velocity vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      tau.size(), model.nv, "The joint torque vector is not of right size");

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;
    typedef typename ConstraintModel::Matrix36 Matrix36;

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

    for (std::size_t i = 0; i < constraint_models.size(); ++i)
    {
      ConstraintData & cdata = constraint_datas[i];
      const ConstraintModel & cmodel = constraint_models[i];

      typedef internal::LCABAConstraintCalcStep<ConstraintModel> CalcStep;
      CalcStep::run(model, data, cmodel, cdata, mu);
    }

    typedef LCABABackwardStep<Scalar, Options, JointCollectionTpl> Pass2;

    const auto & elimination_order = data.elimination_order;

    for (JointIndex i : elimination_order)
    {
      Pass2::run(model.joints[i], data.joints[i], typename Pass2::ArgsType(model, data));
    }

    assert(settings.mu > 0 && "constrainedABA requires mu > 0.");

    typedef LCABAForwardStep2<Scalar, Options, JointCollectionTpl> Pass3;

    for (int it = int(elimination_order.size()) - 1; it >= 0; --it)
    {
      const JointIndex i = elimination_order[size_t(it)];
      if (data.constraints_supported_dim[i] > 0)
        Pass3::run(model.joints[i], data.joints[i], typename Pass3::ArgsType(model, data));
    }

    typedef LCABAReducedBackwardStep<Scalar, Options, JointCollectionTpl> ReducedPass2;
    typedef LCABAReducedForwardStep<Scalar, Options, JointCollectionTpl> ReducedPass3;
    data.g.setZero();
    int iter = 1;
    for (; iter < settings.max_iter; iter++)
    {
      settings.absolute_residual = Scalar(0);
      for (JointIndex j = 1; j < (JointIndex)model.njoints; ++j)
      {
        if (data.constraints_supported_dim[j] > 0)
          data.of[j].setZero();
      }
      for (std::size_t j = 0; j < constraint_models.size(); ++j)
      {
        const ConstraintModel & cmodel = constraint_models[j];
        ConstraintData & cdata = constraint_datas[j];
        const JointIndex joint1_id = cmodel.joint1_id;
        const JointIndex joint2_id = cmodel.joint2_id;
        typename ConstraintData::Motion & contact_acc_err = cdata.contact_acceleration_error;

        if (cmodel.type == CONTACT_6D)
        {
          if (joint2_id > 0)
            contact_acc_err = cdata.oMc1.actInv((data.oa[joint1_id] - data.oa[joint2_id]))
                              - cdata.contact_acceleration_desired;
          else
            contact_acc_err =
              cdata.oMc1.actInv((data.oa[joint1_id])) - cdata.contact_acceleration_desired;

          const auto mu_lambda = Force(mu * contact_acc_err.toVector());
          cdata.contact_force += mu_lambda;

          if (joint1_id > 0)
            data.of[joint1_id] += cdata.oMc1.act(mu_lambda);

          if (joint2_id > 0)
            data.of[joint2_id] -= cdata.oMc1.act(mu_lambda);
        }
        else if (cmodel.type == CONTACT_3D)
        {
          contact_acc_err.linear() = -cdata.contact_acceleration_desired.linear();
          if (joint1_id > 0)
            contact_acc_err.linear() += cdata.oMc1.actInv(data.oa[joint1_id]).linear();
          if (joint2_id > 0)
            contact_acc_err.linear() -=
              cdata.c1Mc2.rotation() * cdata.oMc2.actInv(data.oa[joint2_id]).linear();

          const auto mu_lambda = Force(mu * contact_acc_err.toVector());
          cdata.contact_force.linear() += mu_lambda.linear();

          if (joint1_id > 0)
            data.of[joint1_id] += cdata.oMc1.act(mu_lambda);

          if (joint2_id > 0)
          {
            const Matrix36 A2 =
              -cdata.c1Mc2.rotation() * (cdata.oMc2.toActionMatrixInverse().template topRows<3>());

            data.of[joint2_id].toVector().noalias() += A2.transpose() * mu_lambda.linear();
          }
        }

        const Scalar constraint_residual_norm =
          contact_acc_err.toVector().template lpNorm<Eigen::Infinity>();
        if (settings.absolute_residual < constraint_residual_norm)
          settings.absolute_residual = constraint_residual_norm;
      }

      if (settings.absolute_residual < settings.absolute_accuracy)
        break;

      // reduced backward sweep
      for (JointIndex j : elimination_order)
      {
        if (data.constraints_supported_dim[j] > 0)
          ReducedPass2::run(
            model.joints[j], data.joints[j], typename ReducedPass2::ArgsType(model, data));
      }

      // reduced forward sweep
      data.oa_gf[0].setZero();
      for (int it = int(elimination_order.size()) - 1; it >= 0; it--)
      {
        const JointIndex j = elimination_order[size_t(it)];
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
    for (int it = int(elimination_order.size()) - 1; it >= 0; --it)
    {
      const JointIndex j = elimination_order[size_t(it)];
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
