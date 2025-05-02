//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_visitors_hxx__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_visitors_hxx__

namespace pinocchio
{

  template<typename DelassusOperator, typename ConfigVectorType>
  struct DelassusOperatorRigidBodySystemsComputeForwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodySystemsComputeForwardPass<DelassusOperator, ConfigVectorType>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;

    typedef boost::fusion::vector<const Model &, Data &, const ConfigVectorType &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      Data & data,
      const Eigen::MatrixBase<ConfigVectorType> & q)
    {
      typedef typename Model::JointIndex JointIndex;

      const JointIndex i = jmodel.id();
      jmodel.calc(jdata.derived(), q.derived());

      const JointIndex parent = model.parents[i];
      data.liMi[i] = model.jointPlacements[i] * jdata.M();
      auto & oMi = data.oMi[i];
      if (parent > 0)
        oMi = data.oMi[parent] * data.liMi[i];
      else
        oMi = data.liMi[i];

      // ABA in WORLD frame requires these quantities
      jmodel.jointCols(data.J) = oMi.act(jdata.S());
    }
  };

  template<typename DelassusOperator, bool damping_compliance_update_only = false>
  struct DelassusOperatorRigidBodySystemsComputeBackwardPass
  : public fusion::JointUnaryVisitorBase<DelassusOperatorRigidBodySystemsComputeBackwardPass<
      DelassusOperator,
      damping_compliance_update_only>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
    typedef typename Model::Scalar Scalar;

    typedef boost::fusion::vector<const Model &, Data &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      Data & data)
    {
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Matrix6 Matrix6;
      typedef typename JointModel::JointDataDerived JointData;
      typedef std::pair<JointIndex, JointIndex> JointPair;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      // ApplyOnTheRight
      if (!damping_compliance_update_only)
      {
        auto & Ia = data.Yaba[i];
        jmodel.calc_aba(
          jdata.derived(), jmodel.jointVelocitySelector(model.armature), Ia, parent > 0);
        if (parent > 0)
        {
          data.Yaba[parent] += impl::internal::SE3actOn<Scalar>::run(data.liMi[i], Ia);
        }
      }

      // SolveInPlace
      {
        JointData & _jdata_augmented = boost::get<JointData>(data.joints_augmented[i]);
        JointDataBase<JointData> & jdata_augmented =
          static_cast<JointDataBase<JointData> &>(_jdata_augmented);

        auto Jcols = jmodel.jointCols(data.J);
        auto & Ia_augmented = data.oYaba_augmented[i];

        jdata_augmented.U().noalias() = Ia_augmented * Jcols;
        jdata_augmented.StU().noalias() = Jcols.transpose() * jdata_augmented.U();

        // Account for the rotor inertia contribution
        jdata_augmented.StU().diagonal() += jmodel.jointVelocitySelector(model.armature);

        pinocchio::internal::PerformStYSInversion<Scalar>::run(
          jdata_augmented.StU(), jdata_augmented.Dinv());

        jdata_augmented.UDinv().noalias() = jdata_augmented.U() * jdata_augmented.Dinv();

        if (parent > 0)
        {
          Ia_augmented.noalias() -= jdata_augmented.UDinv() * jdata_augmented.U().transpose();
          data.oYaba_augmented[parent] += Ia_augmented;
        }

        // End of the classic ABA backward pass - beginning of cross-coupling handling
        const auto & neighbours = data.neighbour_links;
        auto & joint_cross_coupling = data.joint_cross_coupling;
        const auto & joint_neighbours = neighbours[i];

        if (joint_neighbours.size() == 0)
          return; // We can return from this point as this joint has no neighbours

        //        return;
        using Matrix6xNV = typename std::remove_reference<typename JointData::UDTypeRef>::type;
        typedef Eigen::Map<Matrix6xNV> MapMatrix6xNV;
        MapMatrix6xNV mat1_tmp = MapMatrix6xNV(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, 6, jmodel.nv()));
        MapMatrix6xNV mat2_tmp = MapMatrix6xNV(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, 6, jmodel.nv()));

        auto & JDinv = mat1_tmp;
        JDinv.noalias() = Jcols * jdata_augmented.Dinv();

        // oL == data.oL[i]
        Matrix6 oL = -JDinv * jdata_augmented.U().transpose();
        oL += Matrix6::Identity();

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
          crosscoupling_ij_Jcols_Dinv.noalias() = crosscoupling_ix_Jcols * jdata_augmented.Dinv();

          data.oYaba_augmented[vertex_j].noalias() -=
            crosscoupling_ij_Jcols_Dinv
            * crosscoupling_ix_Jcols.transpose(); // Warning: UDinv() is actually edge_ij * J, U()
                                                  // is actually edge_ij * J_cols * Dinv
                                                  //          data.of[vertex_j].toVector().noalias()
                                                  //          += crosscoupling_ij * a_tmp;

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
                * crosscoupling_ix_Jcols.transpose(); // Warning: UDinv() is actually edge_ik *
                                                      // J_col, U() is edge_ij * J_col * Dinv
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
    }
  };

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodySystemsTplApplyOnTheRightBackwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodySystemsTplApplyOnTheRightBackwardPass<DelassusOperator>>
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
      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      jmodel.jointVelocitySelector(custom_data.u) -= jdata.S().transpose() * custom_data.f[i];

      if (parent > 0)
      {
        auto & pa = custom_data.f[i];
        pa.toVector().noalias() += jdata.UDinv() * jmodel.jointVelocitySelector(custom_data.u);
        custom_data.f[parent] += data.liMi[i].act(pa);
      }
    }
  };

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass<DelassusOperator>>
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
      typedef typename Model::JointIndex JointIndex;

      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      //      typename JointData::TangentVector_t ddq_joint;
      auto ddq_joint = jmodel.jointVelocitySelector(custom_data.ddq);
      if (parent > 0)
      {
        custom_data.a[i] += data.liMi[i].actInv(custom_data.a[parent]);
        ddq_joint = jdata.Dinv() * jmodel.jointVelocitySelector(custom_data.u)
                    - jdata.UDinv().transpose() * custom_data.a[i].toVector();
        custom_data.a[i] += jdata.S() * ddq_joint;
      }
      else
      {
        ddq_joint.noalias() = jdata.Dinv() * jmodel.jointVelocitySelector(custom_data.u);
        custom_data.a[i] = jdata.S() * ddq_joint;
      }
    }

  }; // struct DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass

  template<typename DelassusOperator>
  struct AugmentedMassMatrixOperatorSolveInPlaceBackwardPass
  : public fusion::JointUnaryVisitorBase<
      AugmentedMassMatrixOperatorSolveInPlaceBackwardPass<DelassusOperator>>
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
  struct AugmentedMassMatrixOperatorSolveInPlaceForwardPass
  : public fusion::JointUnaryVisitorBase<
      AugmentedMassMatrixOperatorSolveInPlaceForwardPass<DelassusOperator>>
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
