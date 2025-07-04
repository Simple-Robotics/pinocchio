//
// Copyright (c) 2019-2024 INRIA
//

#ifndef __pinocchio_algorithm_contact_cholesky_hxx__
#define __pinocchio_algorithm_contact_cholesky_hxx__

#include "pinocchio/algorithm/check.hpp"
#include "pinocchio/multibody/data.hpp"

#include <algorithm>

namespace pinocchio
{

  // TODO Remove when API is stabilized
  PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
  PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS

  template<typename Scalar, int Options>
  template<
    typename S1,
    int O1,
    template<typename, int> class JointCollectionTpl,
    template<typename T> class Holder,
    class ConstraintModel,
    class ConstraintAllocator>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::resize(
    const ModelTpl<S1, O1, JointCollectionTpl> & model,
    const std::vector<Holder<const ConstraintModel>, ConstraintAllocator> & contact_models)
  {
    typedef ModelTpl<S1, O1, JointCollectionTpl> Model;
    typedef typename Model::JointModel JointModel;

    static_assert(
      std::is_base_of<ConstraintModelBase<ConstraintModel>, ConstraintModel>::value,
      "ConstraintModel is not a ConstraintModelBase");

    assert(model.check(MimicChecker()) && "Function does not support mimic joints");

    nv = model.nv;

    Eigen::DenseIndex num_total_constraints = 0;
    for (std::size_t i = 0; i < contact_models.size(); i++)
    {
      const ConstraintModel & cmodel = contact_models[(std::size_t)i];
      num_total_constraints += cmodel.activeSize();
    }

    const Eigen::DenseIndex total_dim = nv + num_total_constraints;

    // Compute first parents_fromRow for all the joints.
    // This code is very similar to the code of Data::computeParents_fromRow,
    // but shifted with a value corresponding to the number of constraints.
    parents_fromRow.resize(total_dim);
    parents_fromRow.fill(-1);

    nv_subtree_fromRow.resize(total_dim);
    //      nv_subtree_fromRow.fill(0);

    // Fill nv_subtree_fromRow for model
    for (JointIndex joint_id = 1; joint_id < (JointIndex)(model.njoints); joint_id++)
    {
      const JointIndex parent_id = model.parents[joint_id];

      const JointModel joint = model.joints[joint_id];
      const JointModel parent_joint = model.joints[parent_id];
      const int nvj = joint.nv();
      const int idx_vj = joint.idx_v();

      if (parent_id > 0)
        parents_fromRow[idx_vj + num_total_constraints] =
          parent_joint.idx_v() + parent_joint.nv() - 1 + num_total_constraints;
      else
        parents_fromRow[idx_vj + num_total_constraints] = -1;

      const JointIndex last_child =
        model.subtrees[joint_id].size() > 0 ? model.subtrees[joint_id].back() : JointIndex(0);
      nv_subtree_fromRow[idx_vj + num_total_constraints] =
        model.joints[last_child].idx_v() + model.joints[last_child].nv() - idx_vj;

      for (int row = 1; row < nvj; ++row)
      {
        parents_fromRow[idx_vj + num_total_constraints + row] =
          idx_vj + row - 1 + num_total_constraints;
        nv_subtree_fromRow[idx_vj + num_total_constraints + row] =
          nv_subtree_fromRow[idx_vj + num_total_constraints] - row;
      }
    }

    Eigen::DenseIndex row_id = 0;
    for (std::size_t i = 0; i < contact_models.size(); i++)
    {
      const ConstraintModel & cmodel = contact_models[(std::size_t)i];
      for (Eigen::DenseIndex k = 0; k < cmodel.activeSize(); ++k, row_id++)
      {
        const auto & row_active_indexes = cmodel.getRowActiveIndexes(k);
        nv_subtree_fromRow[row_id] =
          num_total_constraints - row_id + 1
          + (row_active_indexes.size() > 0 ? row_active_indexes.back() : 0);
      }
    }
    assert(row_id == num_total_constraints);

    // Fill the sparsity pattern for each Row of the Cholesky decomposition (matrix U)
    /*
          static const Slice default_slice_value(1,1);
          static const SliceVector default_slice_vector(1,default_slice_value);

          rowise_sparsity_pattern.clear();
          rowise_sparsity_pattern.resize((size_t)num_total_constraints,default_slice_vector);
          row_id = 0; size_t ee_id = 0;
          for(typename RigidConstraintModelVector::const_iterator it = contact_models.begin();
              it != contact_models.end();
              ++it, ++ee_id)
          {
            const RigidConstraintModel & cmodel = *it;
            const BooleanVector & joint1_indexes_ee = cmodel.colwise_joint1_sparsity;
            const Eigen::DenseIndex contact_dim = cmodel.size();

            for(Eigen::DenseIndex k = 0; k < contact_dim; ++k)
            {
              SliceVector & slice_vector = rowise_sparsity_pattern[(size_t)row_id];
              slice_vector.clear();
              slice_vector.push_back(Slice(row_id,num_total_constraints-row_id));

              bool previous_index_was_true = true;
              for(Eigen::DenseIndex joint1_indexes_ee_id = num_total_constraints;
                  joint1_indexes_ee_id < total_dim;
                  ++joint1_indexes_ee_id)
              {
                if(joint1_indexes_ee[joint1_indexes_ee_id])
                {
                  if(previous_index_was_true) // no discontinuity
                    slice_vector.back().size++;
                  else // discontinuity; need to create a new slice
                  {
                    const Slice new_slice(joint1_indexes_ee_id,1);
                    slice_vector.push_back(new_slice);
                  }
                }

                previous_index_was_true = joint1_indexes_ee[joint1_indexes_ee_id];
              }

              row_id++;
            }
          }
     */

    // Allocate Eigen memory if needed
    compliance_storage.resize(num_total_constraints);
    compliance.setZero();
    damping_storage.resize(num_total_constraints);
    damping.setZero();

    D_storage.resize(total_dim);
    Dinv_storage.resize(total_dim);
    U_storage.resize(total_dim, total_dim);
    U.setIdentity();
    DUt_storage.resize(total_dim);
  }

  template<typename Scalar, int Options>
  template<
    typename S1,
    int O1,
    template<typename, int> class JointCollectionTpl,
    template<typename T> class Holder,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename VectorLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::compute(
    const ModelTpl<S1, O1, JointCollectionTpl> & model,
    DataTpl<S1, O1, JointCollectionTpl> & data,
    const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & contact_models,
    std::vector<Holder<ConstraintData>, ConstraintDataAllocator> & contact_datas,
    const Eigen::MatrixBase<VectorLike> & mus)
  {
    static_assert(
      std::is_base_of<ConstraintModelBase<ConstraintModel>, ConstraintModel>::value,
      "ConstraintModel is not a ConstraintModelBase");
    static_assert(
      std::is_base_of<ConstraintDataBase<ConstraintData>, ConstraintData>::value,
      "ConstraintData is not a ConstraintDataBase");

    assert(model.check(data) && "data is not consistent with model.");
    assert(model.check(MimicChecker()) && "Function does not support mimic joints");

    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      contact_models.size() == contact_datas.size(),
      "The number of constraints between contact_models and contact_datas vectors is different.");
    PINOCCHIO_ONLY_USED_FOR_DEBUG(model);

    const Eigen::DenseIndex total_dim = size();
    const Eigen::DenseIndex total_constraints_dim = total_dim - nv;

    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
    const typename Data::MatrixXs & M = data.M;

    const size_t num_ee = contact_models.size();

    // Update frame placements if needed
    for (size_t ee_id = 0; ee_id < num_ee; ++ee_id)
    {
      const ConstraintModel & cmodel = contact_models[ee_id].get();
      ConstraintData & cdata = contact_datas[ee_id].get();
      // TODO: should we call resize on cmodel ?
      cmodel.calc(model, data, cdata);
    }

    D.tail(model.nv) = M.diagonal();
    U.bottomRightCorner(model.nv, model.nv).template triangularView<Eigen::StrictlyUpper>() =
      M.template triangularView<Eigen::StrictlyUpper>();

    // Constraint filling
    Eigen::DenseIndex current_row = 0;
    U.topRightCorner(total_constraints_dim, model.nv).setZero();
    for (size_t ee_id = 0; ee_id < num_ee; ++ee_id)
    {
      const ConstraintModel & cmodel = contact_models[ee_id].get();
      ConstraintData & cdata = contact_datas[ee_id].get();

      const Eigen::DenseIndex constraint_dim = cmodel.activeSize();
      cmodel.jacobian(
        model, data, cdata, U.block(current_row, total_constraints_dim, constraint_dim, model.nv));
      current_row += constraint_dim;
    }

    // Cholesky
    for (Eigen::DenseIndex j = nv - 1; j >= 0; --j)
    {
      // Classic Cholesky decomposition related to the mass matrix
      const Eigen::DenseIndex jj = total_constraints_dim + j; // shifted index
      const Eigen::DenseIndex NVT = nv_subtree_fromRow[jj] - 1;
      auto DUt_partial = DUt.head(NVT);

      if (NVT)
        DUt_partial.noalias() =
          U.row(jj).segment(jj + 1, NVT).transpose().cwiseProduct(D.segment(jj + 1, NVT));

      D[jj] -= U.row(jj).segment(jj + 1, NVT).dot(DUt_partial);
      assert(
        check_expression_if_real<Scalar>(D[jj] != Scalar(0))
        && "The diagonal element is equal to zero.");
      Dinv[jj] = Scalar(1) / D[jj];

      for (Eigen::DenseIndex _ii = parents_fromRow[jj]; _ii >= total_constraints_dim;
           _ii = parents_fromRow[_ii])
      {
        U(_ii, jj) -= U.row(_ii).segment(jj + 1, NVT).dot(DUt_partial);
        U(_ii, jj) *= Dinv[jj];
      }

      // Constraint part
      Eigen::DenseIndex current_row = total_constraints_dim - 1;
      for (size_t ee_id = 0; ee_id < num_ee; ++ee_id)
      {
        const ConstraintModel & cmodel = contact_models[num_ee - 1 - ee_id];
        const Eigen::DenseIndex constraint_dim = cmodel.activeSize();

        for (Eigen::DenseIndex constraint_row_id = constraint_dim - 1; constraint_row_id >= 0;
             --constraint_row_id, --current_row)
        {
          const auto & colwise_sparsity = cmodel.getRowActiveSparsityPattern(constraint_row_id);
          if (colwise_sparsity[j])
          {
            U(current_row, jj) -= U.row(current_row).segment(jj + 1, NVT).dot(DUt_partial);
            U(current_row, jj) *= Dinv[jj];
          }
        }
      }
    }

    // Setting physical compliance
    int cindex = 0;
    for (std::size_t ee_id = 0; ee_id < num_ee; ee_id++)
    {
      const ConstraintModel & cmodel = contact_models[ee_id];
      // TODO use active compliance
      const int cdim = cmodel.activeSize();
      compliance.segment(cindex, cdim) = cmodel.getActiveCompliance();
      cindex += cdim;
    }

    // Setting numerical damping
    updateDamping(mus);
  }

  template<typename Scalar, int Options>
  template<typename VectorLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::updateCompliance(
    const Eigen::MatrixBase<VectorLike> & vec)
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
    compliance = vec;

    // The diagonal term of the KKT should be updated with the new compliance
    updateDamping(getDamping());
  }

  template<typename Scalar, int Options>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::updateCompliance(const Scalar & compliance)
  {
    const Eigen::DenseIndex total_dim = size();
    const Eigen::DenseIndex total_constraints_dim = total_dim - nv;
    updateCompliance(Vector::Constant(total_constraints_dim, compliance));
  }

  template<typename Scalar, int Options>
  template<typename VectorLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::updateDamping(
    const Eigen::MatrixBase<VectorLike> & vec)
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
    damping = vec;
    const Eigen::DenseIndex total_dim = size();
    const Eigen::DenseIndex total_constraints_dim = total_dim - nv;

    // Upper left triangular part of U
    for (Eigen::DenseIndex j = total_constraints_dim - 1; j >= 0; --j)
    {
      const Eigen::DenseIndex slice_dim = total_dim - j - 1;
      auto DUt_partial = DUt.head(slice_dim);
      DUt_partial.noalias() =
        U.row(j).segment(j + 1, slice_dim).transpose().cwiseProduct(D.segment(j + 1, slice_dim));

      D[j] = -damping[j] - compliance[j] - U.row(j).segment(j + 1, slice_dim).dot(DUt_partial);
      assert(
        check_expression_if_real<Scalar>(D[j] != Scalar(0))
        && "The diagonal element is equal to zero.");
      Dinv[j] = Scalar(1) / D[j];

      for (Eigen::DenseIndex _i = j - 1; _i >= 0; _i--)
      {
        U(_i, j) = -U.row(_i).segment(j + 1, slice_dim).dot(DUt_partial) * Dinv[j];
      }
    }
  }

  template<typename Scalar, int Options>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::updateDamping(const Scalar & mu)
  {
    //      PINOCCHIO_CHECK_INPUT_ARGUMENT(check_expression_if_real<Scalar>(mu >= 0), "mu should be
    //      positive.");

    const Eigen::DenseIndex total_dim = size();
    const Eigen::DenseIndex total_constraints_dim = total_dim - nv;
    updateDamping(Vector::Constant(total_constraints_dim, mu));
  }

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::solveInPlace(
    const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    MatrixLike & mat_ = PINOCCHIO_EIGEN_CONST_CAST(MatrixLike, mat);

    Uiv(mat_);
    mat_.array().colwise() *= Dinv.array();
    Utiv(mat_);
  }

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  typename ContactCholeskyDecompositionTpl<Scalar, Options>::Matrix
  ContactCholeskyDecompositionTpl<Scalar, Options>::solve(
    const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    Matrix res(mat);
    solveInPlace(res);
    return res;
  }

  template<typename Scalar, int Options>
  template<typename S1, int O1, template<typename, int> class JointCollectionTpl>
  ContactCholeskyDecompositionTpl<Scalar, Options>
  ContactCholeskyDecompositionTpl<Scalar, Options>::getMassMatrixChoeslkyDecomposition(
    const ModelTpl<S1, O1, JointCollectionTpl> & model) const
  {
    typedef ContactCholeskyDecompositionTpl<Scalar, Options> ReturnType;
    ReturnType res(model);

    res.D = D.tail(nv);
    res.Dinv = Dinv.tail(nv);
    res.U = U.bottomRightCorner(nv, nv);

    return res;
  }

  namespace details
  {
    template<typename MatrixLike, int ColsAtCompileTime>
    struct UvAlgo
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<MatrixLike> & mat)
      {
        MatrixLike & mat_ = PINOCCHIO_EIGEN_CONST_CAST(MatrixLike, mat);

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          mat.rows() == chol.size(), "The input matrix is of wrong size");

        for (Eigen::DenseIndex col_id = 0; col_id < mat_.cols(); ++col_id)
          UvAlgo<typename MatrixLike::ColXpr>::run(chol, mat_.col(col_id));
      }
    };

    template<typename VectorLike>
    struct UvAlgo<VectorLike, 1>
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<VectorLike> & vec)
      {
        EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
        VectorLike & vec_ = PINOCCHIO_EIGEN_CONST_CAST(VectorLike, vec);

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          vec.size() == chol.size(), "The input vector is of wrong size");
        const Eigen::DenseIndex num_total_constraints = chol.size() - chol.nv;

        // TODO: exploit the Sparsity pattern of the first rows of U
        for (Eigen::DenseIndex k = 0; k < num_total_constraints; ++k)
        {
          const Eigen::DenseIndex slice_dim = chol.size() - k - 1;
          vec_[k] += chol.U.row(k).tail(slice_dim).dot(vec_.tail(slice_dim));
        }

        for (Eigen::DenseIndex k = num_total_constraints; k <= chol.size() - 2; ++k)
          vec_[k] += chol.U.row(k)
                       .segment(k + 1, chol.nv_subtree_fromRow[k] - 1)
                       .dot(vec_.segment(k + 1, chol.nv_subtree_fromRow[k] - 1));
      }
    };
  } // namespace details

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::Uv(
    const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    details::UvAlgo<MatrixLike>::run(*this, mat.const_cast_derived());
  }

  namespace details
  {
    template<typename MatrixLike, int ColsAtCompileTime>
    struct UtvAlgo
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<MatrixLike> & mat)
      {
        MatrixLike & mat_ = mat.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          mat.rows() == chol.size(), "The input matrix is of wrong size");

        for (Eigen::DenseIndex col_id = 0; col_id < mat_.cols(); ++col_id)
          UtvAlgo<typename MatrixLike::ColXpr>::run(chol, mat_.col(col_id));
      }
    };

    template<typename VectorLike>
    struct UtvAlgo<VectorLike, 1>
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<VectorLike> & vec)
      {
        EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
        VectorLike & vec_ = vec.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          vec.size() == chol.size(), "The input vector is of wrong size");
        const Eigen::DenseIndex num_total_constraints = chol.constraintDim();

        for (Eigen::DenseIndex k = chol.size() - 2; k >= num_total_constraints; --k)
          vec_.segment(k + 1, chol.nv_subtree_fromRow[k] - 1) +=
            chol.U.row(k).segment(k + 1, chol.nv_subtree_fromRow[k] - 1).transpose() * vec_[k];

        // TODO: exploit the Sparsity pattern of the first rows of U
        for (Eigen::DenseIndex k = num_total_constraints - 1; k >= 0; --k)
        {
          const Eigen::DenseIndex slice_dim = chol.size() - k - 1;
          vec_.tail(slice_dim) += chol.U.row(k).tail(slice_dim).transpose() * vec_[k];
        }
      }
    };
  } // namespace details

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::Utv(
    const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    details::UtvAlgo<MatrixLike>::run(*this, mat.const_cast_derived());
  }

  namespace details
  {
    template<typename MatrixLike, int ColsAtCompileTime>
    struct UivAlgo
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<MatrixLike> & mat)
      {
        MatrixLike & mat_ = mat.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          mat.rows() == chol.size(), "The input matrix is of wrong size");

        for (Eigen::DenseIndex col_id = 0; col_id < mat_.cols(); ++col_id)
          UivAlgo<typename MatrixLike::ColXpr>::run(chol, mat_.col(col_id));
      }
    };

    template<typename VectorLike>
    struct UivAlgo<VectorLike, 1>
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<VectorLike> & vec)
      {
        EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
        VectorLike & vec_ = vec.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          vec.size() == chol.size(), "The input vector is of wrong size");

        const Eigen::DenseIndex num_total_constraints = chol.size() - chol.nv;
        for (Eigen::DenseIndex k = chol.size() - 2; k >= num_total_constraints; --k)
          vec_[k] -= chol.U.row(k)
                       .segment(k + 1, chol.nv_subtree_fromRow[k] - 1)
                       .dot(vec_.segment(k + 1, chol.nv_subtree_fromRow[k] - 1));

        // TODO: exploit the Sparsity pattern of the first rows of U
        for (Eigen::DenseIndex k = num_total_constraints - 1; k >= 0; --k)
        {
          const Eigen::DenseIndex slice_dim = chol.size() - k - 1;
          vec_[k] -= chol.U.row(k).tail(slice_dim).dot(vec_.tail(slice_dim));
        }
      }
    };
  } // namespace details

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::Uiv(
    const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    details::UivAlgo<MatrixLike>::run(*this, mat.const_cast_derived());
  }

  namespace details
  {
    template<typename MatrixLike, int ColsAtCompileTime>
    struct UtivAlgo
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<MatrixLike> & mat)
      {
        MatrixLike & mat_ = mat.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          mat.rows() == chol.size(), "The input matrix is of wrong size");

        for (Eigen::DenseIndex col_id = 0; col_id < mat_.cols(); ++col_id)
          UtivAlgo<typename MatrixLike::ColXpr>::run(chol, mat_.col(col_id));
      }
    };

    template<typename VectorLike>
    struct UtivAlgo<VectorLike, 1>
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<VectorLike> & vec)
      {
        EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
        VectorLike & vec_ = vec.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          vec.size() == chol.size(), "The input vector is of wrong size");
        const Eigen::DenseIndex num_total_constraints = chol.constraintDim();

        // TODO: exploit the Sparsity pattern of the first rows of U
        for (Eigen::DenseIndex k = 0; k < num_total_constraints; ++k)
        {
          const Eigen::DenseIndex slice_dim = chol.size() - k - 1;
          vec_.tail(slice_dim) -= chol.U.row(k).tail(slice_dim).transpose() * vec_[k];
        }

        for (Eigen::DenseIndex k = num_total_constraints; k <= chol.size() - 2; ++k)
          vec_.segment(k + 1, chol.nv_subtree_fromRow[k] - 1) -=
            chol.U.row(k).segment(k + 1, chol.nv_subtree_fromRow[k] - 1).transpose() * vec_[k];
      }
    };
  } // namespace details

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::Utiv(
    const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    details::UtivAlgo<MatrixLike>::run(*this, PINOCCHIO_EIGEN_CONST_CAST(MatrixLike, mat));
  }

  template<typename Scalar, int Options>
  typename ContactCholeskyDecompositionTpl<Scalar, Options>::Matrix
  ContactCholeskyDecompositionTpl<Scalar, Options>::matrix() const
  {
    Matrix res(size(), size());
    matrix(res);
    return res;
  }

  template<typename Scalar, int Options>
  template<typename MatrixType>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::matrix(
    const Eigen::MatrixBase<MatrixType> & res) const
  {
    MatrixType & res_ = PINOCCHIO_EIGEN_CONST_CAST(MatrixType, res);
    res_.noalias() = U * D.asDiagonal() * U.transpose();
  }

  template<typename Scalar, int Options>
  typename ContactCholeskyDecompositionTpl<Scalar, Options>::Matrix
  ContactCholeskyDecompositionTpl<Scalar, Options>::inverse() const
  {
    Matrix res(size(), size());
    inverse(res);
    return res;
  }

  template<typename Scalar, int Options>
  template<typename S1, int O1>
  bool ContactCholeskyDecompositionTpl<Scalar, Options>::operator==(
    const ContactCholeskyDecompositionTpl<S1, O1> & other) const
  {
    bool is_same = true;

    if (nv != other.nv)
      return false;

    if (
      D.size() != other.D.size() || Dinv.size() != other.Dinv.size() || U.rows() != other.U.rows()
      || U.cols() != other.U.cols())
      return false;

    is_same &= (D == other.D);
    is_same &= (Dinv == other.Dinv);
    is_same &= (U == other.U);

    is_same &= (parents_fromRow == other.parents_fromRow);
    is_same &= (nv_subtree_fromRow == other.nv_subtree_fromRow);
    //        is_same &= (rowise_sparsity_pattern == other.rowise_sparsity_pattern);

    return is_same;
  }

  template<typename Scalar, int Options>
  template<typename S1, int O1>
  bool ContactCholeskyDecompositionTpl<Scalar, Options>::operator!=(
    const ContactCholeskyDecompositionTpl<S1, O1> & other) const
  {
    return !(*this == other);
  }

  namespace details
  {

    template<typename Scalar, int Options, typename VectorLike>
    PINOCCHIO_DONT_INLINE VectorLike & inverseAlgo(
      const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
      const Eigen::DenseIndex col,
      const Eigen::MatrixBase<VectorLike> & vec)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike);

      typedef ContactCholeskyDecompositionTpl<Scalar, Options> ContactCholeskyDecomposition;

      const Eigen::DenseIndex & chol_dim = chol.size();
      PINOCCHIO_CHECK_INPUT_ARGUMENT(col < chol_dim && col >= 0);
      PINOCCHIO_CHECK_INPUT_ARGUMENT(vec.size() == chol_dim);

      const typename ContactCholeskyDecomposition::EigenIndexVector & nvt = chol.nv_subtree_fromRow;
      VectorLike & vec_ = PINOCCHIO_EIGEN_CONST_CAST(VectorLike, vec);

      const Eigen::DenseIndex last_col =
        std::min(col - 1, chol_dim - 2); // You can start from nv-2 (no child in nv-1)
      vec_[col] = Scalar(1);
      vec_.tail(chol_dim - col - 1).setZero();

      // TODO: exploit the sparsity pattern of the first rows of U
      for (Eigen::DenseIndex k = last_col; k >= 0; --k)
      {
        const Eigen::DenseIndex nvt_max = std::min(col - k, nvt[k] - 1);
        const auto U_row = chol.U.row(k);
        vec_[k] = -U_row.segment(k + 1, nvt_max).dot(vec_.segment(k + 1, nvt_max));
        //          if(k >= chol_constraint_dim)
        //          {
        //            vec_[k] = -U_row.segment(k+1,nvt_max).dot(vec_.segment(k+1,nvt_max));
        //          }
        //          else
        //          {
        //            typedef typename ContactCholeskyDecomposition::SliceVector SliceVector;
        //            typedef typename ContactCholeskyDecomposition::Slice Slice;
        //            const SliceVector & slice_vector = chol.rowise_sparsity_pattern[(size_t)k];
        //
        //            const Slice & slice_0 = slice_vector[0];
        //            assert(slice_0.first_index == k);
        //            Eigen::DenseIndex last_index1 = slice_0.first_index + slice_0.size;
        //            const Eigen::DenseIndex last_index2 = k + nvt_max;
        //            Eigen::DenseIndex slice_dim = std::min(last_index1,last_index2) - k;
        //            vec_[k] =
        //            -U_row.segment(slice_0.first_index+1,slice_dim-1).dot(vec_.segment(slice_0.first_index+1,slice_dim-1));
        //
        //            typename SliceVector::const_iterator slice_it = slice_vector.begin()++;
        //            for(;slice_it != slice_vector.end(); ++slice_it)
        //            {
        //              const Slice & slice = *slice_it;
        //              last_index1 = slice.first_index + slice.size;
        //              slice_dim = std::min(last_index1,last_index2+1) - slice.first_index;
        //              if(slice_dim <= 0) break;
        //
        //              vec_[k] -=
        //              U_row.segment(slice.first_index,slice_dim).dot(vec_.segment(slice.first_index,slice_dim));
        //            }
        //          }
      }

      vec_.head(col + 1).array() *= chol.Dinv.head(col + 1).array();

      for (Eigen::DenseIndex k = 0; k < col + 1; ++k) // You can stop one step before nv.
      {
        const Eigen::DenseIndex nvt_max = nvt[k] - 1;
        vec_.segment(k + 1, nvt_max) -= chol.U.row(k).segment(k + 1, nvt_max).transpose() * vec_[k];
      }

      return vec_;
    }
  } // namespace details

  template<typename Scalar, int Options>
  template<typename MatrixType>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::inverse(
    const Eigen::MatrixBase<MatrixType> & res) const
  {
    PINOCCHIO_CHECK_INPUT_ARGUMENT(res.rows() == size());
    PINOCCHIO_CHECK_INPUT_ARGUMENT(res.cols() == size());

    MatrixType & res_ = PINOCCHIO_EIGEN_CONST_CAST(MatrixType, res);

    for (Eigen::DenseIndex col_id = 0; col_id < size(); ++col_id)
      details::inverseAlgo(*this, col_id, res_.col(col_id));

    res_.template triangularView<Eigen::StrictlyLower>() =
      res_.transpose().template triangularView<Eigen::StrictlyLower>();
  }

  PINOCCHIO_COMPILER_DIAGNOSTIC_POP
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_contact_cholesky_hxx__
