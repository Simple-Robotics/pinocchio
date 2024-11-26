//
// Copyright (c) 2024 INRIA
//

#include <eigenpy/memory.hpp>
#include <eigenpy/eigen-from-python.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/bindings/python/fwd.hpp"

#include "pinocchio/algorithm/admm-solver.hpp"
#include "pinocchio/algorithm/contact-cholesky.hpp"
#include "pinocchio/algorithm/delassus-operator-dense.hpp"
#include "pinocchio/algorithm/delassus-operator-sparse.hpp"

#include "pinocchio/bindings/python/algorithm/contact-solver-base.hpp"
#include "pinocchio/bindings/python/utils/std-vector.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    typedef ADMMContactSolverTpl<context::Scalar> Solver;
    typedef Solver::SolverStats SolverStats;
    typedef context::Scalar Scalar;
    typedef context::VectorXs VectorXs;
    typedef const Eigen::Ref<const VectorXs> ConstRefVectorXs;
    typedef ContactCholeskyDecompositionTpl<context::Scalar, context::Options>
      ContactCholeskyDecomposition;

#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE

    template<typename DelassusDerived, typename ConstraintModel, typename ConstraintModelAllocator>
    static bool solve_wrapper(
      Solver & solver,
      DelassusDerived & delassus,
      const context::VectorXs & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const context::VectorXs & R,
      const boost::optional<ConstRefVectorXs> primal_solution = boost::none,
      const boost::optional<ConstRefVectorXs> dual_solution = boost::none,
      bool solve_ncp = true,
      bool compute_largest_eigen_values = true,
      ADMMUpdateRule admm_update_rule = ADMMUpdateRule::SPECTRAL,
      bool stat_record = false)
    {
      return solver.solve(
        delassus, g, constraint_models, R, primal_solution, dual_solution, solve_ncp,
        compute_largest_eigen_values, admm_update_rule, stat_record);
    }

    template<typename DelassusDerived, typename ConstraintModel, typename ConstraintModelAllocator>
    static bool solve_wrapper2(
      Solver & solver,
      DelassusDerived & delassus,
      const context::VectorXs & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      Eigen::Ref<context::VectorXs> x,
      bool solve_ncp = true)
    {
      return solver.solve(delassus, g, constraint_models, x, solve_ncp);
    }
#endif

#ifndef PINOCCHIO_PYTHON_SKIP_CASADI_UNSUPPORTED

    template<typename ConstraintModel, typename ConstraintModelAllocator>
    static context::VectorXs computeConeProjection_wrapper(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const context::VectorXs & forces)
    {
      context::VectorXs res(forces.size());

      typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
      typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

      WrappedConstraintModelVector wrapped_constraint_models(
        constraint_models.cbegin(), constraint_models.cend());

      ::pinocchio::internal::computeConeProjection(wrapped_constraint_models, forces, res);
      return res;
    }

    template<typename ConstraintModel, typename ConstraintModelAllocator>
    static context::VectorXs computeDualConeProjection_wrapper(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const context::VectorXs & velocities)
    {
      context::VectorXs res(velocities.size());

      typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
      typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

      WrappedConstraintModelVector wrapped_constraint_models(
        constraint_models.cbegin(), constraint_models.cend());

      ::pinocchio::internal::computeDualConeProjection(wrapped_constraint_models, velocities, res);
      return res;
    }

    template<typename ConstraintModel, typename ConstraintModelAllocator>
    static context::Scalar computePrimalFeasibility_wrapper(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const context::VectorXs & forces)
    {
      return ::pinocchio::internal::computePrimalFeasibility(constraint_models, forces);
    }

    template<typename ConstraintModel, typename ConstraintModelAllocator>
    static context::Scalar computeReprojectionError_wrapper(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const context::VectorXs & forces,
      const context::VectorXs & velocities)
    {
      return ::pinocchio::internal::computeReprojectionError(constraint_models, forces, velocities);
    }

    template<typename ConstraintModel, typename ConstraintModelAllocator>
    static context::VectorXs computeComplementarityShift_wrapper(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const context::VectorXs & velocities)
    {
      context::VectorXs res(velocities.size());
      typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
      typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

      WrappedConstraintModelVector wrapped_constraint_models(
        constraint_models.cbegin(), constraint_models.cend());

      ::pinocchio::internal::computeComplementarityShift(
        wrapped_constraint_models, velocities, res);
      return res;
    }
#endif // PINOCCHIO_PYTHON_SKIP_CASADI_UNSUPPORTED

    struct SolveMethodExposer
    {
      SolveMethodExposer(bp::class_<Solver> & class_)
      : class_(class_)
      {
      }

      template<class T>
      void operator()(T)
      {
        run(static_cast<typename T::type *>(nullptr));
      }

      template<typename ConstraintModel>
      void run(ConstraintModelBase<ConstraintModel> * ptr = 0)
      {
        PINOCCHIO_UNUSED_VARIABLE(ptr);
        typedef Eigen::aligned_allocator<ConstraintModel> ConstraintModelAllocator;

        if (!eigenpy::register_symbolic_link_to_registered_type<::pinocchio::ADMMUpdateRule>())
        {
          bp::enum_<::pinocchio::ADMMUpdateRule>("ADMMUpdateRule")
            .value("SPECTRAL", ::pinocchio::ADMMUpdateRule::SPECTRAL)
            .value("LINEAR", ::pinocchio::ADMMUpdateRule::LINEAR)
            // .export_values()
            ;
        }

        class_
          .def(
            "solve",
            solve_wrapper<
              ContactCholeskyDecomposition::DelassusCholeskyExpression, ConstraintModel,
              ConstraintModelAllocator>,
            (bp::args("self", "delassus", "g", "constraint_models", "R"),
             bp::arg("primal_solution") = boost::none, bp::arg("dual_solution") = boost::none,
             bp::arg("solve_ncp") = true, bp::arg("compute_largest_eigen_values") = true,
             bp::arg("admm_update_rule") = ADMMUpdateRule::SPECTRAL,
             bp::arg("stat_record") = false),
            "Solve the constrained conic problem, starting from the optional initial guess.")
          .def(
            "solve",
            solve_wrapper<
              context::DelassusOperatorDense, ConstraintModel, ConstraintModelAllocator>,
            (bp::args("self", "delassus", "g", "constraint_models", "R"),
             bp::arg("primal_solution") = boost::none, bp::arg("dual_solution") = boost::none,
             bp::arg("solve_ncp") = true, bp::arg("compute_largest_eigen_values") = true,
             bp::arg("admm_update_rule") = ADMMUpdateRule::SPECTRAL,
             bp::arg("stat_record") = false),
            "Solve the constrained conic problem, starting from the optional initial guess.")
          .def(
            "solve",
            solve_wrapper<
              context::DelassusOperatorSparse, ConstraintModel, ConstraintModelAllocator>,
            (bp::args("self", "delassus", "g", "constraint_models", "R"),
             bp::arg("primal_solution") = boost::none, bp::arg("dual_solution") = boost::none,
             bp::arg("solve_ncp") = true, bp::arg("compute_largest_eigen_values") = true,
             bp::arg("admm_update_rule") = ADMMUpdateRule::SPECTRAL,
             bp::arg("stat_record") = false),
            "Solve the constrained conic problem, starting from the optional initial guess.");
#ifdef PINOCCHIO_WITH_ACCELERATE_SUPPORT
        {
          typedef Eigen::AccelerateLLT<context::SparseMatrix> AccelerateLLT;
          typedef DelassusOperatorSparseTpl<context::Scalar, context::Options, AccelerateLLT>
            DelassusOperatorSparseAccelerate;
          class_.def(
            "solve",
            solve_wrapper<
              DelassusOperatorSparseAccelerate, ConstraintModel, ConstraintModelAllocator>,
            (bp::args("self", "delassus", "g", "constraint_models", "R"),
             bp::arg("primal_solution") = boost::none, bp::arg("dual_solution") = boost::none,
             bp::arg("compute_largest_eigen_values") = true, bp::arg("solve_ncp") = true,
             bp::arg("admm_update_rule") = pinocchioADMMUpdateRule::SPECTRAL,
             bp::arg("stat_record") = false),
            "Solve the constrained conic problem, starting from the optional initial guess.");
        }
#endif

        bp::def(
          "computeConeProjection",
          computeConeProjection_wrapper<ConstraintModel, ConstraintModelAllocator>,
          bp::args("constraint_models", "forces"),
          "Project a vector on the cartesian product of the constraint set associated with each "
          "constraint model.");

        bp::def(
          "computeDualConeProjection",
          computeDualConeProjection_wrapper<ConstraintModel, ConstraintModelAllocator>,
          bp::args("constraint_models", "velocities"),
          "Project a vector on the cartesian product of dual cones.");

        // TODO(jcarpent): restore these two next signatures
        //        bp::def(
        //                "computePrimalFeasibility", computePrimalFeasibility_wrapper,
        //                bp::args("constraint_models", "forces"), "Compute the primal
        //                feasibility.");

        //        bp::def(
        //                "computeReprojectionError", computeReprojectionError_wrapper,
        //                bp::args("constraint_models", "forces", "velocities"), "Compute the
        //                reprojection error.");

        bp::def(
          "computeComplementarityShift",
          computeComplementarityShift_wrapper<ConstraintModel, ConstraintModelAllocator>,
          bp::args("constraint_models", "velocities"),
          "Compute the complementarity shift associated to the De Saxé function.");
      }
      //
      //      template<typename S, int O>
      //      void run(FictiousConstraintModelTpl<S, O> * ptr = 0)
      //      {
      //        PINOCCHIO_UNUSED_VARIABLE(ptr);
      //      }
      //
      void run(boost::blank * ptr = 0)
      {
        PINOCCHIO_UNUSED_VARIABLE(ptr);
      }

      bp::class_<Solver> & class_;
    };

    template<typename ConstraintModel>
    static void expose_solve(bp::class_<Solver> & class_)
    {
      SolveMethodExposer expose(class_);
      expose.run(static_cast<ConstraintModel *>(nullptr));
    }

    void exposeADMMContactSolver()
    {
#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE

      //      bp::enum_<::pinocchio::ADMMUpdateRule>("ADMMUpdateRule")
      //        .value("SPECTRAL", ::pinocchio::ADMMUpdateRule::SPECTRAL)
      //        .value("LINEAR", ::pinocchio::ADMMUpdateRule::LINEAR)
      //        // .export_values()
      //        ;

      bp::class_<Solver> cl(
        "ADMMContactSolver",
        "Alternating Direction Method of Multi-pliers solver for contact dynamics.",
        bp::init<int, Scalar, Scalar, Scalar, Scalar, Scalar, Scalar, int>(
          (bp::arg("self"), bp::arg("problem_dim"), bp::arg("mu_prox") = Scalar(1e-6),
           bp::arg("tau") = Scalar(0.5), bp::arg("rho_power") = Scalar(0.2),
           bp::arg("rho_power_factor") = Scalar(0.05),
           bp::arg("linear_update_rule_factor") = Scalar(10),
           bp::arg("ratio_primal_dual") = Scalar(10),
           bp::arg("max_it_largest_eigen_value_solver") = 20),
          "Default constructor."));

      cl.def(ContactSolverBasePythonVisitor<Solver>())

        .def("setRho", &Solver::setRho, bp::args("self", "rho"), "Set the ADMM penalty value.")
        .def("getRho", &Solver::getRho, bp::arg("self"), "Get the ADMM penalty value.")

        .def(
          "setRhoPower", &Solver::setRhoPower, bp::args("self", "rho_power"),
          "Set the power associated to the ADMM spectral update rule.")
        .def(
          "getRhoPower", &Solver::getRhoPower, bp::arg("self"),
          "Get the power associated to the ADMM spectral update rule.")

        .def(
          "setRhoPowerFactor", &Solver::setRhoPowerFactor, bp::args("self", "rho_power_factor"),
          "Set the power factor associated to the ADMM spectral update rule.")
        .def(
          "getRhoPowerFactor", &Solver::getRhoPowerFactor, bp::arg("self"),
          "Get the power factor associated to the ADMM spectral update rule.")

        .def(
          "setLinearUpdateRuleFactor", &Solver::setLinearUpdateRuleFactor,
          bp::args("self", "linear_update_rule_factor"),
          "Set the factor associated with the ADMM linear update rule.")
        .def(
          "getLinearUpdateRuleFactor", &Solver::getLinearUpdateRuleFactor, bp::arg("self"),
          "Get the factor associated with the ADMM linear update rule.")

        .def(
          "setTau", &Solver::setTau, bp::args("self", "tau"), "Set the tau linear scaling factor.")
        .def("getTau", &Solver::getTau, bp::arg("self"), "Get the tau linear scaling factor.")

        .def(
          "setProximalValue", &Solver::setProximalValue, bp::args("self", "mu"),
          "Set the proximal value.")
        .def(
          "getProximalValue", &Solver::getProximalValue, bp::arg("self"), "Get the proximal value.")

        .def(
          "setRatioPrimalDual", &Solver::setRatioPrimalDual, bp::args("self", "ratio_primal_dual"),
          "Set the primal/dual ratio.")
        .def(
          "getRatioPrimalDual", &Solver::getRatioPrimalDual, bp::arg("self"),
          "Get the primal/dual ratio.")

        .def(
          "getPrimalSolution", &Solver::getPrimalSolution, bp::arg("self"),
          "Returns the primal solution of the problem.", bp::return_internal_reference<>())

        .def(
          "getDualSolution", &Solver::getDualSolution, bp::arg("self"),
          "Returns the dual solution of the problem.", bp::return_internal_reference<>())

        .def(
          "getCholeskyUpdateCount", &Solver::getCholeskyUpdateCount, bp::arg("self"),
          "Returns the number of updates of the Cholesky factorization due to rho updates.")

        //    .def("computeRho",&Solver::computeRho,bp::args("L","m","rho_power"),
        //         "Compute the penalty ADMM value from the current largest and lowest eigenvalues
        //         and the scaling spectral factor.")
        //    .staticmethod("computeRho")
        //    .def("computeRhoPower",&Solver::computeRhoPower,bp::args("L","m","rho"),
        //         "Compute the  scaling spectral factor of the ADMM penalty term from the current
        //         largest and lowest eigenvalues and the ADMM penalty term.")
        //    .staticmethod("computeRhoPower")

        .def("getStats", &Solver::getStats, bp::arg("self"), bp::return_internal_reference<>());

      typedef context::ConstraintModel::ConstraintModelVariant ConstraintModelVariant;

      SolveMethodExposer solve_exposer(cl);
      boost::mpl::for_each<
        ConstraintModelVariant::types, boost::mpl::make_identity<boost::mpl::_1>>(solve_exposer);
      expose_solve<context::ConstraintModel>(cl);

      {
        bp::class_<SolverStats>(
          "SolverStats", "",
          bp::init<int>((bp::arg("self"), bp::arg("max_it")), "Default constructor"))
          .def("reset", &SolverStats::reset, bp::arg("self"), "Reset the stasts.")
          .def(
            "size", &SolverStats::size, bp::arg("self"),
            "Size of the vectors stored in the structure.")

          .PINOCCHIO_ADD_PROPERTY_READONLY(SolverStats, primal_feasibility, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(SolverStats, dual_feasibility, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(SolverStats, dual_feasibility_ncp, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(SolverStats, complementarity, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(SolverStats, rho, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(
            SolverStats, it, "Number of iterations performed by the algorithm.")
          .PINOCCHIO_ADD_PROPERTY_READONLY(
            SolverStats, cholesky_update_count,
            "Number of Cholesky updates performed by the algorithm.");
      }

      {
        typedef PowerIterationAlgoTpl<context::VectorXs> PowerIterationAlgo;
        bp::class_<PowerIterationAlgo>(
          "PowerIterationAlgo", "",
          bp::init<Eigen::DenseIndex, int, Scalar>(
            (bp::arg("self"), bp::arg("size"), bp::arg("max_it") = 10, bp::arg("rel_tol") = 1e-8),
            "Default constructor"))
          .def("run", &PowerIterationAlgo::run<context::MatrixXs>, bp::arg("self"), "")
          .def(
            "lowest", &PowerIterationAlgo::lowest<context::MatrixXs>,
            (bp::arg("self"), bp::arg("compute_largest") = true), "")
          .def("reset", &PowerIterationAlgo::reset, bp::arg("self"))

          .PINOCCHIO_ADD_PROPERTY(PowerIterationAlgo, principal_eigen_vector, "")
          .PINOCCHIO_ADD_PROPERTY(PowerIterationAlgo, lowest_eigen_vector, "")
          .PINOCCHIO_ADD_PROPERTY(PowerIterationAlgo, max_it, "")
          .PINOCCHIO_ADD_PROPERTY(PowerIterationAlgo, rel_tol, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(PowerIterationAlgo, largest_eigen_value, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(PowerIterationAlgo, lowest_eigen_value, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(PowerIterationAlgo, it, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(PowerIterationAlgo, convergence_criteria, "");
      }
#endif
    }

  } // namespace python
} // namespace pinocchio
