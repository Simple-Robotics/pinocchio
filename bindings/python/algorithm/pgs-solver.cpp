//
// Copyright (c) 2022-2025 INRIA
//

#include "pinocchio/algorithm/pgs-solver.hpp"
#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/bindings/python/fwd.hpp"

#include "pinocchio/bindings/python/algorithm/contact-solver-base.hpp"
#include "pinocchio/bindings/python/utils/std-vector.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"
#include <eigenpy/eigen-from-python.hpp>

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    typedef PGSContactSolverTpl<context::Scalar> Solver;
    typedef Solver::SolverStats SolverStats;
    typedef typename Solver::RefConstVectorXs RefConstVectorXs;

#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
    template<typename DelassusMatrixType, typename ConstraintModel>
    static bool solve_wrapper(
      Solver & solver,
      const DelassusOperatorDense & delassus,
      const context::VectorXs & g,
      const context::ConstraintModelVector & constraint_models,
      const context::Scalar dt,
      const boost::optional<RefConstVectorXs> x = boost::none,
      const context::Scalar over_relax = 1,
      const bool solve_ncp = true,
      const bool stat_record = false)
    {
      return solver.solve(
        delassus, g, constraint_models, dt, x, over_relax, solve_ncp, stat_record);
    }
#endif

    template<typename Solver>
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
#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
        class_
          .def(
            "solve", solve_wrapper<context::MatrixXs, ConstraintModel>,
            (bp::args("self", "delassus", "g", "constraint_models", "dt"),
             bp::arg("primal_solution") = boost::none, bp::arg("over_relax") = context::Scalar(1),
             bp::arg("solve_ncp") = true, bp::arg("stat_record") = false),
            "Solve the constrained conic problem composed of problem data (G,g,cones) and starting "
            "from the initial guess.")
          .def(
            "solve", solve_wrapper<context::SparseMatrix, ConstraintModel>,
            (bp::args("self", "delassus", "g", "constraint_models", "dt"),
             bp::arg("primal_solution") = boost::none, bp::arg("over_relax") = context::Scalar(1),
             bp::arg("solve_ncp") = true, bp::arg("stat_record") = false),
            "Solve the constrained conic problem composed of problem data (G,g,cones) and starting "
            "from the initial guess.");
#endif // ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
      }

      //      template<typename S, int O>
      //      void run(FictiousConstraintModelTpl<S, O> * ptr = 0)
      //      {
      //        PINOCCHIO_UNUSED_VARIABLE(ptr);
      //      }

      void run(boost::blank * ptr = 0)
      {
        PINOCCHIO_UNUSED_VARIABLE(ptr);
      }

      bp::class_<Solver> & class_;
    };

    template<typename ConstraintModel>
    static void expose_solve(bp::class_<Solver> & class_)
    {
      SolveMethodExposer<Solver> expose(class_);
      expose.run(static_cast<ConstraintModel *>(nullptr));
    }

    void exposePGSContactSolver()
    {
#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
      bp::class_<Solver> class_(
        "PGSContactSolver", "Projected Gauss Siedel solver for contact dynamics.",
        bp::init<int>(bp::args("self", "problem_dim"), "Default constructor."));
      class_.def(ContactSolverBasePythonVisitor<Solver>())
        .def(
          "getPrimalSolution", &Solver::getPrimalSolution, bp::arg("self"),
          "Returns the primal solution of the problem.", bp::return_internal_reference<>())
        .def(
          "getDualSolution", &Solver::getDualSolution, bp::arg("self"),
          "Returns the dual solution of the problem.", bp::return_internal_reference<>())

        .def("getStats", &Solver::getStats, bp::arg("self"), bp::return_internal_reference<>());

      //      typedef context::ConstraintModel::ConstraintModelVariant ConstraintModelVariant;
      //
      //       SolveMethodExposer<Solver> solve_exposer(class_);
      //       boost::mpl::for_each<
      //         ConstraintModelVariant::types,
      //         boost::mpl::make_identity<boost::mpl::_1>>(solve_exposer);
      expose_solve<context::ConstraintModel>(class_);

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
          .PINOCCHIO_ADD_PROPERTY_READONLY(
            SolverStats, it, "Number of iterations performed by the algorithm.");
      }

#endif
    }

  } // namespace python
} // namespace pinocchio
