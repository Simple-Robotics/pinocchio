//
// Copyright (c) 2022-2024 INRIA
//

#include "pinocchio/algorithm/pgs-solver.hpp"
#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/bindings/python/fwd.hpp"

#include "pinocchio/bindings/python/algorithm/contact-solver-base.hpp"
#include "pinocchio/bindings/python/utils/std-vector.hpp"
#include <eigenpy/eigen-from-python.hpp>

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    typedef PGSContactSolverTpl<context::Scalar> Solver;

#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
    template<typename DelassusMatrixType, typename ConstraintModel>
    static bool solve_wrapper(
      Solver & solver,
      const DelassusMatrixType & G,
      const context::VectorXs & g,
      const PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(ConstraintModel) & constraint_models,
      Eigen::Ref<context::VectorXs> x,
      const context::Scalar over_relax = 1)
    {
      return solver.solve(G, g, constraint_models, x, over_relax);
    }
#endif

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
        class_
          .def(
            "solve", solve_wrapper<context::MatrixXs, ConstraintModel>,
            (bp::args("self", "G", "g", "constraint_sets", "x"),
             (bp::arg("over_relax") = context::Scalar(1))),
            "Solve the constrained conic problem composed of problem data (G,g,cones) and starting "
            "from the initial guess.")
          .def(
            "solve", solve_wrapper<context::SparseMatrix, ConstraintModel>,
            (bp::args("self", "G", "g", "constraint_sets", "x"),
             (bp::arg("over_relax") = context::Scalar(1))),
            "Solve the constrained conic problem composed of problem data (G,g,cones) and starting "
            "from the initial guess.");
      }

      template<typename S, int O>
      void run(FictiousConstraintModelTpl<S, O> * ptr = 0)
      {
        PINOCCHIO_UNUSED_VARIABLE(ptr);
      }

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
          "Returns the dual solution of the problem.", bp::return_internal_reference<>());

      typedef context::ConstraintModel::ConstraintModelVariant ConstraintModelVariant;

      SolveMethodExposer solve_exposer(class_);
      boost::mpl::for_each<
        ConstraintModelVariant::types, boost::mpl::make_identity<boost::mpl::_1>>(solve_exposer);
      expose_solve<context::ConstraintModel>(class_);

#endif
    }

  } // namespace python
} // namespace pinocchio
