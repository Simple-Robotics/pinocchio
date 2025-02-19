import unittest
from pathlib import Path

import numpy as np
import pinocchio as pin
from test_case import ContactSolverTestCase as TestCase

import importlib.util

coal_spec = importlib.util.find_spec("coal")
coal_found = coal_spec is not None

matplotlib_spec = importlib.util.find_spec("matplotlib")
matplotlib_found = matplotlib_spec is not None

meshcat_spec = importlib.util.find_spec("meshcat")
meshcat_found = meshcat_spec is not None


class TestADMM(TestCase):
    def test_box(self):
        model, constraint_models = self.buildStackOfCubesModel([1e-3])
        q0 = pin.neutral(model)
        v0 = np.zeros(model.nv)
        tau0 = np.zeros(model.nv)
        fext = [pin.Force.Zero() for i in range(model.njoints)]
        dt = 1e-3
        delassus_matrix, g = self.setupTest(
            model, constraint_models, q0, v0, tau0, fext, dt
        )
        delassus = pin.DelassusOperatorDense(delassus_matrix)
        dim_pb = g.shape[0]
        solver = pin.ADMMContactSolver(dim_pb)
        solver.setAbsolutePrecision(1e-13)
        solver.setRelativePrecision(1e-14)
        solver.setLanczosSize(g.size)
        solver.solve(delassus, g, constraint_models, dt)

    @unittest.skipUnless(coal_found, "Needs Coal.")
    def test_cassie(self, display=False, stat_record=True):
        current_dir = Path(__file__).parent
        model_dir = current_dir / "../models/"
        model_path = model_dir / "closed_chain.xml"
        constraint_models = pin.StdVec_ConstraintModel()

        # Parsing model, constraint models and geometry model from xml description
        model, constraint_models_dict, geom_model, visual_model = (
            pin.buildModelsFromMJCF(model_path)
        )

        # Adding all constraintds would be
        for typed_constraint_models in constraint_models_dict.values():
            for cm in typed_constraint_models:
                constraint_models.append(pin.ConstraintModel(cm))
        # Adding only bilateral constraints to the list of constraints
        # for bpcm in constraint_models_dict['bilateral_point_constraint_models']:
        #     constraint_models.append(pin.ConstraintModel(bpcm))

        # adding joint limit constraints
        active_joints_limits = [i for i in range(1, model.njoints)]
        jlcm = pin.JointLimitConstraintModel(model, active_joints_limits)
        constraint_models.append(pin.ConstraintModel(jlcm))

        # adding friction on joints
        active_joints_friction = [i for i in range(1, model.njoints)]
        fjcm = pin.FrictionalJointConstraintModel(model, active_joints_friction)
        fjcm.set = pin.BoxSet(model.lowerDryFrictionLimit, model.upperDryFrictionLimit)
        constraint_models.append(pin.ConstraintModel(fjcm))

        q0 = model.referenceConfigurations["home"]
        v0 = np.zeros(model.nv)
        tau0 = np.zeros(model.nv)
        fext = [pin.Force.Zero() for i in range(model.njoints)]
        dt = 1e-3
        self.addFloor(geom_model, visual_model)
        self.addSystemCollisionPairs(model, geom_model, q0)

        # Adding constraints from frictional contacts
        contact_constraints = self.computeContactConstraints(model, geom_model, q0)
        for fpcm in contact_constraints:
            constraint_models.append(pin.ConstraintModel(fpcm))

        delassus_matrix, g = self.setupTest(
            model, constraint_models, q0, v0, tau0, fext, dt
        )
        delassus = pin.DelassusOperatorDense(delassus_matrix)

        self.assertTrue(
            delassus.matrix().shape[0]
            == (3 * len(constraint_models_dict["bilateral_point_constraint_models"]))
            + (6 * len(constraint_models_dict["weld_constraint_models"]))
            + (3 * len(contact_constraints))
            + (model.upperPositionLimit != np.inf).sum()
            - 4 * 3
            + (model.lowerPositionLimit != -np.inf).sum()
            - 4 * 3
            + model.nv,
            "constraint problem is of wrong size.",
        )

        dim_pb = g.shape[0]
        solver = pin.ADMMContactSolver(dim_pb)
        solver.setAbsolutePrecision(1e-13)
        solver.setRelativePrecision(1e-14)
        solver.setLanczosSize(g.size)

        has_converged = solver.solve(
            delassus,
            g,
            constraint_models,
            dt,
            None,
            None,
            None,
            True,
            pin.ADMMUpdateRule.SPECTRAL,
            stat_record,
        )
        self.assertTrue(has_converged, "Solver did not converge.")
        print(solver.getIterationCount())
        print(solver.getAbsoluteConvergenceResidual())
        print(solver.getRelativeConvergenceResidual())

        if stat_record and matplotlib_found:
            self.plotContactSolver(solver)

        if display and meshcat_found:
            vizer, viewer = self.createVisualizer(model, geom_model, geom_model)
            vizer.display(q0)


if __name__ == "__main__":
    unittest.main()
