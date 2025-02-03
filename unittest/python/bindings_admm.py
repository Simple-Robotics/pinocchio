import unittest

import numpy as np
import pinocchio as pin
from test_case import PinocchioTestCase as TestCase


class TestADMM(TestCase):
    def buildStackOfCubesModel(self, masses):
        model = pin.Model()
        n_cubes = len(masses)
        box_dims = np.ones((3, 1))
        for i in range(n_cubes):
            box_mass = masses[i]
            box_inertia = pin.Inertia.FromBox(
                box_mass, box_dims[0, 0], box_dims[1, 0], box_dims[2, 0]
            )
            joint_id = model.addJoint(
                0, pin.JointModelFreeFlyer(), pin.SE3.Identity(), "free_flyer_" + str(i)
            )
            model.appendBodyToJoint(joint_id, box_inertia, pin.SE3.Identity())

        friction_coeff = 0.4
        list_cm = pin.StdVec_ConstraintModel()
        for i in range(n_cubes):
            local_trans_box_1 = 0.5 * box_dims
            local_trans_box_2 = 0.5 * box_dims
            local_trans_box_2[2] *= -1
            rot = np.identity(3)
            theta = np.pi / 2
            c, s = np.cos(theta), np.sin(theta)
            R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
            for j in range(4):
                local_placement_1 = pin.SE3(np.identity(3), (rot @ local_trans_box_1))
                local_placement_2 = pin.SE3(np.identity(3), (rot @ local_trans_box_2))
                fpcm = pin.FrictionalPointConstraintModel(
                    model, i, local_placement_1, i + 1, local_placement_2
                )
                fpcm.set = pin.CoulombFrictionCone(friction_coeff)
                list_cm.append(pin.ConstraintModel(fpcm))
                rot = R @ rot

        return model, list_cm

    def setupTest(self, model, constraint_models, q0, v0, tau0, fext, dt):
        data = model.createData()
        constraint_datas = pin.StdVec_ConstraintData()
        for cm in constraint_models:
            constraint_datas.append(cm.createData())
        pin.crba(model, data, q0, pin.Convention.WORLD)
        chol = pin.ContactCholeskyDecomposition(model, constraint_models)
        chol.compute(model, data, constraint_models, constraint_datas, 1e-10)
        delassus_matrix = chol.getDelassusCholeskyExpression().matrix()
        delassus = pin.DelassusOperatorDense(delassus_matrix)
        vfree = v0 + dt * pin.aba(model, data, q0, v0, tau0, fext)
        Jc = pin.getConstraintsJacobian(
            model, data, constraint_models, constraint_datas
        )
        g = Jc @ vfree
        return delassus, g

    def test_box(self):
        model, constraint_models = self.buildStackOfCubesModel([1e-3])
        q0 = pin.neutral(model)
        v0 = np.zeros(model.nv)
        tau0 = np.zeros(model.nv)
        fext = [pin.Force.Zero() for i in range(model.njoints)]
        dt = 1e-3
        delassus, g = self.setupTest(model, constraint_models, q0, v0, tau0, fext, dt)
        compliance = np.zeros_like(g)
        dim_pb = g.shape[0]
        solver = pin.ADMMContactSolver(dim_pb)
        solver.setAbsolutePrecision(1e-13)
        solver.setRelativePrecision(1e-14)
        solver.solve(delassus, g, constraint_models, dt, compliance)


if __name__ == "__main__":
    unittest.main()
