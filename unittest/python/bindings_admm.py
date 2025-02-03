import unittest
from pathlib import Path

import numpy as np
import pinocchio as pin
from test_case import PinocchioTestCase as TestCase

import importlib.util

coal_spec = importlib.util.find_spec("coal")
coal_found = coal_spec is not None
if coal_found:
    import coal

meshcat_spec = importlib.util.find_spec("meshcat")
meshcat_found = meshcat_spec is not None
if meshcat_found:
    import meshcat
    from pinocchio.visualize import MeshcatVisualizer


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
        list_cm = []
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
        constraint_datas = []
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

    @unittest.skipUnless(coal_found, "Needs Coal.")
    def addFloor(self, geom_model: pin.GeometryModel, visual_model: pin.GeometryModel):
        floor_collision_shape = coal.Halfspace(0, 0, 1, 0)
        M = pin.SE3.Identity()
        floor_collision_object = pin.GeometryObject(
            "floor", 0, 0, M, floor_collision_shape
        )
        floor_collision_object.meshColor = np.array([0.5, 0.5, 0.5, 0.5])
        geom_model.addGeometryObject(floor_collision_object)

        h = 0.01
        floor_visual_shape = coal.Box(20, 20, h)
        Mvis = pin.SE3.Identity()
        Mvis.translation = np.array([0.0, 0.0, -h / 2])
        floor_visual_object = pin.GeometryObject(
            "floor", 0, 0, Mvis, floor_visual_shape
        )
        floor_visual_object.meshColor = np.array([0.5, 0.5, 0.5, 0.4])
        visual_model.addGeometryObject(floor_visual_object)

    @unittest.skipUnless(coal_found, "Needs Coal.")
    def addSystemCollisionPairs(self, model, geom_model, qref):
        """
        Add the right collision pairs of a model, given qref.
        qref is here as a `T-pose`. The function uses this pose to determine which objects are in collision
        in this ref pose. If objects are in collision, they are not added as collision pairs, as they are considered
        to always be in collision.
        """
        data = model.createData()
        geom_data = geom_model.createData()
        pin.updateGeometryPlacements(model, data, geom_model, geom_data, qref)
        geom_model.removeAllCollisionPairs()
        num_col_pairs = 0
        for i in range(len(geom_model.geometryObjects)):
            for j in range(i + 1, len(geom_model.geometryObjects)):
                # Don't add collision pair if same object
                if i != j:
                    gobj_i: pin.GeometryObject = geom_model.geometryObjects[i]
                    gobj_j: pin.GeometryObject = geom_model.geometryObjects[j]
                    if gobj_i.name == "floor" or gobj_j.name == "floor":
                        num_col_pairs += 1
                        col_pair = pin.CollisionPair(i, j)
                        geom_model.addCollisionPair(col_pair)
                    else:
                        if gobj_i.parentJoint != gobj_j.parentJoint:
                            # Compute collision between the geometries. Only add the collision pair if there is no collision.
                            M1 = geom_data.oMg[i]
                            M2 = geom_data.oMg[j]
                            colreq = coal.CollisionRequest()
                            colreq.security_margin = 1e-2  # 1cm of clearance
                            colres = coal.CollisionResult()
                            coal.collide(
                                gobj_i.geometry, M1, gobj_j.geometry, M2, colreq, colres
                            )
                            if not colres.isCollision():
                                num_col_pairs += 1
                                col_pair = pin.CollisionPair(i, j)
                                geom_model.addCollisionPair(col_pair)

    def complete_orthonormal_basis(self, ez, joint_placement):
        ex = joint_placement.rotation[:, 0]
        ey = np.cross(ez, ex)
        if np.linalg.norm(ey) < 1e-6:
            ex = joint_placement.rotation[:, 1]
            ey = np.cross(ez, ex)
        ey /= np.linalg.norm(ey)
        ex = np.cross(ey, ez)
        return ex, ey

    @unittest.skipUnless(coal_found, "Needs Coal.")
    def computeContactConstraints(self, model, geom_model, q):
        data = model.createData()
        geom_data = geom_model.createData()
        pin.updateGeometryPlacements(model, data, geom_model, geom_data, q)
        pin.computeCollisions(geom_model, geom_data, False)
        contact_constraints = pin.StdVec_ConstraintModel()
        for i, res in enumerate(geom_data.collisionResults):
            if res.isCollision():
                geom_id1, geom_id2 = (
                    geom_model.collisionPairs[i].first,
                    geom_model.collisionPairs[i].second,
                )
                joint_id1 = geom_model.geometryObjects[geom_id1].parentJoint
                joint_id2 = geom_model.geometryObjects[geom_id2].parentJoint
                contacts = res.getContacts()
                joint_placement_1 = data.oMi[joint_id1]
                joint_placement_2 = data.oMi[joint_id2]
                for contact in contacts:
                    pos_i = contact.pos
                    normal_i = contact.normal
                    ex_i, ey_i = self.complete_orthonormal_basis(
                        contact.normal, joint_placement_1
                    )
                    ex_i = np.expand_dims(ex_i, axis=1)
                    ey_i = np.expand_dims(ey_i, axis=1)
                    normal_i = np.expand_dims(contact.normal, axis=1)
                    R_i = np.concatenate((ex_i, ey_i, normal_i), axis=1)
                    R_i1 = np.dot(joint_placement_1.rotation.T, R_i)
                    R_i2 = np.dot(joint_placement_2.rotation.T, R_i)
                    pos_i1 = joint_placement_1.rotation.T @ (
                        pos_i - joint_placement_1.translation
                    )
                    pos_i2 = joint_placement_2.rotation.T @ (
                        pos_i - joint_placement_2.translation
                    )
                    placement_i1 = pin.SE3(R_i1, pos_i1)
                    placement_i2 = pin.SE3(R_i2, pos_i2)
                    contact_model_i = pin.FrictionalPointConstraintModel(
                        model, joint_id2, placement_i2, joint_id1, placement_i1
                    )
                    contact_constraints.append(contact_model_i)
        return contact_constraints

    def createVisualizer(
        self,
        model: pin.GeometryModel,
        geom_model: pin.GeometryModel,
        visual_model: pin.GeometryModel,
    ):
        viewer = meshcat.Visualizer(zmq_url="tcp://127.0.0.1:6000")
        viewer.delete()
        for obj in visual_model.geometryObjects:
            color = np.random.rand(4)
            color[3] = 1.0
            obj.meshColor = color
        vizer: MeshcatVisualizer = MeshcatVisualizer(model, geom_model, visual_model)
        vizer.initViewer(viewer=viewer, open=False, loadModel=True)
        return vizer, viewer

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

    @unittest.skipUnless(coal_found, "Needs Coal.")
    def test_cassie(self, display=False):
        current_dir = Path(__file__).parent
        model_dir = current_dir / "../models/"
        model_path = model_dir / "closed_chain.xml"
        constraint_models = pin.StdVec_ConstraintModel()

        # Parsing model, constraint models (bilateral constraints) and geometry model from xml description
        model, bilateral_constraint_models, geom_model, visual_model = (
            pin.buildModelsFromMJCF(model_path)
        )

        # Adding bilateral constraints to the list of constraints
        for bpcm in bilateral_constraint_models:
            constraint_models.append(pin.ConstraintModel(bpcm))

        # adding joint limit constraints
        active_joints_limits = [i for i in range(1, model.njoints)]
        jlcm = pin.JointLimitConstraintModel(model, active_joints_limits)
        constraint_models.append(pin.ConstraintModel(jlcm))

        # adding friction on joints
        active_joints_friction = [i for i in range(1, model.njoints)]
        fjcm = pin.FrictionalJointConstraintModel(model, active_joints_friction)
        constraint_models.append(pin.ConstraintModel(fjcm))

        data = model.createData()
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

        delassus, g = self.setupTest(model, constraint_models, q0, v0, tau0, fext, dt)

        compliance = np.zeros_like(g)
        dim_pb = g.shape[0]
        solver = pin.ADMMContactSolver(dim_pb)
        solver.setAbsolutePrecision(1e-13)
        solver.setRelativePrecision(1e-14)

        if display:
            vizer, viewer = self.createVisualizer(model, geom_model, geom_model)
            vizer.display(q0)


if __name__ == "__main__":
    unittest.main()
