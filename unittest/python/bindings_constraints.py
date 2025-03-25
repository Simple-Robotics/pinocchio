import itertools
import unittest

import coal as coal
import numpy as np
import pinocchio as pin
from test_case import PinocchioTestCase as TestCase


class TestJointsAlgo(TestCase):
    def setUp(self):
        self.construct_model()
        self.construct_constraints()

    def construct_model(self):
        model = pin.Model()
        geom_model = pin.GeometryModel()

        world_joint_id = 0

        # Add a plane in world model
        plane_geom_name = "geom_plane"
        plane_geom_shape = coal.Halfspace(0.0, 0.0, 1.0, 0.0)
        plane_geom_placement = pin.SE3.Identity()
        plane_geom_obj = pin.GeometryObject(
            plane_geom_name, world_joint_id, plane_geom_placement, plane_geom_shape
        )
        plane_geom_id = geom_model.addGeometryObject(plane_geom_obj)

        # Add a freeflyer with a ball attach to it
        ball_joint_name = "joint_ball_ff"
        ball_joint_obj = pin.JointModelFreeFlyer()
        ball_joint_placement = pin.SE3.Identity()
        ball_joint_id = model.addJoint(
            world_joint_id, ball_joint_obj, ball_joint_placement, ball_joint_name
        )
        ball_geom_radius = 1.0
        ball_geom_name = "geom_ball_ff"
        ball_geom_shape = coal.Sphere(ball_geom_radius)
        ball_geom_placement = pin.SE3.Identity()
        ball_geom_obj = pin.GeometryObject(
            ball_geom_name, ball_joint_id, ball_geom_placement, ball_geom_shape
        )
        ball_geom_id = geom_model.addGeometryObject(ball_geom_obj)
        ball_body_mass = 1.0
        ball_body_placement = pin.SE3.Identity()
        ball_body_inertia = pin.Inertia.FromSphere(ball_body_mass, ball_geom_radius)
        model.appendBodyToJoint(
            ball_joint_id, ball_body_inertia, ball_body_placement
        )  # Add inertia of the ball to the model
        geom_model.addCollisionPair(
            pin.CollisionPair(plane_geom_id, ball_geom_id)
        )  # Add collision pair with the floor to the floor

        # Add a prismatic with a cylinder attach to it
        cyl_joint_name = "joint_cyl_PZ"
        cyl_joint_obj = pin.JointModelPZ()
        cyl_joint_placement = pin.SE3.Identity()
        cyl_joint_id = model.addJoint(
            world_joint_id, cyl_joint_obj, cyl_joint_placement, cyl_joint_name
        )
        cyl_geom_radius = 1.0
        cyl_geom_length = 1.0
        cyl_geom_name = "geom_cyl_PZ"
        cyl_geom_shape = coal.Cylinder(cyl_geom_radius, cyl_geom_length)
        cyl_geom_placement = pin.SE3.Identity()
        cyl_geom_obj = pin.GeometryObject(
            cyl_geom_name, cyl_joint_id, cyl_geom_placement, cyl_geom_shape
        )
        cyl_geom_id = geom_model.addGeometryObject(cyl_geom_obj)
        cyl_body_mass = 1.0
        cyl_body_placement = pin.SE3.Identity()
        cyl_body_inertia = pin.Inertia.FromCylinder(
            cyl_body_mass, cyl_geom_radius, cyl_geom_length
        )
        model.appendBodyToJoint(
            cyl_joint_id, cyl_body_inertia, cyl_body_placement
        )  # Add inertia of the ball to the model
        geom_model.addCollisionPair(
            pin.CollisionPair(plane_geom_id, cyl_geom_id)
        )  # Add collision pair with the floor to the floor

        # Add a revolute joint with a capsule attach to it
        caps_joint_name = "joint_caps_RX"
        caps_joint_obj = pin.JointModelRX()
        caps_joint_placement = pin.SE3.Identity()
        caps_joint_id = model.addJoint(
            world_joint_id, caps_joint_obj, caps_joint_placement, caps_joint_name
        )
        caps_geom_radius = 1.0
        caps_geom_halflength = 1.0
        caps_geom_name = "geom_caps_RX"
        caps_geom_shape = coal.Capsule(caps_geom_radius, caps_geom_halflength)
        caps_geom_placement = pin.SE3.Identity()
        caps_geom_obj = pin.GeometryObject(
            caps_geom_name, caps_joint_id, caps_geom_placement, caps_geom_shape
        )
        caps_geom_id = geom_model.addGeometryObject(caps_geom_obj)
        caps_body_mass = 1.0
        caps_body_placement = pin.SE3.Identity()
        caps_body_inertia = pin.Inertia.FromCylinder(
            caps_body_mass, caps_geom_radius, caps_geom_halflength * 2
        )
        model.appendBodyToJoint(
            caps_joint_id, caps_body_inertia, caps_body_placement
        )  # Add inertia of the ball to the model
        geom_model.addCollisionPair(pin.CollisionPair(plane_geom_id, caps_geom_id))
        geom_model.addCollisionPair(pin.CollisionPair(ball_geom_id, caps_geom_id))

        # Set up limits
        model.positionLimitMargin = np.array([1.0] * model.nq)
        model.upperPositionLimit = np.array([0.5] * model.nq)
        model.upperPositionLimit[-1] = 2.0
        model.lowerPositionLimit = np.array([-0.5] * model.nq)
        model.lowerDryFrictionLimit = -np.ones(model.nv)
        model.upperDryFrictionLimit = np.ones(model.nv)

        # Store some values
        self.model = model
        self.geom_model = geom_model
        self.world_joint_id = world_joint_id
        self.caps_joint_id = caps_joint_id
        self.cyl_joint_id = cyl_joint_id

    def construct_constraints(self):
        model = self.model
        geom_model = self.geom_model
        world_joint_id = self.world_joint_id
        caps_joint_id = self.caps_joint_id
        cyl_joint_id = self.cyl_joint_id

        # Add all constraints
        constraints_std_vec = pin.StdVec_ConstraintModel()
        constraints_list = []

        # Bilateral constraints
        bilat_joint1_anchor_placement = pin.SE3.Identity()
        bilat1 = pin.BilateralPointConstraintModel(model, caps_joint_id)
        bilat2 = pin.BilateralPointConstraintModel(
            model, caps_joint_id, bilat_joint1_anchor_placement
        )
        bilat3 = pin.BilateralPointConstraintModel(model, caps_joint_id, world_joint_id)
        bilat4 = pin.BilateralPointConstraintModel(
            model,
            caps_joint_id,
            bilat_joint1_anchor_placement,
            world_joint_id,
            pin.SE3.Identity(),
        )
        bilats_list = [bilat1, bilat2, bilat3, bilat4]
        for b in bilats_list:
            gb = pin.ConstraintModel(b)
            constraints_std_vec.append(gb)
            constraints_list.append(b)

        # weld constraints
        weld_joint1_anchor_placement = pin.SE3.Identity()
        weld1 = pin.WeldConstraintModel(model, cyl_joint_id)
        weld2 = pin.WeldConstraintModel(
            model, cyl_joint_id, weld_joint1_anchor_placement
        )
        weld3 = pin.WeldConstraintModel(model, cyl_joint_id, world_joint_id)
        weld4 = pin.WeldConstraintModel(
            model,
            cyl_joint_id,
            weld_joint1_anchor_placement,
            world_joint_id,
            pin.SE3.Identity(),
        )
        welds_list = [weld1, weld2, weld3, weld4]
        for w in welds_list:
            gw = pin.ConstraintModel(w)
            constraints_std_vec.append(gw)
            constraints_list.append(w)

        # Joint friction
        jfc = pin.FrictionalJointConstraintModel(model, [1, 3])
        jfc.set = pin.BoxSet(
            np.array([model.lowerDryFrictionLimit[i] for i in jfc.getActiveDofs()]),
            np.array([model.upperDryFrictionLimit[i] for i in jfc.getActiveDofs()]),
        )
        constraints_std_vec.append(pin.ConstraintModel(jfc))
        constraints_list.append(jfc)

        # Limit joint
        model.positionLimitMargin = np.array([1.0] * model.nq)
        model.upperPositionLimit = np.array([0.5] * model.nq)
        model.upperPositionLimit[-1] = 2.0
        model.lowerPositionLimit = np.array([-0.5] * model.nq)
        jlc_raw = pin.JointLimitConstraintModel(model, [1, 2, 3])
        constraints_std_vec.append(pin.ConstraintModel(jlc_raw))
        constraints_list.append(jlc_raw)

        # Create some contacts and active limits
        data = model.createData()
        geom_data = geom_model.createData()
        q = np.array(
            [
                -0.49999217,
                -0.36846221,
                0.25560532,
                -0.72022615,
                0.66439729,
                0.13124927,
                -0.1504133,
                -0.45295538,
                1.19716179,
            ]
        )
        pin.forwardKinematics(model, data, q)
        pin.updateGeometryPlacements(model, data, geom_model, geom_data)
        pin.computeCollisions(model, data, geom_model, geom_data, q, False)
        pin.computeContactPatches(geom_model, geom_data)

        jlc = pin.ConstraintModel(jlc_raw).extract()
        jlc.resize(model, data, jlc.createData())

        frictional_points_list = []
        for col_pair, col_res, patch_res in zip(
            geom_model.collisionPairs,
            geom_data.collisionResults,
            geom_data.contactPatchResults,
        ):
            if col_res.isCollision() and patch_res.numContactPatches() > 0:
                geom_id1 = col_pair.first
                geom_id2 = col_pair.second
                geom1 = geom_model.geometryObjects[geom_id1]
                geom2 = geom_model.geometryObjects[geom_id2]
                joint_id1 = geom1.parentJoint
                joint_id2 = geom2.parentJoint
                oMi1 = data.oMi[joint_id1]
                oMi2 = data.oMi[joint_id2]
                patch = patch_res.getContactPatch(0)
                contact_normal = patch.getNormal()
                oMc = pin.SE3.Identity()
                contact_normal = contact_normal / np.linalg.norm(contact_normal)
                e_ref = (
                    np.array([0, 1, 0])
                    if np.isclose(contact_normal[0], 1.0)
                    else np.array([1, 0, 0])
                )
                comp = np.cross(e_ref, contact_normal)
                comp = comp / np.linalg.norm(comp)
                oMc.rotation[:, 2] = contact_normal
                oMc.rotation[:, 0] = comp
                oMc.rotation[:, 1] = np.cross(contact_normal, comp)
                for i in range(patch.size()):
                    oMc.translation = patch.getPointShape1(i)
                    i1Mc = oMi1.actInv(oMc)
                    oMc.translation = patch.getPointShape2(i)
                    i2Mc = oMi2.actInv(oMc)
                    fp = pin.FrictionalPointConstraintModel(
                        model, joint_id1, i1Mc, joint_id2, i2Mc
                    )
                    gfp = pin.ConstraintModel(fp)
                    constraints_std_vec.append(gfp)
                    constraints_list.append(fp)
                    frictional_points_list.append(fp)

        # Sotre some values
        self.data = data
        self.geom_data = geom_data
        self.constraints_std_vec = constraints_std_vec
        self.constraints_list = constraints_list
        self.bilats_list = bilats_list
        self.bilat = bilat1
        self.welds_list = welds_list
        self.weld = weld1
        self.frictional_points_list = frictional_points_list
        self.fp = frictional_points_list[0]
        self.jlc_raw = jlc_raw
        self.jlc = jlc
        self.jfc = jfc
        self.one_of_each = [self.bilat, self.weld, self.fp, self.jlc, self.jfc]

    def test_bilateral(self):
        # Coherence between all inits
        for b_1, b_2 in itertools.product(self.bilats_list, self.bilats_list):
            self.assertTrue(b_1 == b_2)

        # Check size
        for b in self.bilats_list:
            self.assertTrue(b.size() == b.activeSize() == 3)

    def test_weld(self):
        # Coherence between all inits
        for w_1, w_2 in itertools.product(self.welds_list, self.welds_list):
            self.assertTrue(w_1 == w_2)

        # Check size
        for w in self.welds_list:
            self.assertTrue(w.size() == w.activeSize() == 6)

    def test_frictional_point(self):
        for fp in self.frictional_points_list:
            self.assertTrue(fp.size() == fp.activeSize() == 3)

    def test_joint_frictional(self):
        self.assertTrue(self.jfc.activeSize() == self.jfc.size() <= self.model.nv)

    def test_joint_limit(self):
        self.assertTrue(
            self.jlc_raw.activeSize() <= self.jlc_raw.size() <= 2 * self.model.nq
        )
        self.assertTrue(self.jlc_raw.activeSize() == 0)
        self.assertTrue(self.jlc.activeSize() == (self.jlc.size() - 1))
        self.assertTrue(self.jlc.activeSize() <= self.jlc.size() <= 2 * self.model.nq)

    def test_generic_methods(self):
        ref_set = set(range(self.model.nv))

        for gcm, ccm in zip(self.constraints_std_vec, self.constraints_list):
            # Test hierarchy
            self.assertTrue(gcm.extract() == ccm)
            self.assertTrue(gcm.shortname() == ccm.classname() == ccm.shortname())
            for i in range(gcm.activeSize()):
                self.assertTrue(
                    np.all(
                        np.where(gcm.getRowActiveSparsityPattern(i))[0]
                        == np.array(gcm.getRowActiveIndexes(i))
                    )
                )
                self.assertTrue(
                    np.all(
                        np.where(gcm.getRowActivableSparsityPattern(i))[0]
                        == np.array(gcm.getRowActivableIndexes(i))
                    )
                )
                self.assertTrue(
                    set(gcm.getRowActiveIndexes(i))
                    <= set(gcm.getRowActivableIndexes(i))
                    <= ref_set
                )
                self.assertTrue(gcm.activeSize() <= gcm.size())
            dummy_compliance = 0.1 * np.ones(gcm.size())
            gcm.compliance = dummy_compliance
            self.assertTrue(np.all(dummy_compliance == gcm.compliance))
            self.assertTrue(len(ccm.getActiveCompliance()) == ccm.activeSize())
            self.assertTrue(ccm.activeSize() == ccm.set.size() == ccm.set.dim())
            if not hasattr(ccm, "baumgarte_corrector_parameters"):
                self.assertTrue(isinstance(ccm, pin.FrictionalJointConstraintModel))
                do_except = False
                self.assertTrue("baumgarte_corrector_parameters" in dir(gcm))
                try:
                    gcm.baumgarte_corrector_parameters
                except Exception as _:
                    do_except = True
                self.assertTrue(do_except)
            else:
                ccm.baumgarte_corrector_parameters.Kp = 1.0
                ccm.baumgarte_corrector_parameters.Kd = 1.0

    def test_jacobians_methods(self):
        model = self.model
        data = self.data

        v = np.stack([np.ones(model.nv), 2 * np.ones(model.nv)], axis=1)
        for cmodel in self.one_of_each:
            gcmodel = pin.ConstraintModel(cmodel)
            cdata = cmodel.createData()
            gcdata = gcmodel.createData()
            g2cdata = pin.ConstraintData(cdata)

            self.assertTrue(g2cdata == gcdata)
            self.assertTrue(gcdata.extract() == cdata)
            self.assertTrue(g2cdata.extract() == cdata)

            self.assertTrue(
                gcdata.shortname() == cdata.classname() == cdata.shortname()
            )

            lamb = np.stack(
                [np.ones(cmodel.activeSize()), 2 * np.ones(cmodel.activeSize())], axis=1
            )

            cmodel.calc(model, data, cdata)
            jac = cmodel.jacobian(model, data, cdata)
            sig = cmodel.jacobianMatrixProduct(model, data, cdata, v)
            tau = cmodel.jacobianTransposeMatrixProduct(model, data, cdata, lamb)

            gcmodel.calc(model, data, gcdata)
            gjac = gcmodel.jacobian(model, data, gcdata)
            gsig = gcmodel.jacobianMatrixProduct(model, data, gcdata, v)
            gtau = gcmodel.jacobianTransposeMatrixProduct(model, data, gcdata, lamb)

            self.assertTrue(np.all(jac == gjac))
            self.assertTrue(np.all(sig == gsig))
            self.assertTrue(np.all(tau == gtau))
            self.assertTrue(np.all(jac @ v == sig))
            self.assertTrue(np.all(jac.T @ lamb == tau))

    def test_sets_of_constraints(self):
        # Test projection
        # bilat
        force = np.array([1] * 3)
        p_force = self.bilat.set.project(force)
        self.assertTrue(self.bilat.set.isInside(force))
        self.assertTrue(np.all(p_force == force))

        # weld
        force = np.array([1] * 6)
        p_force = self.weld.set.project(force)
        self.assertTrue(self.weld.set.isInside(force))
        self.assertTrue(np.all(p_force == force))

        # jlc
        force_out = np.array([1] * self.jlc.activeSize())
        p_force_out = self.jlc.set.project(force_out)
        self.assertTrue(not self.jlc.set.isInside(force_out))
        p2_force_out = self.jlc.set.project(p_force_out)
        self.assertTrue(self.jlc.set.isInside(p_force_out))
        self.assertTrue(np.all(p_force_out == p2_force_out))
        force_in = np.array([0] * self.jlc.activeSize())
        p_force_in = self.jlc.set.project(force_in)
        self.assertTrue(self.jlc.set.isInside(p_force_in))
        self.assertTrue(np.all(p_force_in == force_in))

        # jfc
        force_out = np.array([2] * self.jfc.activeSize())
        p_force_out = self.jfc.set.project(force_out)
        self.assertTrue(not self.jfc.set.isInside(force_out))
        p2_force_out = self.jfc.set.project(p_force_out)
        self.assertTrue(self.jfc.set.isInside(p_force_out))
        self.assertTrue(np.all(p2_force_out == p_force_out))
        force_in = np.array([0] * self.jfc.activeSize())
        p_force_in = self.jfc.set.project(force_in)
        self.assertTrue(self.jfc.set.isInside(p_force_in))
        self.assertTrue(np.all(p_force_in == force_in))

        # fp
        self.fp.set.mu = 2.0
        force_out = np.array([0.0, 3.0, 1.0])
        p_force_out = self.fp.set.project(force_out)
        self.assertTrue(not self.fp.set.isInside(force_out))
        p2_force_out = self.fp.set.project(p_force_out)
        self.assertTrue(self.fp.set.isInside(p_force_out))
        self.assertTrue(np.all(p2_force_out == p_force_out))

        force_in = np.array([0.0, 1.0, 1.0])
        p_force_in = self.fp.set.project(force_in)
        self.assertTrue(self.fp.set.isInside(force_in))
        self.assertTrue(np.all(p_force_in == force_in))


if __name__ == "__main__":
    unittest.main()
