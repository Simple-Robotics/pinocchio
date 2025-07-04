import unittest

import numpy as np
import pinocchio as pin
from test_case import PinocchioTestCase as TestCase


class TestJointsAlgo(TestCase):
    def setUp(self):
        self.model = pin.buildSampleModelHumanoidRandom()

        qmax = np.full((self.model.nq, 1), np.pi)
        self.q = pin.randomConfiguration(self.model, -qmax, qmax)
        self.v = np.random.rand(self.model.nv)

    def test_basic(self):
        model = self.model

        q_ones = np.ones(model.nq)
        self.assertFalse(pin.isNormalized(model, q_ones))
        self.assertFalse(pin.isNormalized(model, q_ones, 1e-8))
        self.assertTrue(pin.isNormalized(model, q_ones, 1e2))

        q_rand = np.random.rand(model.nq)
        q_rand = pin.normalize(model, q_rand)
        self.assertTrue(pin.isNormalized(model, q_rand))
        self.assertTrue(pin.isNormalized(model, q_rand, 1e-8))

        self.assertTrue(abs(np.linalg.norm(q_rand[3:7]) - 1.0) <= 1e-8)

        q_next = pin.integrate(model, self.q, np.zeros(model.nv))
        self.assertApprox(q_next, self.q)

        v_diff = pin.difference(model, self.q, q_next)
        self.assertApprox(v_diff, np.zeros(model.nv))

        q_next = pin.integrate(model, self.q, self.v)
        q_int = pin.interpolate(model, self.q, q_next, 0.5)

        self.assertApprox(q_int, q_int)

        value = pin.squaredDistance(model, self.q, self.q)
        self.assertTrue((value <= 1e-8).all())

        dist = pin.distance(model, self.q, self.q)
        self.assertTrue(dist <= 1e-8)

        q_neutral = pin.neutral(model)
        self.assertApprox(q_neutral, q_neutral)

        q_rand1 = pin.randomConfiguration(model)
        q_rand2 = pin.randomConfiguration(model, -np.ones(model.nq), np.ones(model.nq))

        self.assertTrue(pin.isSameConfiguration(model, self.q, self.q, 1e-8))

        self.assertFalse(pin.isSameConfiguration(model, q_rand1, q_rand2, 1e-8))

        lgo1 = pin.lieGroup(model)
        self.assertTrue(model.nq == lgo1.nq)
        self.assertTrue(model.nv == lgo1.nv)

        lgo2 = pin.LieGroup()
        for j in model.joints[1:]:
            lgo2 *= j.lieGroup()

        self.assertTrue(lgo1.name == lgo2.name)
        self.assertTrue(lgo1 == lgo2)
        self.assertApprox(lgo1.neutral, pin.neutral(model))

    def test_derivatives(self):
        model = self.model

        q = self.q
        v = self.v

        J0, J1 = pin.dIntegrate(model, q, v)
        res_0 = pin.dIntegrate(model, q, v, pin.ARG0)
        res_1 = pin.dIntegrate(model, q, v, pin.ARG1)

        self.assertApprox(J0, res_0)
        self.assertApprox(J1, res_1)

        q_next = pin.integrate(model, q, v)

        J0, J1 = pin.dDifference(model, q, q_next)
        res_0 = pin.dDifference(model, q, q_next, pin.ARG0)
        res_1 = pin.dDifference(model, q, q_next, pin.ARG1)

        self.assertApprox(J0, res_0)
        self.assertApprox(J1, res_1)

        TM = pin.tangentMap(model, q)

        TMv1 = TM.reshape(model.nq, model.nv) @ v.reshape(model.nv, 1)
        TMv2 = pin.tangentMapProduct(model, q, v.reshape(model.nv, 1))
        self.assertApprox(TMv1, TMv2)

        TMq1 = TM.reshape(model.nq, model.nv).T @ q.reshape(model.nq, 1)
        TMq2 = pin.tangentMapTransposeProduct(model, q, q.reshape(model.nq, 1))
        self.assertApprox(TMq1, TMq2)

        nvs, idx_vs = pin.indexvInfo(model)
        TMc = pin.compactTangentMap(model, q)
        TM_recons = np.zeros((model.nq, model.nv))
        for i in range(model.nq):
            TM_recons[i, idx_vs[i] : (idx_vs[i] + nvs[i])] = TMc[i, : nvs[i]]
        self.assertApprox(TM, TM_recons)


if __name__ == "__main__":
    unittest.main()
