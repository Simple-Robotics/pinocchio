import unittest
import pinocchio as pin
from test_case import PinocchioTestCase as TestCase


class TestBilateralBindings(TestCase):
    def test_vector(self):
        m = pin.Model()
        placement1 = pin.SE3.Identity()
        list_bpcm = pin.StdVec_BilateralPointConstraintModel()
        bpcm = pin.BilateralPointConstraintModel(m, 0, placement1)
        list_bpcm.append(bpcm)
        assert list_bpcm[-1].joint1_id == 0


if __name__ == "__main__":
    unittest.main()
