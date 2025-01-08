import unittest
from pathlib import Path

import pinocchio as pin

import importlib.util

mujoco_spec = importlib.util.find_spec("mujoco")
mujoco_found = mujoco_spec is not None


class TestMJCFBindings(unittest.TestCase):
    def test_load(self):
        model = pin.Model()
        current_dir = Path(__file__).parent
        model_dir = current_dir / "../models/"
        model_path = model_dir / "closed_chain.xml"
        model = pin.buildModelFromMJCF(model_path, model)
        constraint_models = pin.buildConstraintModelsFromMJCF(model_path, model)
        self.assertEqual(len(constraint_models), 4)


@unittest.skipUnless(mujoco_found, "Needs MuJoCo.")
class TestMJCFBindingsWithMujoco(unittest.TestCase):
    def test_cassie(self):
        import mujoco
        import numpy as np

        model_pin = pin.Model()
        current_dir = Path(__file__).parent
        model_dir = current_dir / "../models/"
        model_path = model_dir / "closed_chain.xml"
        model_pin = pin.buildModelFromMJCF(model_path, model_pin)
        constraint_models_pin = pin.buildConstraintModelsFromMJCF(model_path, model_pin)
        data_pin = model_pin.createData()
        q0_pin = model_pin.referenceConfigurations["home"]
        pin.forwardKinematics(model_pin, data_pin, q0_pin)
        model_mj = mujoco.MjModel.from_xml_path(str(model_path))
        data_mj = mujoco.MjData(model_mj)
        q0_mujoco = model_mj.key_qpos[0].copy()
        data_mj.qpos = q0_mujoco
        mujoco.mj_fwdPosition(model_mj, data_mj)
        for joint_id in range(model_mj.njnt):  # Total number of joints
            assert (
                np.linalg.norm(
                    data_mj.xanchor[joint_id] - data_pin.oMi[joint_id + 1].translation
                )
                < 1e-6
            )
        return


if __name__ == "__main__":
    unittest.main()
