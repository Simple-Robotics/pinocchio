import unittest
from pathlib import Path

import pinocchio as pin


class TestMJCFBindings(unittest.TestCase):
    def test_load(self):
        model = pin.Model()
        current_dir = Path(__file__).parent
        model_dir = current_dir / "../models/"
        model_path = model_dir / "closed_chain.xml"
        model = pin.buildModelFromMJCF(model_path, model)
        constraint_models = pin.buildConstraintModelsFromMJCF(model_path, model)
        self.assertEqual(len(constraint_models), 4)


if __name__ == "__main__":
    unittest.main()
