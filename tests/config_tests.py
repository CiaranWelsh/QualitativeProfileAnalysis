import unittest, os
import yaml
import tellurium as te
from collections import ChainMap

import site
site.addsitedir(os.path.join(os.path.dirname(__file__), 'qualitative_profile_analysis'))

from .test_models import *
from qualitative_profile_analysis.config import Config


class ConfigTests(unittest.TestCase):

    def setUp(self) -> None:
        self.sbml = te.antimonyToSBML(TEST_MODEL1)
        self.yaml_file = os.path.join(os.path.dirname(__file__), 'config.yaml')

    # def tearDown(self) -> None:
    #     if os.path.isfile(self.yaml_file):
    #         os.remove(self.yaml_file)

    def test_isa(self):
        config = Config(self.sbml)
        self.assertIsInstance(config, Config)
        self.assertIsInstance(config, ChainMap)

    def test_change_setting(self):
        config = Config(self.sbml, settings_kwargs={'sampling_mode': 'lhs'})
        self.assertEqual(config['settings']['sampling_mode'], 'lhs')

    def test_write_to_file(self):
        config = Config(self.sbml, settings_kwargs={'sampling_mode': 'lhs'})
        config.to_yaml(filename=self.yaml_file)
        self.assertTrue(os.path.isfile(self.yaml_file))





















if __name__ == '__main__':
    unittest.main()
