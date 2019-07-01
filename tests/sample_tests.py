import unittest
import os, glob
import numpy as np
import scipy.stats as stats
import tellurium as te

import site

site.addsitedir(os.path.join(os.path.dirname(__file__), 'qualitative_profile_analysis'))

from .test_models import *
from qualitative_profile_analysis.sample import *
from qualitative_profile_analysis.config import *


class SamplerTests(unittest.TestCase):

    def setUp(self) -> None:
        self.sbml = te.antimonyToSBML(TEST_MODEL2)
        self.config = Config(self.sbml)

    def test_montecarlo_sampler(self):
        s = MonteCarloSampler(self.sbml, self.config)
        self.assertEqual(len(s._sample()), 6)


if __name__ == '__main__':
    unittest.main()
