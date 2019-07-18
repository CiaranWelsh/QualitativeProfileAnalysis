import pandas
import unittest
import tellurium as te
import site
import os

site.addsitedir(os.path.join(os.path.dirname(__file__), 'qualitative_profile_analysis'))

from .test_models import *
from qualitative_profile_analysis.analysis import *
from qualitative_profile_analysis.config import *



class BaseClassTests(unittest.TestCase):

    def setUp(self) -> None:
        self.sbml = te.antimonyToSBML(SIMPLE_AKT_MODEL_SBML)
        self.config = Config(
            self.sbml,
            settings_kwargs=dict(
                n_iterations=100,
                integration_options=dict(
                    start_time=0,
                    end_time=100,
                    intervals=100
                ),
                loc=0.001, scale=10
            )
        )
        self.qpa = QualitativeProfileAnalysis(self.sbml, self.config)


    def test(self):
        print(self.qpa.results)
        p = Plotter(self.qpa.results)
        p.plot_all_one_canvas()
        plt.show()



























