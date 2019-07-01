import pandas
import unittest
import tellurium as te
import site
import os
site.addsitedir(os.path.join(os.path.dirname(__file__), 'qualitative_profile_analysis'))

from .test_models import *
from qualitative_profile_analysis.analysis import *
from qualitative_profile_analysis.config import *


class AnalysisTests(unittest.TestCase):

    def setUp(self) -> None:
        self.sbml1 = te.antimonyToSBML(TEST_MODEL1)
        self.sbml2 = te.antimonyToSBML(TEST_MODEL2)

        self.sbml = [self.sbml1, self.sbml2]
        self.config = Config(
            self.sbml1,
            settings_kwargs=dict(
                n_iterations=10,
                integration_options=dict(
                    start_time=0,
                    end_time=100,
                    intervals=1000
                )
            ))

    def test_run(self):
        qpa = QualitativeProfileAnalysis(self.sbml1, self.config)
        results = qpa.run()
        df = pandas.DataFrame(results, columns=qpa.result_labels())
        df = df.set_index(['simulation_id', 'time'])
        self.assertEqual((110, 3), df.shape)

    def test_visualisation(self):
        qpa = QualitativeProfileAnalysis(self.sbml1, self.config)
        results = qpa.run()
        p = Plotter(results)
        p.plot_all_one_canvas()


    def test_visualisation2(self):
        qpa = QualitativeProfileAnalysis(self.sbml1, self.config)
        results = qpa.run()
        p = Plotter(results)
        p.plot_variable('C')




if __name__ == '__main__':
    unittest.main()


































