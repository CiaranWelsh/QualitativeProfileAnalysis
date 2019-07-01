import os, glob
import yaml
import tellurium as te
from collections import ChainMap
import re
import scipy.stats as stats

class Config(ChainMap):

    def __init__(self, sbml, settings_kwargs={}, species_kwargs={}, conditions_kwargs={}):
        self.model = te.loadSBMLModel(sbml)
        self._config = ChainMap()
        self.settings_kwargs = settings_kwargs
        self.species_kwargs = species_kwargs
        self.conditions_kwargs = conditions_kwargs

        self.config = self._create_config()
        super().__init__(self.config)

    def _settings_defaults(self):
        return dict(
            sampler='monte_carlo',
            n_iterations=10,
            integration_options=dict(
                start_time=0,
                end_time=100,
                intervals=101
            ),
            quantity_type='concentration',           # or particle numbers
        )

    def _settings(self, **kwargs):
        settings = self._settings_defaults()
        for k, v in kwargs.items():
            settings[k] = v
        return settings

    def _conditions_defaults(self):
        return dict()

    def _conditions(self, **kwargs):
        conditions = self._conditions_defaults()
        for k, v in kwargs.items():
            conditions[k] = v
        return conditions

    def _species_defaults(self):
        dct = {}
        for i in self.model.getFloatingSpeciesConcentrationIds():
            species = re.findall('\[(.*)\]', i)
            assert len(species) == 1
            species = species[0]
            dct[species] = dict(distribution=stats.uniform,
                                loc=0.1, scale=99.9)
        return dct

    def _species(self, **kwargs):
        species = self._species_defaults()
        for k, v in kwargs.items():
            species[k] = v
        return species

    def _create_config(self):
        settings = self._settings(**self.settings_kwargs)
        conditions = self._conditions(**self.conditions_kwargs)
        species = self._species(**self.species_kwargs)
        return dict(
            settings=settings,
            conditions=conditions,
            species=species
        )

    def to_yaml(self, filename):
        with open(filename, 'w') as f:
            yaml_string = yaml.dump(self.config, f, default_flow_style=False)
        return yaml_string

    @staticmethod
    def from_yaml(yaml_file):
        pass
