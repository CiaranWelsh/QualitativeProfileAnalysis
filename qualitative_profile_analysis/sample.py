import os, glob
import numpy as np
import scipy.stats as stats
import pandas


class _Sampler:

    def __init__(self, sbml, config):
        self.sbml = sbml
        self.config = config

    def sample(self):
        raise NotImplementedError

    def _sample1dist(self):
        return NotImplementedError


class MonteCarloSampler(_Sampler):

    def sample(self):
        if self.config['settings']['sampler'] != 'monte_carlo':
            raise ValueError('Must set "settings.sampler = monte_carlo" to use '
                             'monte carlo sampling')
        species = self.config['species']

        # check for special case when all distributions and parameters to the distributions are equal
        df = pandas.DataFrame(species).transpose()
        if len(set(df['distribution'] == 1)) and len(set(df['loc'])) == 1 and len(set(df['scale'])) == 1:
            return list(self._sample1dist())[0]

        sample_list = np.array()
        for k, v in species.items():
            dist = v['distribution']
            v.pop('distribution')
            np.append(dist(**v, size=1), sample_list)
        return list(sample_list)[0]

    def _sample1dist(self):
        """
        Sample the required number of time from a single distribution.

        For the special case when only one distribution with one set of parameters is being used
        it is more efficient to sample using the inbuilt 'size' argument.

        Returns:

        """
        species = self.config['species']
        n = len(species)
        keys = list(species.keys())
        dist = species[keys[0]]['distribution']
        loc = species[keys[0]]['loc']
        scale = species[keys[0]]['scale']

        samples = dist(loc=loc, scale=scale).rvs(n)
        yield samples
