import os, glob
import tellurium as te
import numpy as np
import pandas as pd
import re
from .sample import MonteCarloSampler

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Qt5Agg')


class Result:

    def __init__(self, data, labels, meta=None):
        self.data = data
        self.labels = labels
        self.meta = meta


class QualitativeProfileAnalysis:

    def __init__(self, sbml, config):
        self.sbml = sbml
        self.config = config
        self.model = te.loadSBMLModel(sbml)

        self.nspecies = len(self.config['species'])
        self.n_iterations = self.config['settings']['n_iterations']
        self.start = self.config['settings']['integration_options']['start_time']
        self.end = self.config['settings']['integration_options']['end_time']
        self.intervals = self.config['settings']['integration_options']['intervals']

        self.ics = np.zeros(self.n_iterations * self.nspecies).reshape(self.n_iterations, self.nspecies)
        width = self.nspecies + 1 + 1  # +1 for time +1 for iteration id
        self.results = np.zeros(width * self.intervals * self.n_iterations) \
            .reshape(self.intervals * self.n_iterations, width)

        self.current_iter = 0
        self.sampler = MonteCarloSampler(self.sbml, self.config)

    def result_labels(self):
        ic_names = self.model.getFloatingSpeciesConcentrationIds()
        ic_names = [re.findall('\[(.*)\]', i)[0] for i in ic_names]
        return ['simulation_id', 'time'] + ic_names

    def run1(self):
        self.model.reset()
        ic_names = self.model.getFloatingSpeciesConcentrationIds()
        ics_vals = self.sampler.sample()
        for i in range(len(ic_names)):
            specie = re.findall('\[(.*)\]', ic_names[i])
            if len(specie) != 1:
                raise ValueError
            specie = specie[0]
            setattr(self.model, specie, ics_vals[i])
        # print(dict(zip(ic_names, ics_vals)))
        try:
            return self.model.simulate(self.start, self.end, self.intervals)
        except Exception:
            return None

    def run(self):
        for i in range(self.n_iterations):
            simulation_data = self.run1()
            simulation_id = np.array([i] * self.intervals)
            current_index = i * self.intervals
            self.results[current_index:current_index + self.intervals, 0] = simulation_id
            self.results[current_index:current_index + self.intervals, 1:] = simulation_data
        return Result(self.results, self.result_labels(), meta=self.config)


class Plotter:

    def __init__(self, results):
        if not isinstance(results, Result):
            raise ValueError('Expected class of type "Result" but got {} '.format(type(results)))
        self.results = results

        self.long_form_df = self._long_form()
        self.short_form_df = self._short_form()

    def _long_form(self):
        df = pd.DataFrame(self.results.data, columns=self.results.labels)
        df = df.set_index(['simulation_id', 'time'])
        df = pd.DataFrame(df.stack())
        df.columns = ['value']
        df.index.names = ['simulation_id', 'time', 'species']
        df = df.reset_index()
        return df

    def _short_form(self):
        data = self.results.data
        df = pd.DataFrame(data)
        df.columns = self.results.labels
        df = df.set_index(['simulation_id', 'time'])
        return df

    def plot_all_one_canvas(self, fname=None, linewidth=1, **kwargs):
        sns.set_context(context='talk')

        print(self.long_form_df)

        fig = plt.figure()
        sns.lineplot(x='time', y='value', hue='species', data=self.long_form_df,
                     units='simulation_id', estimator=None,
                     linewidth=linewidth, **kwargs
                     )
        sns.despine(fig=fig, top=True, right=True)
        plt.legend(loc=(1, 0.1))

        if fname is None:
            plt.show()
        else:
            plt.savefig(fname, dpi=300, bbox_per_inches='tight')

    def plot_variable(self, variable, fname=None, **kwargs):
        sns.set_context(context='talk')
        print(self.short_form_df.head())
        fig = plt.figure()
        sns.lineplot(x='time', y=variable, data=self.short_form_df.reset_index(), **kwargs,
                     units='simulation_id', estimator=None,
                     linewidth=1)
        sns.despine(fig=fig, top=True, right=True)

        if fname is None:
            plt.show()
        else:
            plt.savefig(fname, dpi=300, bbox_per_inches='tight')

    def plot_as_matrix(self, vars, ncols=3, fname=None, wspace=0.5,
                       hspace=0.5, figsize=(12, 12), **kwargs):
        sns.set_context(context='talk')
        vars = [i for i in vars if i not in ['simulation_id', 'time']]
        df = self.short_form_df[vars]
        n_to_plot = len(vars)
        nrows = int(np.floor(n_to_plot / ncols))
        remainder = n_to_plot % ncols  # number in last column
        if remainder != 0:
            nrows = nrows + 1
        vars_iter = iter(vars)
        fig, ax = plt.subplots(ncols=ncols, nrows=nrows, sharex=True, figsize=figsize)
        print('len vars', len(vars))
        for i in range(nrows):
            for j in range(ncols):
                try:
                    current_var = vars[nrows*i+j]
                    sns.lineplot(
                        ax=ax[i, j], x='time', y=current_var,
                        data=df.reset_index(),
                        units='simulation_id', estimator=None,
                        linewidth=1,
                        **kwargs)
                except IndexError:
                    break

        sns.despine(fig=fig, top=True, right=True)
        plt.subplots_adjust(wspace=wspace, hspace=hspace)

        if fname is None:
            plt.show()
        else:
            plt.savefig(fname, dpi=300, bbox_per_inches='tight')

    # build a way to subsetthe simulations
