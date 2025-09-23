#  QSARpy, a knowledge discovery toolkit for automatic QSAR rules induction.
#  Copyright (C) 2014-2021  Thomas Ferrari

#    This file is part of QSARpy.
#
#    QSARpy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation version 2 of the License.
#
#    QSARpy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with QSARpy.  If not, see <http://www.gnu.org/licenses/>.

import itertools
import statistics

import networks as net
import QSAR_networks
from datasets import Testset, QSARTrainingset
from parameters import DECIMAL

FORMAT = '{:.{d}f}'
SGN_FORMAT = '{:+{d}f}'


class QSARModel(dict):

    def __init__(self, trainingset: QSARTrainingset, **kwargs):
        """ kwargs:  min_samples, max_atoms, max_error, deep_learning """
        super().__init__(kwargs)
        # self['train_path'] = str(trainingset.path)
        self['train_name'] = trainingset.path
        self['pop_outliers'] = trainingset.excludeOutliers
        self['speed_up'] = trainingset.fast_fragmentation
        self['depth'] = trainingset.fragmentation_level
        self['type'] = 'QSAR'
        self.__dict__ = self.copy()  # SET keys AS ATTRIBUTES
        self.trainingset = trainingset
        self.ruleset = []  # None
        self.neuron_log = None  # to export database

        # MONKEY PATCH to share networks with SAR
        net.Network = QSAR_networks.QSARNetwork
        net.LearningNetwork = QSAR_networks.QSARLearningNetwork

    def __repr__(self):
        return "{} structures + {} modulators".format(len(self.trainingset), len(self.ruleset))

    def getParameters(self):
        return dict(self)

    def summary(self, verbose=False):
        if not verbose or not self.ruleset:
            return str(self)
        text = "\n  Modulators:\n\n"
        text += "  ID\tabs err\tshift\tSMILES\n\n"
        for mod in self.ruleset:
            text += f"  {mod.getID()}\t{mod.error:.{DECIMAL}f} \t{mod.bias:+.{DECIMAL}f}\t{mod}\n"
        return text.expandtabs(8) + '\n  MODEL: ' + str(self)

    def extract(self):
        self.neuron_log = {}
        modulators = {}
        for level in range(1, self.depth + 1):
            modulators.update(self._getNeurons(level, modulators, net.learningNetworks))
            if self.deep_learning and level > 1:
                modulators.update(self._getNeurons(level, modulators, net.deepLearningNetworks))
        modulators = self.getPruned(modulators)
        self.ruleset = sorted(modulators.values(), key=lambda mod: (mod.sample_size, mod.cansmiles), reverse=True)
        for i, mod in enumerate(self.ruleset):
            mod.index = i + 1

    def _getNeurons(self, level, modulators, network_gen):
        trained = set()
        for structure in self.trainingset.structures():
            # propagating max_atoms to the generators saves time, but it also excludes some obs from self.neuron_log...
            # AT PRESENT: the check is implemented for deep networks only (NO DATABASES!)
            for network in network_gen(structure, self.trainingset, modulators, level, max_atoms=self.max_atoms):
                if network.error <= self.max_error:  # NETWORK error
                    # keep neurons UNIQUE (check and swap from neuron_log)
                    network.learner = self.neuron_log.setdefault(network.learner.cansmiles, network.learner)
                    network.trainLearner()
                    trained.add(network.learner)
        # assert trained.isdisjoint(self.modulators.values())
        return self._selectModulators(trained)

    def _selectModulators(self, neurons):
        # NOTE: approved modulators are not trained further (excluded in learningNetworks())
        modulators = {}
        for neuron in neurons:
            if neuron.unique_obs() >= self.min_samples and neuron.fragment.atoms <= self.max_atoms:
                neuron.activate()
                if neuron.error <= self.max_error:  # NEURON error
                    modulators[neuron.cansmiles] = neuron
        return modulators

    def getPruned(self, neuron_dict):  # redundant modulators: keep biggest
        modulators = {}
        for key, group in itertools.groupby(
                sorted(neuron_dict.values(), key=lambda x: self.trainingset.getBySub(x.cansmiles)),
                key=lambda x: self.trainingset.getBySub(x.cansmiles)):
            best_neuron = max(group, key=lambda x: x.fragment.size())
            modulators[best_neuron.cansmiles] = best_neuron
        return modulators

    def validateOnTrain(self):
        trainingset = self.trainingset
        if trainingset.excludeOutliers:
            testset = Testset(trainingset.copy())  # copy because I re-fragment
            if count := len(testset.predictOutliers(self.trainingset)):
                print("\n  * Detected %s structures SIMILAR to OUTLIERS" % count)
            testset.grind(self.depth, self.speed_up)  # fragmenting outliers left unfragmented
        else:
            testset = Testset(self.trainingset)
        self.predict(testset, testset.fragmentation_level)
        self.printStatistics(testset)
        trainingset.predictions = testset.predictions  # to output training predictions
        assert testset.path == trainingset.path  # DEBUG

    def predictFurther(self, dataset, depth):
        return self.predict(dataset, depth)

    def predictSameFragmets(self, testset):
        for test_structure in list(testset.structures()):
            matches = {}
            for train_structure in self.trainingset.structures():
                if train_structure.cansmiles == test_structure.cansmiles:  # If validating on train
                    continue
                if test_structure.checkSameFragments(train_structure, self.depth):
                    matches[train_structure] = train_structure.target
            if len(matches) > 0:
                prediction = statistics.median(matches.values())
                # prediction = sum(matches.values()) / len(matches)
                MAE = sum([abs(s - prediction) for s in matches.values()]) / len(matches)
                if MAE <= self.max_error:
                    testset._setPrediction(test_structure.cansmiles, "Same fragments", matches.keys(), prediction)

    def predict(self, dataset: Testset, depth=None):
        trainingset = self.trainingset
        # unpredicted = dataset.unpredicted
        # modulators = self.modulators
        modulators = {mod.cansmiles: mod for mod in self.ruleset}
        if depth is None:
            depth = self.depth
        level = 1
        onTrain = dataset.path == self.trainingset.path

        self.predictSameFragmets(dataset)

        while dataset.unpredicted and level <= depth:  # x-networks can predict deeper

            network_dict = {}
            # Direct Networks
            # for structure in unpredicted.values():
            for structure in dataset.structures():
                network_dict.setdefault(structure.cansmiles, [])\
                    .extend(net.networkGen(structure, trainingset, modulators, level))  # , dataset.fragments))
            # Reverse Networks
            if depth <= self.depth:  # disable reverse prediction if depth > trainingset depth (in predict_further)
                for train_structure in trainingset.structures():
                    for network in net.networkGen(train_structure, dataset.unpredicted, modulators, level, reverse=True):  # , trainingset.fragments, reverse=True):
                        network_dict.setdefault(network.structure.cansmiles, []).append(network)

            # Prediction
            for cansmiles, networks in network_dict.items():
                if networks:
                    best_network = min(networks)
                    if best_network.error <= self.max_error:
                        dataset.setPrediction(cansmiles, best_network)

            # IF UNPREDICTED: X-NETWORK PREDICTION
            if level > 1:
                # for structure in list(unpredicted.values()):
                for structure in list(dataset.structures()):
                    networks = []
                    for te_depth, tr_depth in net.weak_compositions(level):
                        if tr_depth <= self.depth:  # disable if depth > trainingset depth (in predict_further)
                            networks.extend(net.deepNetworks(structure, trainingset, modulators, (te_depth, tr_depth), onTrain))
                            # break first avail prediction: do go deeper de-indent and remove break
                            if networks:
                                best_network = min(networks)
                                if best_network.error <= self.max_error:
                                    dataset.setPrediction(structure.cansmiles, best_network)
                                    break
            level += 1
        return f"Predicted {len(dataset.predictions)} molecular structures"

    def printStatistics(self, dataset, verbose=False):
        """ Can raise ValueError """
        dataset.validate(numeric=True)  # can raise ValueError
        print("\n\n  PREDICTION STATISTICS:")

        def printMAE(description, errors, total=None):
            if errors:
                samples = len(errors)
                sum_abs_errors = sum((abs(error) for error in errors))
                if total:
                    err_desc = "  MAE (Mean Absolute Error) = "
                    perc = " ({0:.1f}%)".format(samples * 100 / total)
                else:
                    err_desc = "   MAE = "
                    perc = ""
                print(description + str(samples) + perc)
                print(err_desc, round(sum_abs_errors / samples, DECIMAL))

        error_dict = {}
        for record in dataset.predictions.values():
            error_dict.setdefault(record['mode'], []).append(record['error'])

        if verbose or error_dict.get('similar to outliers'):
            printMAE("\n  Known structures (identities): ", error_dict.get('identity'))
            printMAE("\n  Predicted outliers: ", error_dict.get('similar to outliers'))
            if verbose:
                printMAE("\n  Predicted by same fragments: ", error_dict.get('Same fragments'))
                printMAE("\n  Predicted by direct decomposition: ", error_dict.get('direct'))
                printMAE("\n  Predicted by reverse decomposition: ", error_dict.get('reverse'))
                printMAE("\n  Predicted by substitution: ", error_dict.get('generic'))
            else:
                printMAE("\n  Predicted by decomposition: ", error_dict.get('direct', [])
                         + error_dict.get('reverse', []) + error_dict.get('generic', []))

        all_errors = [error for sublist in error_dict.values() for error in sublist]
        printMAE("\n  TOTAL predicted: ", all_errors, len(dataset))

        if R2 := dataset.r2():
            print("  R2 = {}".format(R2))
            return R2
