#  DeltaQSAR, a knowledge discovery toolkit for automatic SAR/QSAR rules induction.
#  Copyright (C) 2021  Thomas Ferrari


from copy import copy
from operator import attrgetter
from cluster import HierarchicalClustering as Cluster

from QSAR_model import QSARModel
from datasets import Dataset
from parameters import DECIMAL


class ClusterModel(QSARModel):

    def __init__(self, model):
        super().__init__(model.trainingset, **dict(model))
        self.original_model = model
        self.ruleset = [copy(neuron) for neuron in self.original_model.ruleset]
        # self.ruleset = deepcopy(self.original_model.ruleset)
        for mod in self.ruleset:  # FRAGMENTING MODULATORS
            mod.decompositions = [{(mod.cansmiles,)}]
        mod_dataset = Dataset({mod.cansmiles: mod for mod in self.ruleset})
        mod_dataset.grind()
        self.clustered_view = None

    def __repr__(self):
        return "{} structures + {} modulator CLUSTERS".format(len(self.trainingset), self.clustered_view.countClusters())

    def summary(self, verbose=False):
        cv = self.clustered_view
        if not verbose or not self.ruleset:
            return str(cv)
        text = '\n\n CLUSTERED VIEW:\n\n'
        for key in cv.labels(sort=True):
            tab_len = 12
            part = ''
            tab = '\t'
            if len(key) + 1 < tab_len:
                tab = tab * 2
            base = ' {}'.format(key) + tab
            if len(base) >= tab_len * 2:
                tab_len = 2
            part += base.expandtabs(tab_len)
            for c, cluster in enumerate(cv[key]):
                if c == 0:
                    spacer = ''
                else:
                    spacer = ' ' * len(base.expandtabs(tab_len))
                bias_string = '{shift:+.{d}f}'.format(shift=cluster[0].bias, d=DECIMAL)
                smiles_string = repr(cluster)[1:-1].replace(',', ' ')
                part += spacer + bias_string + ' ' + smiles_string + '\n'
            if len(cv[key]) > 1:
                text = text.rstrip() + '\n\n' + part + '\n'
            else:
                text += part
        if cv['PLACEHOLDERS']:
            text += "\n DUMMY MODULATORS: {}  (shift = 0)\n".format(len(cv['PLACEHOLDERS']))
            text += ' ' + repr(cv['PLACEHOLDERS'])[1:-1].replace(',', ' ') + '\n'
        if cv['REDUNDANT']:
            text += "\n REDUNDANT MODULATORS: {}\n".format(len(cv['REDUNDANT']))
            text += "  (omitted)\n"
        return text + '\n\n  MODEL: {}'.format(self)  # .decode('utf-8')

    def clusterize(self, tolerance):
        self.tolerance = self['tolerance'] = tolerance  # the dict is for getParam
        # RESET
        for i in range(len(self.ruleset)):
            self.ruleset[i].original_bias = self.original_model.ruleset[i].bias
            self.ruleset[i].generalized_bias = self.original_model.ruleset[i].bias

        unprocessed = {mod.cansmiles: mod for mod in self.ruleset}
        self.clustered_view = cv = ClusteredView(tolerance)
        # FILTER REDUNDANCY: identifico i GHOST ma non aggiorno ancora il bias
        #   ATTENZIONE: che la learning impara usando anche i neuroni ghost prima che siano generalizzati,
        #   dunque impara usando dei numeri che poi non esistono piu' (x es. ModAnalisys)
        redundant_modulators = _popRedundant(unprocessed, tolerance)
        # PLACEHOLDERS
        cv.updateKey('PLACEHOLDERS', _popPlaceHolders(unprocessed, tolerance))
        # CLUSTERING
        cv.update(unprocessed)
        cv.clustering()
        # CONFIRMING REDUNDANT MODULATORS
        # check if potential redundant are still so after bias changes made by ph and clustering
        rejected = {}
        for key, mod in redundant_modulators.copy().items():
            generalized_bias = sum([agg.generalized_bias for agg in mod.aggregates])
            if abs(generalized_bias - mod.original_bias) < tolerance:
                mod.generalized_bias = generalized_bias
            else:
                rejected[key] = redundant_modulators.pop(key)
                # rejected[mod.cansmiles] = fragment
        cv.updateKey('REDUNDANT', redundant_modulators)
        # REJECTED REDUNDANT: check for PH and CLUSTER
        cv.updateKey('PLACEHOLDERS', _popPlaceHolders(rejected, tolerance))
        cv.updateCluster(rejected)
        # DEBUG: test unit
        # clustered.check()

        # recalculating individual error
        for mod in self.ruleset:
            mod.bias = mod.generalized_bias
            mod.calcError()
            del mod.generalized_bias, mod.original_bias

        # cluster error as weighted average
        # 1 cluster must have 1 shift and 1 error, here I recalculate the global error
        for group_label in cv.labels():
            for cluster in cv[group_label]:
                tot_samples = 0
                sum_errors = 0
                for mod in cluster:
                    sum_errors += mod.error * mod.sample_size
                    tot_samples += mod.sample_size
                weighted_error = round(sum_errors / tot_samples, DECIMAL)
                for mod in cluster:
                    # if abs(mod.error - weighted_error) > tolerance:  # DEBUG
                    #     print(mod.error - weighted_error)
                    mod.error = weighted_error


class ClusteredView(dict):

    def __init__(self, tolerance):
        self.tolerance = tolerance
        self['REDUNDANT'] = []
        self['PLACEHOLDERS'] = []

    # get cluster labels (NO REDUNDANT OR PLACEHOLDERS)
    def labels(self, sort=False):  # OKKIO puoi ordinare solo dopo clustering
        labels = self.keys() - set(('REDUNDANT', 'PLACEHOLDERS'))  # PYTHON BUG IN WIN!! if not set removes 'C', 'O', etc
        if sort:
            labels = sorted(labels, key=lambda key: (self.cardinality(key), key), reverse=True)
        return labels

    def update(self, modulators):
        for mod in modulators.values():
            key = mod.fragment.cansmiles_unstarred  # fragment_dict[mod.cansmiles].cansmiles_unstarred
            if key in self:
                self[key].append(mod)
            else:
                self[key] = [mod]

    def updateKey(self, key, modulators):
        self[key].extend(modulators.values())
        self[key].sort(key=attrgetter('samples'), reverse=True)

    def updateCluster(self, modulators):
        for mod in modulators.values():
            key = mod.fragment.cansmiles_unstarred  # fragment_dict[mod.cansmiles].cansmiles_unstarred
            if key in self:
                clusters = self[key]
                delta = self.tolerance
                best_cluster = None
                for cluster in clusters:
                    generalized_bias = cluster[0].generalized_bias
                    new_delta = abs(mod.original_bias - generalized_bias)
                    if new_delta < delta:
                        delta = new_delta
                        best_cluster = cluster
                if best_cluster:
                    mod.generalized_bias = best_cluster[0].generalized_bias
                    best_cluster.append(mod)
                else:
                    self[key].append([mod])
            else:
                self[key] = [[mod]]

    def clustering(self):
        for key in self.labels():
            group = self[key]

            if len(group) == 1:
                self[key] = [group]
                continue

            level = self.tolerance * 2
            distance = lambda mod1, mod2: abs(mod1.original_bias - mod2.original_bias)
            clusters = Cluster(group, distance, linkage='complete').getlevel(level)
            self[key] = clusters
            for cluster in clusters:
                cluster.sort(key=attrgetter('samples'), reverse=True)
                if len(cluster) == 1:
                    continue
                ## MEDIA PESATA
                ## since the new bias is not the middle point of the cluster it *could* happen
                ## that the delta between new and old bias is slightly bigger than bias_thres
                w_mean = sum([e.original_bias * e.sample_size for e in cluster]) / \
                                   sum([e.sample_size for e in cluster])
                ## PUNTO MEDIO
                # cluster_biases = [m.original_bias for m in cluster]
                # w_mean = (max(cluster_biases) + min(cluster_biases)) / 2
                for mod in cluster:
                    mod.generalized_bias = round(w_mean, DECIMAL)

    def countClusters(self):
        return sum(len(self[key]) for key in self.labels())

    def cardinality(self, key):
        return sum(mod.sample_size for cluster in self[key] for mod in cluster)

    # ## DEBUG
    # def check(self):    #
    #     def checkDistance(modulator):
    #         assert abs(modulator.generalized_bias - modulator.original_bias) < self.tolerance    #
    #     for key in self.labels():
    #         for cluster in self[key]:
    #             for mod in cluster:
    #                 assert mod.generalized_bias == cluster[0].generalized_bias
    #                 # if abs(mod.original_bias - mod.generalized_bias) > self.tolerance:
    #                 # print "ALARM! clustering"
    #                 # print mod, abs(mod.original_bias - mod.generalized_bias)
    #     for ph in self['PLACEHOLDERS']:
    #         assert ph.generalized_bias == 0
    #         checkDistance(ph)
    #     for combo in self['REDUNDANT']:
    #         assert combo.generalized_bias == sum([m.generalized_bias for m in combo.aggregates])
    #         checkDistance(combo)


def _popRedundant(modulators, bias_thres):
    basic_fragments = {}
    ghost_modulators = {}  # combination of basic modulators
    # Prima metto in basic i modulatori di 1 solo atomo, poi guardo quelli di 2 atomi: se non possono
    # essere generalizzati dai basic gia` esistenti li aggiungo... e  cosi` via
    for key, mod in sorted(modulators.items(), key=lambda item: item[1].fragment.atoms):
        mod.aggregates = None
        bestFactors = None
        delta_b = bias_thres
        basic_fragments[key] = mod  # dict of basic (non redundant) modulators (see mod.aggregates)
        for level in range(1, len(mod.decompositions)):
            for factors in mod.decompositions[level]:
                if set(basic_fragments.keys()).issuperset(factors):
                    generalized_bias = sum([basic_fragments[f].generalized_bias for f in factors])
                    new_delta = abs(generalized_bias - mod.original_bias)
                    # non metto <= senno` prende quella del livello piu` profondo
                    if new_delta < delta_b:
                        delta_b = new_delta
                        bestFactors = factors
            if bestFactors:
                mod.aggregates = [basic_fragments[f] for f in bestFactors]
                basic_fragments.pop(key)
                ghost_modulators[key] = modulators.pop(key)
                break
    return ghost_modulators


def _popPlaceHolders(modulators, bias_thres):
    placeholders = {}
    for key, mod in modulators.copy().items():
        if abs(mod.original_bias) < bias_thres:
            mod.generalized_bias = 0
            placeholders[key] = modulators.pop(key)
    return placeholders
