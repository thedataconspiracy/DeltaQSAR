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


from collections import Counter
from grinder import cansmiles_dict


#               #### INTERFACE ONLY (must be subclassed) ####

class Neuron:

    def __init__(self, cansmiles):
        self.cansmiles = cansmiles
        self.fragment = cansmiles_dict[cansmiles]  # should be a method?
        self.index = None
        self.learningNetworks = []  # to output DataBase AND unique_obs(), this takes a lot of space i guess...
        self.sample_size = 0

    def __repr__(self):
        return self.cansmiles

    def __lt__(self, other):  # sorted by selector (in cluster.py)
        return str(self) < str(other)

    def getID(self):
        return 'M' + str(self.index)


class Network:
    # Subclasses must compute self.error and self.prediction

    def __init__(self, outcome, apply_to, minus_neurons=(), plus_neurons=()):
        # Network shouldn't have 'base' and 'structure' but just (id, target, cansmiles)
        self.base = apply_to
        self.structure = outcome
        self.neuron_dict = {'MINUS': Counter(minus_neurons), 'PLUS': Counter(plus_neurons)}
        self.error = None
        self.prediction = None
        self.output_log = None  # Used for output csv files

    def __repr__(self):
        return f"{self.structure} = {self.base} - {self.neuron_dict['MINUS'] or None} + {self.neuron_dict['PLUS'] or None}"

    # WARNING: defining only __lt__ is BAD for future python versions, but if I override __eq__ is too slow...)
    def __lt__(self, other):
        return str(self) < str(other)

    def __len__(self):
        return sum(self.neuron_dict['MINUS'].values()) + sum(self.neuron_dict['PLUS'].values())

    def type(self):
        if not self.neuron_dict['MINUS']:
            return 'direct'
        elif not self.neuron_dict['PLUS']:
            return 'reverse'
        else:
            return 'generic'

    def neurons(self, mono_key=None):
        assert mono_key is None or mono_key in ('MINUS', 'PLUS')  # DEBUG
        for key in ('MINUS', 'PLUS'):
            if mono_key is not None and mono_key != key:
                continue
            for neuron in self.neuron_dict[key].elements():
                yield neuron


class LearningNetwork(Network):

    def __init__(self, outcome, base, learner_smiles, common_smiles, minus_neurons=(), plus_neurons=()):
        super().__init__(outcome, base, minus_neurons, plus_neurons)
        self.learner = Neuron(learner_smiles)
        self.common_smiles = common_smiles

    def __repr__(self):
        return super().__repr__() + ' + LEARNER: {}'.format(self.learner)

    def type(self):
        assert self.base.atoms != self.structure.atoms  # DEBUG
        if self.base.atoms > self.structure.atoms:
            return 'reverse'
        else:
            return 'direct'

    def trainLearner(self):
        pass


#     #############################################################################################


# LEARNING ###

def learningNetworks(structure, sub_dataset, modulators, level, max_atoms=None, reverse=False):  # non uso max_atoms (solo deep)
    log = []
    for dec in sorted(structure.decompositions[level]):  # sorting.... DEBUG! (?)
        # what if I prescreen dec: only 1 fragment can be bigger than max_atoms...
        for fragment in _candidates(dec, modulators, learn=True):
            base = sub_dataset.get(cansmiles_dict[fragment].cansmiles_unstarred)
            if not base:
                continue
            fragments = [*dec]
            fragments.remove(fragment)  # removes the base (1st occurrence)
            for smiles in fragments:
                if smiles not in modulators:  # then it is the learner AND DO NOT train already approved modulators
                    fragments.remove(smiles)
                    if reverse:
                        net = LearningNetwork(base, structure, smiles, fragment,
                                                   minus_neurons=(modulators[key] for key in fragments))
                    else:
                        net = LearningNetwork(structure, base, smiles, fragment,
                                                   plus_neurons=(modulators[key] for key in fragments))
                    if str(net) not in log:  # different fragments, once removed the *, can be the same base (even if _candidates returns a set)
                        log.append(str(net))  # I use __repr__ cause overriding __eq__ takes too much time
                        yield net
                    break  # learner found


def deepLearningNetworks(structure, trainingset, modulators, level, max_atoms=None):
    log = []
    for depth_1, depth_2 in weak_compositions(level):
        for dec in structure.decompositions[depth_1]:
            for fragment in _candidates(dec, modulators, learn=True):
                # if trainingset.fragments[fragment].atoms < structure.atoms / 2:
                #     continue  # SKIP if MINUS modulators are subtracting too much
                plus = [*dec]
                plus.remove(fragment)  # common fragment
                # LEARNER
                for learner_smiles in plus:  # there is only 1 unknown in PLUS: the learner
                    if learner_smiles not in modulators:
                        plus.remove(learner_smiles)
                        break
                # assert learner_smiles not in modulators
                if max_atoms is not None and cansmiles_dict[learner_smiles].atoms > max_atoms:
                    continue  # collateral effect: these obs are excluded from neuron_log

                # find base structures with fragment in common
                for other_structure in trainingset.getBySub(fragment):
                    if cansmiles_dict[fragment].atoms < other_structure.atoms / 2:
                        continue  # SKIP if MINUS modulators are subtracting too much
                    if other_structure.cansmiles == structure.cansmiles:
                        continue
                    for dec_2 in other_structure.decompositions[depth_2]:
                        if fragment not in dec_2:
                            continue
                        if fragment not in _candidates(dec_2, modulators):
                            continue
                        minus = [*dec_2]
                        minus.remove(fragment)  # common fragment
                        net = LearningNetwork(structure, other_structure, learner_smiles, fragment,
                                              (modulators[key] for key in minus), (modulators[key] for key in plus))
                        if str(net) not in log:
                            log.append(str(net))
                            yield net


# PREDICTION ###

def networkGen(structure, sub_dataset, modulators, level, reverse=False):
    for dec in structure.decompositions[level]:
        for fragment in _candidates(dec, modulators):
            base = sub_dataset.get(cansmiles_dict[fragment].cansmiles_unstarred)
            if not base:
                continue
            neurons = [*dec]
            neurons.remove(fragment)  # removes the base (1st occurrence)
            if reverse:
                yield Network(base, structure, minus_neurons=(modulators[key] for key in neurons))
            else:
                yield Network(structure, base, plus_neurons=(modulators[key] for key in neurons))


def deepNetworks(test_structure, trainingset, modulators, levels, onTrain=False):
    networks = []
    for test_dec in test_structure.decompositions[levels[0]]:
        for fragment in _candidates(test_dec, modulators):
            plus_fragments = [*test_dec]
            plus_fragments.remove(fragment)  # 1st occurrence only
            for train_structure in trainingset.getBySub(fragment):
                # fragment is a COMMON fragment (test_structure, train_structure)
                # if frag_dict[fragment].atoms < train_structure.atoms / 2:
                #     continue  # SKIP if MINUS modulators are subtracting too much
                if onTrain and train_structure.cansmiles == test_structure.cansmiles:
                    continue
                for train_dec in train_structure.decompositions[levels[1]]:
                    if fragment not in train_dec:
                        continue
                    minus_fragments = [*train_dec]
                    minus_fragments.remove(fragment)  # 1st occurrence only
                    try:  # here I could SPEED UP using _candidates for filtering before/instead_of the try
                        networks.append(Network(test_structure, train_structure,
                                                (modulators[key] for key in minus_fragments),
                                                 # if key not in plus_fragments),
                                                (modulators[key] for key in plus_fragments)))
                                                 # if key not in minus_fragments)))
                                                # PLUS/MINUS same fragment is ignored!
                    except KeyError:
                        continue
                    except UserWarning:
                        continue
    # if len(networks) > 1:
    #     raise
    return networks


# UTILITIES

def _candidates(factorization, modulators, learn=False):
    # OPTIMIZATION: A factorization can be predictive only if it contains a fragment identical to a training structure
    # (called "base"), and all the other fragments are known modulators. This function checks if these preconditions are
    # matched and it returns a subset of "base candidates". Empty iterable otherwise.
    # count := unknown fragments (NOT in modulators):
    # count = 0 : proceed normally (return all)
    # count = 1 : that fragment must be the base, return it
    fragment_set = set(factorization)
    unknowns = fragment_set.difference(modulators)
    count = sum((factorization.count(elem) for elem in unknowns))  # using set speeds up, but then I have to count...
    if learn:  # For learning networks there is 1 extra degree of freedom (AND if count=0 skip)
        if count == 2:
            return unknowns
        elif count == 1:
            return fragment_set - unknowns  # bases = fragments - learner
    else:
        if count == 1:
            return unknowns
        elif count == 0:
            return fragment_set
    return ()


# returns same-depth pairs of levels
def weak_compositions(depth, width=2, parent=[]):
    depth_levels = []
    if width > 1:
        for i in range(1, depth):
            for x in weak_compositions(i, width - 1, parent + [depth - i, ]):
                # if depth <= 3 or max(x) <= 3:
                depth_levels.append(x)
    else:
        # if depth <= 3 or max(parent + [depth, ]) <= 3:
        depth_levels.append(parent + [depth, ])
    return depth_levels
