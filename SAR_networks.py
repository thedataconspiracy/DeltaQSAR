import statistics

from networks import *
from parameters import DECIMAL


class SARNeuron(Neuron):

    learning_modes = ('direct', 'reverse')

    def __init__(self, cansmiles):
        super().__init__(cansmiles)
        labels = [key for key in SARNetwork.binarizer if isinstance(key, str)]
        self.sample_matrix = {key1: {key2: [] for key2 in self.learning_modes} for key1 in labels}
        self.bias_matrix = {key1: {key2: None for key2 in self.learning_modes} for key1 in labels}
        # INIZIALIZZALO SENZA KEYS
        # self.sample_size e` implementato in activate() come somma di tutti i samples, serve solo per ordinare

    def activate(self, label, mode):
        samples = self.sample_matrix[label][mode]
        self.bias_matrix[label][mode] = round(statistics.mean(samples), DECIMAL)
        self.sample_size += len(samples)


class SARNetwork(Network):

    binarizer = None

    @classmethod
    def setBinaryLabels(cls, label_0, label_1):
        cls.binarizer = {label_0: 0, label_1: 1, 0: label_0, 1: label_1}

    def __init__(self, outcome, apply_to, minus_neurons=(), plus_neurons=()):
        Network.__init__(self, outcome, apply_to, minus_neurons, plus_neurons)
        dir_biases = [neuron.bias_matrix[self.base.target]['direct'] for neuron in self.neurons('PLUS')]
        rev_biases = [neuron.bias_matrix[self.base.target]['reverse'] for neuron in self.neurons('MINUS')]
        biases = rev_biases + dir_biases
        if None in biases:
            # print("None in biases!")
            raise UserWarning
        assert all(value <= 0 for value in biases) or all(value >= 0 for value in biases)  # DEBUG
        self.output_log = total_probability([abs(bias) for bias in biases])
        switch_probability = sum(self.output_log)
        apply_to = self.binarizer[self.base.target]
        if apply_to == 0:
            self.num_prediction = switch_probability
        elif apply_to == 1:
            self.num_prediction = 1 - switch_probability
        self.prediction = self.binarizer[round(self.num_prediction)]  # round: 0 or 1
        self.switch_probability = round(switch_probability, DECIMAL)
        self.error = round(min(self.switch_probability, 1 - self.switch_probability), DECIMAL)

    def __lt__(self, other):  # to get BEST NETWORK
        return (self.error, self.switch_probability, str(self)) < (other.error, other.switch_probability, str(other))


class SARLearningNetwork(SARNetwork, LearningNetwork):  # Multiple inheritance: it gets methods from LearningNetworks

    def __init__(self, outcome, apply_to, learner_smiles, common_smiles, minus_neurons=(), plus_neurons=()):
        SARNetwork.__init__(self, outcome, apply_to, minus_neurons, plus_neurons)
        self.learner = SARNeuron(learner_smiles)
        self.common_smiles = common_smiles

    # "+ LEARNER" MA IN SAR POTREBBE ANCHE ESSERE "- LEARNER"
    def __repr__(self):
        # return super().__repr__() + ' + LEARNER: {}'.format(self.learner) +\
        return super().__repr__() +\
               f"\nBASE = {self.base.target}" +\
               f"\n + MODS = {[m.bias_matrix[self.base.target]['direct'] for m in self.neuron_dict['PLUS']]}" +\
               f"\n - MODS = {[m.bias_matrix[self.base.target]['reverse'] for m in self.neuron_dict['MINUS']]}" +\
               f"\nOUTCOME = {self.structure.target}\n"

    def trainLearner(self):
        mode = self.type()
        label = self.base.target
        input = self.num_prediction
        output = self.binarizer[self.structure.target]
        assert 0 <= input <= 1 and (output == 0 or output == 1)  # DEBUG
        sample = round(output - input, DECIMAL)

        # Never move backward (non tenere in considerazione cambiamenti di stato del livello precedente)
        # e se invece il bias fosse sempre 0 o 1 e avesse un "error" che si propaga?
        if self.binarizer[label] == 0 and sample < 0 or self.binarizer[label] == 1 and sample > 0:
            sample = 0

        self.learner.sample_matrix[label][mode].append(sample)
        self.learner.learningNetworks.append(self)


def total_probability(probabilities, idx=0):
    if idx > len(probabilities) - 1:
        # return 0
        return []
    prob = probabilities[idx]
    coeff = 1
    for prev in probabilities[:idx]:
        coeff *= (1 - prev)
    idx += 1
    # return coeff * prob + total_probability(probabilities, idx)
    return [coeff * prob] + total_probability(probabilities, idx)
