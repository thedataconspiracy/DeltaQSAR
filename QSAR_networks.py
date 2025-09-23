import statistics

from networks import *
from parameters import DECIMAL


class QSARNeuron(Neuron):

    def __init__(self, cansmiles):
        super().__init__(cansmiles)
        self.samples = []
        self.bias = self.prediction = None
        self.error = None
        self.bases = []

    def activate(self):
        self.sample_size = len(self.samples)
        self.bias = round(statistics.mean(self.samples), DECIMAL)  # MSD
        self.prediction = self.bias  # compatibility SAR
        self.calcError()

    def calcError(self):
        # why don't you use standard deviation?
        MAE = statistics.mean([abs(s - self.bias) for s in self.samples])  # Mean Absolute Error
        self.error = round(MAE, DECIMAL)

    def modulate(self, value, reverse=False):
        if reverse:
            return value - self.bias
        return value + self.bias

    def unique_obs(self):
        hits = set((net.structure.cansmiles for net in self.learningNetworks))
        return len(hits)  # count only structure and not base (the one actually containing the fragment)


class QSARNetwork(Network):

    def __init__(self, outcome, apply_to, minus_neurons=(), plus_neurons=()):
        Network.__init__(self, outcome, apply_to, minus_neurons, plus_neurons)
        self.error = sum((neuron.error for neuron in self.neurons()))
        self.prediction = round(self._update(), DECIMAL)

    def __lt__(self, other):  # to get BEST NETWORK
        return (self.error, str(self)) < (other.error, str(other))

    # def confidence(self):
    #     mod_atoms = sum((neuron.fragment.atoms for neuron in self.neurons()))
    #     return 1 - mod_atoms / self.structure.atoms

    def _update(self):  # can be simplified
        output = self.base.target
        self.output_log = [output]
        _input = output
        for neuron in self.neuron_dict['MINUS'].elements():
            output = neuron.modulate(_input, reverse=True)
            self.output_log.append(output)
            _input = output
        for neuron in self.neuron_dict['PLUS'].elements():
            output = neuron.modulate(_input, reverse=False)
            self.output_log.append(output)
            _input = output
        return output


class QSARLearningNetwork(QSARNetwork, LearningNetwork):  # Multiple inheritance: it gets methods from LearningNetworks

    def __init__(self, outcome, base, learner_smiles, common_smiles, minus_neurons=(), plus_neurons=()):
        QSARNetwork.__init__(self, outcome, base, minus_neurons, plus_neurons)
        self.learner = QSARNeuron(learner_smiles)
        self.common_smiles = common_smiles

    def trainLearner(self):
        # input = self._update()
        input = self.prediction  # rounded!
        output = self.structure.target
        self.learner.samples.append(round(output - input, DECIMAL))
        self.learner.learningNetworks.append(self)


        # ERROR
        # tanimoto_coefficient = len(self.base.fp_bits & self.structure.fp_bits) \
        #                        / len(self.base.fp_bits | self.structure.fp_bits)
        # # Peso l'errore con la similitudine base/structure: se e` molto simile, l' errore e` dimezzato
        # return sum((neuron.error for neuron in self.neurons()))  * decimal.Decimal(1 - tanimoto_coefficient / 2)
        # if self.type() == 'direct':
        #     k = sum((n.fragment.atoms for n in self.neuron_dict['PLUS'])) / self.structure.atoms
        # elif self.type() == 'reverse':
        #     k = sum((n.fragment.atoms for n in self.neuron_dict['MINUS'])) / self.base.atoms
        # else:
        #     common = self.base.atoms - sum((n.fragment.atoms for n in self.neuron_dict['MINUS']))
        #     common_str = self.structure.atoms - sum((n.fragment.atoms for n in self.neuron_dict['PLUS']))
        #     incognita = sum((n.fragment.atoms for n in self.neuron_dict['PLUS'])) + sum((n.fragment.atoms for n in self.neuron_dict['MINUS']))
        #     k = incognita / (common + incognita)
        # # print(k, self.type(), self, self.structure)
        # return float(sum((neuron.error for neuron in self.neurons()))) * ((k/3) + 0.66)
