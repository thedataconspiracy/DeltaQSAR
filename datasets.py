#  DeltaQSAR, a knowledge discovery toolkit for automatic SAR/QSAR rules induction.
#  Copyright (C) 2021  Thomas Ferrari


import csv
import statistics
import string
from operator import attrgetter
from grinder import Grinder, cansmiles_dict
import chem

from parameters import DECIMAL


class Structure:

    def __init__(self, mol, smiles=None):
        super().__init__()
        self.cansmiles = chem.mol2smi(mol)
        self.atoms = chem.numAtoms(mol)
        self.fingerprint = chem.mol2fp(mol)  # ONLY for similarity and pre-screen SMARTS matching (SAR)
        self.mol = mol  # required for SAR (SMARTS matching) and for bond_map
        self.smiles = smiles or self.cansmiles  # if SDF and bond_map, I should make sure that the atoms are in canonical order
        self.target = None
        self.S_id = None
        self.decompositions = None  # (populated by the Grinder)

    def __repr__(self):
        return self.smiles

    # used only in split dataset
    def __lt__(self, other: 'Structure'):
        return (self._chemicalClass(), self.target, self.smiles) < (other._chemicalClass(), other.target, other.smiles)

    def getID(self):
        return self.S_id

    def getSimilar(self, dataset):
        """ Return the structures with the same fingerprint of self (ergo Tanimoto = 1)
            Structures whose atom number differs for more than 10% are excluded. """
        similar_dict = {}
        for other in dataset.structures():
            if self.cansmiles != other.cansmiles:
                atom_diff = abs(self.atoms - other.atoms)
                if atom_diff <= 0.1 * self.atoms and self.fingerprint == other.fingerprint:
                    similar_dict.setdefault(atom_diff, []).append(other)
        if similar_dict:  # return the group with lowest difference in atoms number
            return similar_dict.get(min(similar_dict.keys()))

    def checkSameFragments(self, other, max_depth):
        """ Check whether self and other have a decomposition in common (share the same fragments)
            at a given max_depth. """
        if self.atoms != other.atoms: return False
        for level in range(1, max_depth + 1):
            if self.decompositions[level] & other.decompositions[level]:
                return True

    def _chemicalClass(self):
        smiles_letters = set(self.cansmiles.upper()).intersection(string.ascii_letters)
        if smiles_letters.issubset('CNH'):
            return 'CN'
        if smiles_letters.issubset('COH'):
            return 'CO'
        if smiles_letters.issubset('CNOH'):
            return 'CNO'
        return 'other'


class Dataset(dict):
    """ A dictionary of chemical structures with their canonical SMILES as key. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.path = None
        self.fragmentation_level = 0
        self.fast_fragmentation = None
        self.fragments = {}  # {fragment: [structures]}

    def copy(self):
        new = type(self)(self)
        new.__dict__.update(self.__dict__)
        new.fragments = self.fragments.copy()
        for key,value in new.fragments.items():
            new.fragments[key] = [*value]
        return new

    def add(self, structure):
        """ Add a structure to the dataset. """
        # TO DO: give some feedback in case of duplicated structures !!!!!!!!!!!!!!!!!!!!!!!!!
        self[structure.cansmiles] = structure

    # def addHs(self):
    #     """ Add hydrogens. """
    #     for structure in self.structures():
    #         if structure.mol_addHs is not None:
    #             break
    #         structure.mol_addHs = chem.addedHs(structure.mol)

    def _try_make_numeric(self):
        """ Use for QSAR (numeric target): cast all structure.target to float. Can raise errors. """
        numeric_targets = [float(structure.target) for structure in self.values()]
        for structure, value in zip(self.values(), numeric_targets):
            structure.target = value

    def clear_fragmentation(self):
        """ Clear the Dataset and the contained Structures from fragmentation. """
        self.fragmentation_level = 0
        # self.max_atoms_fragmented = None
        self.fragments = {}
        for structure in self.values():
            structure.decompositions = None

    def getSplit(self, testFold=4, foldsNum=5):
        """ Split self into training set and test set (based on chemical class) and return them as new datasets. """
        trainingset = Dataset()
        testset = Dataset()
        structures = sorted(self.values())
        for c, structure in enumerate(structures):
            fold = (c % foldsNum)
            if fold == testFold:
                testset.add(structure)
            else:
                trainingset.add(structure)
        return trainingset, testset

    def percentile(self, attribute, percent):
        """ Return the percentile in the dataset of the given attribute of a structure. """
        getattribute = attrgetter(attribute)
        return statistics.quantiles((getattribute(structure) for structure in self.values()), n=100)[percent - 1]

    def get(self, cansmiles):
        """ Return the structure with the given cansmiles if present, else None. """
        return super().get(cansmiles)

    def structures(self):
        """ Return an object providing a view on the dataset structures. """
        return self.values()

    def getBySub(self, sub_cansmiles):
        """ Return an iterable of structures containing the given substructures. """
        return self.fragments.get(cansmiles_dict[sub_cansmiles]) or ()

    def grindProgress(self, max_depth=None, speed_up=False):
        """ Structural fragmentation with progress feedback. Populating dataset.fragments and structure.decompositions. """
        grinder = Grinder()
        for i, structure in enumerate(self.structures()):
            fragments = grinder.grind(structure, max_depth, speed_up)
            for fragment in fragments.values():
                # Binding substructure: superstructures
                self.fragments.setdefault(fragment, []).append(structure)
            yield i
        self.fragmentation_level = max_depth
        self.fast_fragmentation = speed_up

    def grind(self, max_depth=None, speed_up=False):
        """ Structural fragmentation. Populating dataset.fragments and structure.decompositions. """
        for foo in self.grindProgress(max_depth, speed_up):
            pass


class SARTrainingset(Dataset):
    """ Training set for discrete data: the Dataset given as argument must be already populated. """

    def __init__(self, dataset):
        super().__init__(dataset)
        self.__dict__.update(dataset.__dict__)
        self.type = 'SAR'
        self.target = None
        self.excludeOutliers = False
        self.labels = {structure.target for structure in self.values()}
        self.subset = {label: Dataset() for label in self.labels}
        for structure in self.values():
            # structure.fingerprint = chem.mol2fp(structure.mol)
            self.subset[structure.target].add(structure)
        self.labels = sorted(self.labels, key=lambda key: -len(self.subset[key]))  # sorted for reproducibility
        # self.addHs()  # always?

    def __repr__(self):
        return str(self.counter())
        # return str({label: len(self.subset[label]) for label in self.labels})

    def counter(self):
        return {label: len(self.subset[label]) for label in self.labels}

    def copy(self):
        new = super().copy()
        new.subset = {label: Dataset() for label in new.labels}
        for structure in new.values():
            new.subset[structure.target].add(structure)
        return new

    def add(self, structure):
        super().add(structure)
        self.subset[structure.target].add(structure)
        assert sum(len(dataset) for dataset in self.subset.values()) == len(self)  # DEBUG

    def pop(self, cansmiles):
        self.subset[self[cansmiles].target].pop(cansmiles)  # pop from sub-dataset
        return super().pop(cansmiles)

    def clear_fragmentation(self):
        """ Clear the Dataset and the contained Structures from fragmentation. """
        for dataset in self.subset.values():
            dataset.fragmentation_level = 0
            dataset.fragments = {}
        super().clear_fragmentation()

    def grindProgress(self, max_depth=None, speed_up=False):
        """ Populating dataset.fragments and structure.decompositions. """
        grinder = Grinder()
        for i, structure in enumerate(self.structures()):
            fragments = grinder.grind(structure, max_depth, speed_up)
            for fragment in fragments.values():
                # Binding substructure: superstructures
                self.fragments.setdefault(fragment, []).append(structure)
                self.subset[structure.target].fragments.setdefault(fragment, []).append(structure)
            yield i
        self.fragmentation_level = max_depth
        self.fast_fragmentation = speed_up

    # def binarize(self, labels_dict):
    #     self.labels = sorted(labels_dict.keys())  #, key=lambda key: -len(self.subset[key]))  # sorted for reproducibility
    #     key0, key1 = self.labels
    #     for structure in self.structures():
    #         if structure.target in labels_dict[key0]:
    #             structure.target = key0
    #         else:
    #             # assert structure.target in labels_dict[key1]  # DEBUG!
    #             structure.target = key1
    #     for new_label, old_labels in labels_dict.items():
    #         self.subset[new_label] = set()
    #         for old_label in old_labels:
    #             self.subset[new_label].update(self.subset.pop(old_label))
    #         self.subset[new_label] = Dataset(self.subset[new_label])


class QSARTrainingset(Dataset):
    """ Training set: the Dataset given as argument must be already populated. """

    def __init__(self, dataset):
        super().__init__(dataset)
        self.__dict__.update(dataset.__dict__)
        self.type = 'QSAR'
        self.excludeOutliers = False
        self.iqr = self.lower_fence = self.upper_fence = 0
        self._bases = None
        self._outliers = None
        self._try_make_numeric()
        if len(self) > 1:
            self._stats()

    def _stats(self):
        q1 = self.percentile('target', 25)
        q3 = self.percentile('target', 75)
        self.iqr = q3 - q1
        self.lower_fence = q1 - 1.5 * self.iqr
        self.upper_fence = q3 + 1.5 * self.iqr
        targets = [structure.target for structure in self.values()]
        self.mean = statistics.mean(targets)
        self.abs_dev = sum((abs(target - self.mean) for target in targets)) / len(targets)
        self.stdev = statistics.stdev(targets)

    def clearMols(self):
        for structure in self.values():
            del structure.mol

    def reset(self):
        """ Undo popOutliers(). """
        self._bases = None
        self._outliers = None
        self.excludeOutliers = False

    def popOutliers(self):
        """ Outliers are virtually popped: structures() and get() will ignore it. Revert with reset(). """
        self._bases = Dataset()
        self._outliers = Dataset()
        for structure in self.values():
            if self.lower_fence <= structure.target <= self.upper_fence:
                self._bases.add(structure)
            else:
                self._outliers.add(structure)
        self.excludeOutliers = True
        return self._outliers

    def get(self, cansmiles):
        """ Return the structure with the given cansmiles if present, else None. """
        if self.excludeOutliers:
            return self._bases.get(cansmiles)
        else:
            return super().get(cansmiles)

    def structures(self):
        """ Return an iterable providing a view on the dataset structures. """
        if self.excludeOutliers:
            return self._bases.structures()
        else:
            return super().structures()


class Testset(Dataset):
    """ Test set: the Dataset given as argument must be already populated. """
    prediction_record = {'mode': None, 'predicted by': None, 'prediction': None, 'target': None, 'error': None}

    def __init__(self, dataset):
        super().__init__(dataset)
        self.__dict__.update(dataset.__dict__)
        self.unpredicted = Dataset(self)
        self.predictions = {}

    def _setPrediction(self, cansmiles, mode, predictor, prediction):
        """ Predicted structures are virtually popped (structures() will ignore them). Revert with reset(). """
        if isinstance(prediction, float):
            prediction = round(prediction, DECIMAL)
        record = Testset.prediction_record.copy()
        record['mode'] = mode
        record['predicted by'] = predictor
        record['prediction'] = prediction
        self.predictions[cansmiles] = record
        self.unpredicted.pop(cansmiles)

    def reset(self):
        """ Clear the test set by any previous prediction. """
        self.unpredicted = Dataset(self)
        self.predictions.clear()

    def predictIdentities(self, trainingset):
        """ Detect, predict and (virtually) pop structures identical to structures in the training set. """
        identities = []
        for cansmiles in self.keys() & trainingset.keys():
            identity = trainingset[cansmiles]
            prediction = identity.target
            self._setPrediction(cansmiles, 'identity', identity, prediction)
            identities.append(identity)
        return identities

    def predictOutliers(self, trainingset):  # QSAR only
        """ Detect, predict and (virtually) pop structures similar to structures in the training set. """
        outliers = []
        for cansmiles, structure in list(self.unpredicted.items()):
            similarOutliers = structure.getSimilar(trainingset._outliers)
            if similarOutliers:
                prediction = statistics.median([structure.target for structure in similarOutliers])
                self._setPrediction(cansmiles, 'similar to outliers', similarOutliers, prediction)
                outliers.append(structure)
        return outliers

    def setPrediction(self, cansmiles, predictor):
        """ Pop and predict the relative structure. (Populate self.predictions) """
        prediction = predictor.prediction
        self._setPrediction(cansmiles, predictor.type(), predictor, prediction)
        return self.get(cansmiles)

    def validate(self, numeric=False):
        if numeric:
            self._try_make_numeric()  # can raise ValueError
        for cansmiles, record in self.predictions.items():
            record['target'] = self.get(cansmiles).target

            if not numeric:  # SAR
                # target = self.binarize[target]
                # prediction = self.binarize[prediction]
                if record['target'] == record['prediction']:
                    error = 0
                else:
                    error = 1
            else:  # QSAR
                error = record['prediction'] - record['target']
            record['error'] = error
        return True

    def structures(self, all=False):
        """ Return an iterable providing a view on the dataset structures.
            Predicted structures are ignored, unless all=True."""
        if all:
            return self.values()
        else:
            return self.unpredicted.structures()

    def r2(self):
        """ Return the coefficient of determination, if possible. """
        if not self.predictions:
            return
        targets = []
        errors = []
        predictions = []
        for structure in self.values():
            targets.append(structure.target)
            if prediction_dict := self.predictions.get(structure.cansmiles):
                predictions.append(prediction_dict['prediction'])
                errors.append(prediction_dict['error'])
        mean_target = sum(targets) / len(targets)  # mean of ALL targets (not only the actually predicted)
        total_sum_of_squares = sum([(pred - mean_target) ** 2 for pred in predictions])
        residual_sum_of_squares = sum([error ** 2 for error in errors])
        try:
            return round(1 - residual_sum_of_squares / total_sum_of_squares, 2)
        except ZeroDivisionError:
            return


def loadDataset(path, idKey=None, smilesKey=None, targetKey=None):  #, isDiscrete=False):
    dataset = Dataset()
    dataset.path = path.name
    ext = path.suffix
    if ext == '.csv':
        generator = _readcsv
        generator.smilesHeader = smilesKey
    elif ext == '.sdf':
        generator = chem.readSDF
    else:
        raise
    index = 0
    errors = 0
    for mol, data in generator(str(path.resolve())):
        index += 1
        if not mol:   # CSV only
            # error_id = str(smilesKey or idKey or 'index') + ' ' + (data.get(smilesKey) or data.get(idKey) or str(index))
            # print(" Chemical Library ERROR at " + error_id)
            errors += 1
            continue
        if smilesKey is not None:
            smiles = data[smilesKey]
        else:
            smiles = None
        mol = chem.removeHs(mol)
        structure = Structure(mol, smiles)
        if idKey is not None:
            structure.S_id = data[idKey]
        else:
            structure.S_id = '#' + str(index)  # str(index)
        if targetKey is not None:
            structure.target = data.get(targetKey)
        dataset.add(structure)
    if errors:
        print(f"\n WARNING: {errors} structure(s) discarded\n")
    return dataset


def _dumpDataset(dataset):
    return {s.smiles: (s.S_id, s.smiles, s.target) for s in dataset.values()}


def _restoreDataset(dumped_dataset):
    dataset = Dataset()
    for cansmiles, record in dumped_dataset.items():
        s_id, smiles, target = record
        mol = chem.smi2mol(smiles)
        structure = Structure(mol, smiles)
        structure.smiles = smiles
        structure.S_id = s_id
        structure.target = target
        dataset.add(structure)
    return dataset

# Use pandas: pandas.read_csv('data.csv',sep=None)  with spe=None tries to detect the delimiter automatically
def _readcsv(path):
    for index, row in enumerate(csv.reader(open(path))):
        if index == 0:
            headers = row
            continue
        data = dict(list(zip(headers, row)))
        smiles = data[_readcsv.smilesHeader]
        try:
            mol = chem.smi2mol(smiles)
            yield mol, data
        except IOError as error:
            print(" #{} Chemical Library ERROR: failed to parse SMILES '{}'".format(index, smiles))
            yield None, error


# def printHierarchy(structure):  # DEBUG only
#     print()
#     print('Hierarchy of %s:\n' % structure)
#     for i in range(len(structure.factorizations)):
#         print(structure.factorizations[i])
#         print()
#     print()