import chem
from parameters import DECIMAL


class StructuralAlert:

    SUMMARY_HEADERS = "ID\tPPV\tTP\tFP\ttarget\tSMILES\n"  # "ID\tPPV\tabs_PPV\tTP\tFP\ttarget\tSMILES\n"

    @staticmethod
    def type():
        return 'SA'

    def __init__(self, fragment, targets, nostar=False, addHs=False, smarts=False):
        self.fragment = fragment
        self.explicit_Hs = False
        self.wildcards = fragment.wildcards
        self.atoms = fragment.atoms  # THE CUSTOM PARAMETER IS atoms - wildcards, shouldn't be like that here?
        # self.fingerprint = chem.mol2fp(self.smarts)  # RDKit ONLY (smarts is mol, but they don't match the same)
        if nostar:
            mol = chem.smi2mol(self.fragment.cansmiles_unstarred)
            self.fingerprint = chem.mol2fp(mol)
            # Chemical libraries (RDKit) can add Hs in 'cansmiles_unstarred' to fix the valence
            # (e.g., CS(=O)=O becomes C[SH](=O)=O ). But SMARTS are just patterns and valence shouldn't be considered.
            # Therefore, in case of added Hs, the SMARTS is derived from 'cansmiles' after removing the wildcards (*).
            # (In RDKit I could probably use MolFromSmarts() instead of MolFromSmiles()...)
            if self.fragment.cansmiles.count('H') == self.fragment.cansmiles_unstarred.count('H'):
                self.smarts_string = self.fragment.cansmiles_unstarred  # canonical SMILES
            else:
                self.smarts_string = chem.removeWildcards(fragment.cansmiles)  # NON-canonical SMILES
            self.smarts_string = chem.remove_nH(self.smarts_string)
            self.atoms -= self.wildcards
            self.wildcards = 0
        elif addHs:
            mol = chem.smi2mol(self.fragment.cansmiles)
            self.fingerprint = chem.mol2fp(mol, remove_star=True)
            self.explicit_Hs = True
            self.smarts_string = chem.mol2smi_H(mol)
        # elif smarts:
        #     assert chem.LIBRARY == 'RDKit'
        #     mol = chem.smi2mol(self.fragment.cansmiles_unstarred)
        #     self.fingerprint = chem.mol2fp(mol, remove_star=True)
        #     self.smarts_string = chem.rdkit.MolToSmarts(mol)
        #     self.atoms -= self.wildcards
        #     self.wildcards = 0
        else:
            self.smarts_string = fragment.cansmiles
        self.smarts = chem.getSmarts(self.smarts_string)
        self.precision = self.train_precision = self.F_score = 0
        self.index = None
        self.prediction = None
        self.train_hits = set()
        self._hits = {label: set() for label in targets}
        self.potential_exclusions = set()
        self.exclusions = []
        self.excluded_hits = {True: set(), False: set()}
        # self.generalizedBy = []

    def __repr__(self):
        if self.exclusions:
            return self.smarts_string + '  - EXC ' + str([ex.smarts_string for ex in self.exclusions])
        return self.smarts_string

    def __lt__(self, other: 'StructuralAlert'):
        return (self.F_score, self.train_precision, self.size(), str(self)) < \
               (other.F_score, other.train_precision, other.size(), str(other))

    def getID(self):
        return 'SA' + str(self.index)

    def summary(self):
        d = DECIMAL
        tab = 8
        # summary = f"  {self.getID()}\t{self.precision:.{d}f} \t{self.train_precision:.{d}f}" \
        summary = f"  {self.getID()}\t{self.precision:.{d}f}"\
                  f"\t{self.countTrue()}\t{self.countFalse()}" \
                  f"\t{self.prediction[:tab - 1]}\t{self}".expandtabs(tab)
        # if self.exclusions:
        #     summary += f"\t(before exclusions  T = {len(self.allTrueHits())}\tF = {len(self.allFalseHits())})"
        return summary

    def size(self):
        return self.atoms, -self.wildcards

    def activate(self, target):
        self.prediction = target
        self.train_precision = self.precision

    def getMatches(self, structures, prescreen_FP=False):
        matched = set()
        for structure in structures:
            # Pre-screen with fingerprints to speed up
            if prescreen_FP and not chem.fpIsSubset(self.fingerprint, structure.fingerprint):
                continue
            mol = structure.mol
            if chem.match(mol, self.smarts):
                if any(chem.match(mol, exclusion.smarts) for exclusion in self.exclusions):
                    continue
                matched.add(structure)
        return matched

    def match(self, structures, prescreen_FP=False):  # called after __init__
        self.train_hits = self.getMatches(structures, prescreen_FP=prescreen_FP)
        for hit in self.train_hits:
            self._hits[hit.target].add(hit.cansmiles)

    def yieldHits(self):
        for target, matched in self._hits.items():
            for cansmiles in matched:
                if cansmiles not in self.excluded_hits[target is self.prediction]:
                    yield cansmiles, target

    def removeHits(self, other: 'StructuralAlert'):
        for key in self._hits:
            self._hits[key] -= other._hits[key] - other.excluded_hits[key is other.prediction]

    def countTrue(self, target=None):
        return len(self.allTrueHits(target)) - len(self.excluded_hits[True])

    def countFalse(self, target=None):
        return len(self.allFalseHits(target)) - len(self.excluded_hits[False])

    def allTrueHits(self, target=None):
        return self._hits[target or self.prediction]

    def allFalseHits(self, target=None):
        return set.union(*(self._hits[key] for key in self._hits if key != (target or self.prediction)))

    def calcMetrics(self, dataset, beta=0.1, target=None):
        target = self.prediction or target
        trueMatches = self.countTrue(target)
        if trueMatches == 0:
            self.precision = self.F_score = 0
            return
        self.precision = trueMatches / (trueMatches + self.countFalse(target))
        recall = trueMatches / len(dataset.subset[target])
        if beta == 0:  # if beta == 0 => F_score = precision: it doesn't take into account sensitivity at all
            self.F_score = (self.precision, recall)
            return
        factor = beta * beta
        self.F_score = (1 + factor) * (self.precision * recall) / (factor * self.precision + recall)

    def setPotentialExcusions(self, alerts):
        for other in alerts:
            # assert other.atoms >= self.atoms  # DEBUG
            if self.explicit_Hs and not other.explicit_Hs:
                continue
            if self.prediction != other.prediction:
                if self.train_hits >= other.train_hits:  # fingerprint slower
                    # if chem.match(other.mol, self.smarts):  # waste of time
                    self.potential_exclusions.add(other)

    def evalAsExclusionsOf(self, base_alert, min_samples, beta):
        self._correctly_removed = self.allFalseHits(base_alert.prediction) - base_alert.excluded_hits[False]
        self._mistakenly_removed = self.allTrueHits(base_alert.prediction) - base_alert.excluded_hits[True]
        F = len(self._mistakenly_removed)
        T = len(self._correctly_removed)
        if T < min_samples or T + F + len(base_alert.excluded_hits[False]) + len(base_alert.excluded_hits[True])\
                >= (len(base_alert.allTrueHits()) + len(base_alert.allFalseHits())) / 2:
                # >= len(base_alert.train_hits) / 2:
            self._exc_F_score = 0
            return False
        precision = T / (T + F)
        recall = T / base_alert.countFalse()
        if beta == 0:  # if beta == 0 => F_score = precision: it doesn't take into account sensitivity at all
            self._exc_F_score = (precision, recall)
        else:
            factor = beta * beta
            self._exc_F_score = (1 + factor) * (precision * recall) / (factor * precision + recall)
        return precision >= base_alert.precision

    def resetExclusions(self):
        self.exclusions.clear()
        self.excluded_hits[True].clear()
        self.excluded_hits[False].clear()

    def assignExclusions(self, min_samples, beta):
        self.resetExclusions()
        if self.countFalse() == 0: return

        potential_exclusions = []
        for alert in sorted(self.potential_exclusions, key=lambda x: x.smarts_string):  # reproducibility
            if alert.evalAsExclusionsOf(self, min_samples, beta):
                potential_exclusions.append(alert)

        while potential_exclusions:
            potential_exclusions.sort(key=lambda ex: (ex._exc_F_score, ex.size()))
            exclusion = potential_exclusions.pop()
            self.excluded_hits[False] |= exclusion._correctly_removed
            self.excluded_hits[True] |= exclusion._mistakenly_removed
            self.exclusions.append(exclusion)
            if self.countFalse() == 0 or len(self.exclusions) == 2:
                break
            # UPDATE OTHERS
            for i in reversed(range(len(potential_exclusions))):
                other = potential_exclusions[i]
                if not other.evalAsExclusionsOf(self, min_samples, beta):
                    potential_exclusions.pop(i)
        # Cleaning
        for alert in self.potential_exclusions:
            del alert._correctly_removed, alert._mistakenly_removed, alert._exc_F_score
        # if self.exclusions:
        #     assert self.excluded_hits[False] == set.union(*(exclusion.allFalseHits(self.prediction) for exclusion in self.exclusions))
        #     assert self.excluded_hits[True] == set.union(*(exclusion.allTrueHits(self.prediction) for exclusion in self.exclusions))

    # def interferenceCheck(self, others):
    #     for other_rule in others:
    #         assert other_rule.index > self.index  # DEBUG
    #         if self.prediction == other_rule.prediction and self.train_hits.issubset(other_rule.train_hits):
    #             self.generalizedBy.append(other_rule.index)
