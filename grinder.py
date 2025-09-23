#  DeltaQSAR, a knowledge discovery toolkit for automatic SAR/QSAR rules induction.
#  Copyright (C) 2021  Thomas Ferrari


import chem
from parameters import MAX_DEC_FAST, MAX_DEC


TwoCharsElements = {'Br', 'Cl'}
OneCharElements = {'C', 'N', 'O', 'P', 'S', 'B', 'F', 'I'}
AromaticElements = {'c', 'n', 'o', 'p', 's'}
Bonds = {'-', '=', '#', ':', '.', '~', '/', '\\'}
SplittingChars = OneCharElements | set([e[0] for e in TwoCharsElements]) | AromaticElements
Break = '*'

# GLOBAL
cansmiles_dict = {}  # {cansmiles: Fragment OR Structure} It's a log to retrieve Structure OR Fragment FROM cansmiles


class Fragment:

    def __init__(self, mol, smiles=None, cansmiles=None):
        super().__init__()
        self.smiles = smiles
        self.cansmiles = cansmiles or chem.mol2smi(mol)
        self.cansmiles_unstarred = chem.getCansmilesUnstarred(mol)  # RDKit add Hs if needed, OB no
        # if self.cansmiles_unstarred.count('H') != self.cansmiles.count('H'):
        #     print('mismatch', self.cansmiles_unstarred, self.cansmiles)
        self.atoms = chem.numAtoms(mol)
        self.wildcards = self.cansmiles.count('*')
        # self.bond_map = chem._get_bond_map(mol, self.smiles or self.cansmiles)  # RDKit ONLY
        cansmiles_dict[self.cansmiles] = self

    def __repr__(self):
        return self.cansmiles

    def __lt__(self, other):
        return str(self) < str(other)

    def size(self):
        return self.atoms, -self.wildcards


class Grinder:

    def __init__(self):
        # self.cansmilesDict = dataset.fragments.copy()  # speeds up if dataset already fragmented
        self.smilesFrag_dict = {}
        self.couples = {}

    # FRAGMENT FACTORY: can raise IOError
    def _makeCouple(self, head, tail):
        couple = []
        for smilesFrag in head, tail:
            try:
                frag = self.smilesFrag_dict[smilesFrag]
            except KeyError:
                mol = chem.smi2mol(smilesFrag)  # CAN RAISE IOError, handled in the function call
                cansmiles = chem.mol2smi(mol)
                try:
                    frag = cansmiles_dict[cansmiles]
                except KeyError:
                    frag = Fragment(mol, smilesFrag, cansmiles)
                    # assert frag.atoms != 0  # DEBUG
                    # self.cansmilesDict[cansmiles] = frag
                self.smilesFrag_dict[smilesFrag] = frag
            couple.append(frag.cansmiles)
        assert len(couple) == 2  # DEBUG
        couple.sort()
        return couple

    def _getCouples(self, cansmiles):

        if cansmiles in self.couples:
            return self.couples[cansmiles]

        # To use bond_map, the original input SMILES must be fragmented (WARNING: can contain H and @)
        # smiles = cansmiles_dict[cansmiles].smiles
        smiles = cansmiles

        # bond_map = cansmiles_dict[cansmiles].bond_map
        # bond_index = len(bond_map)
        couples = []
        inBrackets = False
        parenthesisStack = []
        index = len(smiles)
        while index > 1:  # Traverse SMILES string from right to left seeking for a valid index where to split
            index -= 1
            char = smiles[index]
            if inBrackets and char != '[':
                continue  # don't split inside square brackets
            elif char in SplittingChars:
                pass  # EXIT
            elif char == ')':
                parenthesisStack.append(index)
                continue
            elif char == '(':
                parenthesisStack.pop()
                continue
            elif char == ']':
                inBrackets = True
                continue
            elif char == '[':
                inBrackets = False  # EXIT: split before '['
            elif char == '*':  # I can remove it if I put * in Splitting chars
                # bond_index -= 1
                continue
            else:
                continue

            # bond_index -= 1  # for bond_map
            # if not bond_map[bond_index]:  # RDKit ONLY !!
            #     continue  # do not break

            # Get 2 SMILES fragments (BRAKE 1 BOND)
            left_char = smiles[index - 1]
            if left_char in Bonds:
                bond = left_char
            else:
                bond = ''
            if parenthesisStack:
                closingParenthesesIndex = parenthesisStack[-1]
                head = smiles[:index] + Break + smiles[closingParenthesesIndex:]  # left bit + Break + right bit
                tail = Break + bond + smiles[index:closingParenthesesIndex]  # Break + bond + middle bit
            else:
                head = smiles[:index] + Break  # left bit + Break
                tail = Break + bond + smiles[index:]  # Break + bond + right bit

            if chem.invalidSmiles(head, tail):  # validity check (without using bond_map)
                continue

            try:
                # print("\nMAKING COUPLE:", head, tail)
                couple = self._makeCouple(head, tail)
            except IOError:
                continue
            assert couple  # DEBUG
            if couple not in couples:
                couples.append(couple)

        self.couples[smiles] = couples
        return couples

    def _getDeeper(self, prev_level):  # a decomposition level is a set of tuples of fragments
        next_level = set()
        for decomposition in prev_level:
            for smiles in decomposition:
                unchanged = [*decomposition]
                unchanged.remove(smiles)
                # (CHK self.couples should be moved here)
                for new_decomposition in (unchanged + couple for couple in self._getCouples(smiles)):
                    assert new_decomposition  # DEBUG                       TO BE ROMOVED
                    next_level.add(tuple(sorted(new_decomposition)))  # SORTING (it should be unordered: i.e., Counter)
        return next_level

    def grind(self, structure, max_depth, speed_up=False):
        if structure.decompositions is None:
            cansmiles_dict[structure.cansmiles] = structure
            structure.decompositions = [{(structure.cansmiles,)}]
            # structure.bond_map = chem._get_bond_map(structure.mol, structure.smiles)
        fragments = {}
        depth = len(structure.decompositions) - 1
        while max_depth is None or depth < max_depth:  # There is also depth 0
            last_level = structure.decompositions[-1]
            max_dec = MAX_DEC
            if speed_up:
                max_dec = MAX_DEC_FAST
            # if depth >= 1 and len(structure.decompositions[1]) > 15:  # structure.decompositions[-1] > 15:  len(last_level) > MAX_DEC:
            if depth >= 1 and len(last_level) > max_dec:
                break  # stop if last level already broke too many bonds
            new_level = self._getDeeper(last_level)
            if not new_level:
                break
            structure.decompositions.append(new_level)
            depth += 1
            # non serve un dizionario, potrei semplicemente appiattire la sequenza innestata
            for decomposition in new_level:
                for cansmiles in decomposition:
                    fragments[cansmiles] = cansmiles_dict[cansmiles]
                    # fragments[cansmiles]._superstructures.setdefault(structure, {}).setdefault(depth, []).append(decomposition)
        # fill missing levels
        while max_depth is not None and len(structure.decompositions) < max_depth + 1:
            structure.decompositions.append(set())
        return fragments


def clearContext():
    global cansmiles_dict
    cansmiles_dict.clear()


    # def _getAtomsNum(self, smilesFrag):  # Tested on 10000 smiles
    #     smiles = smilesFrag
    #     c = 0
    #     s = ''
    #     inBrackets = False
    #     #for wildcard in Grinder.WildCards:
    #         #smiles = smiles.replace(wildcard, '')
    #     smiles = smiles.replace(self.Break, '')
    #     for char in smiles:
    #         if char == '[':
    #             inBrackets = True
    #             c += 1
    #             continue
    #         if char == ']':
    #             inBrackets = False
    #             continue
    #         if not inBrackets and char in letters:
    #             s = s + char
    #     for chars in Grinder.TwoCharsElements:
    #         c += s.count(chars)
    #         s = s.replace(chars, '')
    #     atoms = len(s) + c
    #     return atoms
