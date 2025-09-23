#  DeltaQSAR, a knowledge discovery toolkit for automatic SAR/QSAR rules induction.
#  Copyright (C) 2021  Thomas Ferrari


import rdkit.Chem as rdkit  # import only what you need: is the exec any lighter?
from rdkit.Chem import PandasTools
from openbabel import pybel
import re


# Library-specific utils
OB = {'IN_opt': {'a': True},  # 'a': preserve aromaticity
      'OUT_opt': {'c': None, 'i': None},  # 'c': canonical form, 'i': do not include isotopic or chiral markings
      'biphenyl_fp' : set(pybel.readstring('smi', 'c1ccccc1c1ccccc1').calcfp().bits),
     }
# RDKit = {'asterisk': rdkit.MolFromSmiles('*')}

Digits = '123456789'
Bonds = {'-', '=', '#', ':', '.', '~', '/', '\\'}
NonAtomChars = {'*'}.union(Bonds, Digits, '()[]')

wildcard_substrings = []
for bond in Bonds:
    wildcard_substrings.extend(('*' + bond, bond + '*'))
wildcard_substrings.extend(('*', '()'))


def setLibrary(lib_slug='OB'):
    global LIBRARY
    LIBRARY = lib_slug
    if LIBRARY == 'OB':
        pybel.ob.cvar.obErrorLog.StopLogging()  # ob.obErrorLog.SetOutputLevel(0)
    else:  # otherwise RDKit
        from rdkit import RDLogger
        RDLogger.logger().setLevel(RDLogger.CRITICAL)


def removeHs(mol):
    if LIBRARY == 'OB':
        mol.removeh()
        return mol
    return rdkit.RemoveHs(mol, sanitize=False)


def numAtoms(mol):
    if LIBRARY == 'OB':
        return mol.OBMol.NumAtoms()
    return mol.GetNumAtoms()


def removeWildcards(smiles):
    for substring in wildcard_substrings:
        smiles = smiles.replace(substring, '')
    return smiles


def remove_nH(smiles):
    # Here I'd like to remove EVERY H to get more generic SMARTS, but using regex is risky...
    # Should I remove H if there is a charge? What about chirality, isotopic specifications, [2H], etc?
    # For now, I just remove [nH] (the majority)
    return smiles.replace('[nH]', 'n')
    # return re.sub(r"\[(.{1,2})H]", r"\[\1\]", smiles)  # remove square brackets only if organic

# def removeBrackets(smiles):
#     return re.sub(r'\[([B,C,N,O,P,S,F,Cl,Br,I,c,n,o,s,p])\]', r'\1', smiles)


def smi2mol(smiles):
    if LIBRARY == 'OB':
        mol = pybel.readstring('smi', smiles, opt=OB['IN_opt'])  # can raise IOError
        return mol
    mol = rdkit.MolFromSmiles(smiles, sanitize=False)  # with sanitize=False explicit Hs in input are explicit in output
    if mol is None:
        raise IOError
    if mol.NeedsUpdatePropertyCache():
        mol.UpdatePropertyCache()
    return mol


def mol2smi(mol):
    if LIBRARY == 'OB':
        return mol.write(opt=OB['OUT_opt']).rstrip()
    return rdkit.MolToSmiles(mol, isomericSmiles=False)


def mol2smi_H(mol):
    if LIBRARY == 'OB':
        return to_H_smiles(mol.OBMol).rstrip()
    return rdkit.MolToSmiles(mol, allHsExplicit=True, isomericSmiles=False)


def mol2fp(mol, remove_star=False):
    if remove_star:
        mol = getMolUnstarred(mol)
    if LIBRARY == 'OB':
        return set(mol.calcfp().bits)  # FP 140288 FN 64
    return rdkit.RDKFingerprint(mol)  # FP 123596 FN 4
    # PatternFingerprint should be optimized for sub-matching, but its performance sucks
    # return rdkit.PatternFingerprint(mol)   # FP 203186 FN 8


def fpIsSubset(fingerprint_sub, fingerprint_sup):
    if LIBRARY == 'OB':
        # OpenBabel PATCH FROM SARpy 1.0
        # Fingerprints have some FPs, especially with biphenyl (the single bond can be aromatic): since they are used
        # only to pre-screen before SMARTS matching, with biphenyl pre-screen is skipped
        if fingerprint_sub.issuperset(OB['biphenyl_fp']):
            return True
        return fingerprint_sub <= fingerprint_sup
    return fingerprint_sub & fingerprint_sup == fingerprint_sub  # RDkit: type(fingerprint) = ExplicitBitVect


def getSmarts(smiles):
    if LIBRARY == 'OB':
        try:
            smarts = pybel.Smarts(smiles)
        except:  # IOError or OSError ?
            smarts = None
    else:
        smarts = rdkit.MolFromSmarts(smiles)
    if smarts is None:
        print(f" Invalid SMARTS: {smiles}")  # DEBUG
    return smarts


def match(mol, smarts):
    if LIBRARY == 'OB':
        if smarts is None:
            return False
        return smarts.obsmarts.Match(mol.OBMol, True)
    return mol.HasSubstructMatch(smarts)


def getCansmilesUnstarred(mol):
    if LIBRARY == 'OB':
        atoms = mol.atoms
        valid_atoms_indexes = [atoms.index(x) + 1 for x in atoms if x.atomicmass != 0]
        opt = OB['OUT_opt'].copy()
        opt['F'] = ' '.join(str(e) for e in valid_atoms_indexes) + ' '
        return mol.write(opt=opt).rstrip()  # H are NOT added even if needed!
    return mol2smi(getMolUnstarred(mol))  # H are added if needed!
# RDKit adds [.H] to fix the valence, while OB no
# [.H] is good for decompositions (to check if it is a base mol), but it's not good for SMARTS
# Now Hs are removed from SMARTS in the StructuralAlerts constructor, so I could implement OB in the same way of RDKit


def getMolUnstarred(mol):
    if LIBRARY == 'OB':
        return smi2mol(getCansmilesUnstarred(mol))
    return rdkit.DeleteSubstructs(mol, smi2mol('*'))  # H are added if needed!
    # return rdkit.DeleteSubstructs(mol, RDKit['asterisk'])  # H are added if needed!


def readSDF(path):
    # apparently it never raises errors (it just quit parsing)
    if LIBRARY == 'OB':
        for mol in pybel.readfile('sdf', path):
            if mol.atoms:
                yield mol, dict(mol.data)
    else:
        for row in PandasTools.LoadSDF(path).iloc:
            yield row.pop('ROMol'), row
            # rebuilding mol from cansmiles so that mol and SMILES atoms/bonds are in the same order (for bond_map)
            # row['SMILES'] = mol2smi(row.pop('ROMol'))
            # yield smi2mol(row['SMILES']), row


def hasBrokenCycles(smiles):
    # From Daylight: "Structures that require more than 10 ring closures to be open at once are exeedingly rare"
    for d in Digits:  # counting the actual cycles doesn't speed up
        if smiles.count(d) % 2:
            return True


def invalidSmiles(smiles_part_1, smiles_part_2):
    if hasBrokenCycles(smiles_part_1) or \
            NonAtomChars.issuperset(smiles_part_1) or NonAtomChars.issuperset(smiles_part_2):  # dummy smiles
        return True


# def numBonds(mol):
#     if LIBRARY == 'OB':
#         return mol.OBMol.NumBonds()
#     return mol.GetNumBonds()
#
#
# def getBond(mol, idx):
#     if LIBRARY == 'OB':
#         return mol.OBMol.GetBond(idx)
#     return mol.GetBondWithIdx(idx)


# def generalizeAlogens(smiles):
#     return re.sub('(Cl|Br|F|I)', '[Cl,Br,F,I]', smiles)


# OpenBabel only #######################################################################################################
# From Andrew at dalke@dalkescientific.com

# Use Open Babel to generate a SMILES string with the implicit
# hydrogens explicitly included in an [atom] term.

# For example: "[CH3:1]NCF" -> "[CH3:1][NH][CH2][F]"

# Open Babel does not appear to support this directly.

# Instead, assign an atom class to each atom and use the 'a' output
# conversion option to include the atom class in the output. This has
# the side-effect of putting every atom term inside of []s, which
# means the implicit hydrogens are explicitly included. (If the SMILES
# atom did not have an atom class then assign it a class of 0,
# otherwise add one to the existing atom class value.)

# The generated SMILES contains incorrect atom class values so
# post-process the SMILES string (using a regular expression) to
# adjust the atom class either to remove it if the atom didn't
# originally have an atom class, or to decrease it by one.

from openbabel import openbabel as ob
obconv = ob.OBConversion()
obconv.SetInFormat("smi")
obconv.SetOutFormat("smi")
obconv.AddOption("a", obconv.OUTOPTIONS)  # Output atomclass like [C:2], if available


# Return the atom class if present, else return None
def _get_atom_class(atom):
    atom_class = atom.GetData("Atom Class")
    if atom_class is None:
        return None
    else:
        # Can be 0 or larger
        return ob.toPairInteger(atom_class).GetGenericValue()


# Set the atom class. Assumes there is no existing "Atom Class" on the atom.
def _set_atom_class(atom, value):
    atom_class = ob.OBPairInteger()
    atom_class.SetAttribute("Atom Class")
    atom_class.SetValue(value)
    atom.CloneData(atom_class)


# Match the atom class term
_atom_class_pattern = re.compile(r":([0-9]+)\]")


def _adjust_atom_class(m):
    # 'm' is a match from _atom_class_pattern, like ":0]" and ":13]"
    if m.group() == ":0]":
        # The atom didn't start with an atom class so remove it completely.
        return "]"

    # The atom had an atom class so adjust it
    n = int(m.group(1)) - 1
    return f":{n}]"


def to_H_smiles(obmol):
    obmol = ob.OBMol(obmol)  # make a copy because I will modify the molecule
    atoms = [obmol.GetAtom(i + 1) for i in range(obmol.NumAtoms())]

    for atom in atoms:
        atom_class = _get_atom_class(atom)
        if atom_class is None:
            # I'll post-process the output to remove the value.
            atom_class = 0
        else:
            # Need to remove 'Atom Class' so _set_atom_class() works.
            atom.DeleteData("Atom Class")
            # I'll post-process the output to correct the value.
            atom_class += 1

        _set_atom_class(atom, atom_class)

    s = obconv.WriteString(obmol)

    # The output include "\t" followed by the identifier.
    smiles, _, id_and_newline = s.partition("\t")

    # Replace atom class terms in the SMILES
    new_smiles = _atom_class_pattern.sub(_adjust_atom_class, smiles)
    return f"{new_smiles}\t{id_and_newline}"


# RDKIT only ###########################################################################################################

# # ps = rdkit.SmilesParserParams()
# # ps.removeHs = True
# # ps.sanitize = False

# def getMatches(mol, smarts):  # ONLY for bond_map
#     if LIBRARY == 'OB':
#         return smarts.findall(mol)
#     return mol.GetSubstructMatches(smarts)
#
# # heteroatoms = {'N', 'O', 'P', 'S', 'Cl', 'Br', 'F'}
# # FRAGMENTATION EXCLUSIONS
# # BONDS IN hetero-C=[N/O/S]: '[N,O,P,S,Cl,Br,F]~C=[N,O,S]'
# # BONDS BETWEEN HETEROATOMS: '[N,O,P,S,Cl,Br,F]~[N,O,P,S,Cl,Br,F]'
# functional_groups = ()
# # functional_groups = ('[#7,#8,#15,#16,Cl,Br,F]~[#6]=[#7,#8,#16]', '[#7,#8,#15,#16,Cl,Br,F]~[#7,#8,#15,#16,Cl,Br,F]')
# smarts_list = [getSmarts(smarts) for smarts in functional_groups]
#
#
# def isInRing(bond):
#     return bond.IsInRing()
#
#
# def hasWildcard(bond):
#     for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
#         if LIBRARY == 'OB':
#             if atom.GetType() == 'Du':
#                 return True
#         else:
#             if atom.GetSymbol() == '*':
#                 return True
#
#
# def _get_bond_map(mol, smiles):  # RDKit ONY!
#     if LIBRARY == 'OB':
#         return [True for i in range(numBonds(mol))]
#     no_break_bonds = {atom_idx - 1 for smarts in smarts_list
#                                        for match in getMatches(mol, smarts)
#                                            for atom_idx in sorted(match)[1:]}
#     bond_map = []
#     idx = 0
#     # Exclude ring closures (digit): the last bonds in RDkit
#     valid_idx = numBonds(mol) - (sum(c.isdigit() for c in smiles) - smiles.count('%')) // 2
#     while idx < valid_idx:
#         current_bond = getBond(mol, idx)
#         if isInRing(current_bond):  # EXCLUDE CYCLES
#             bond_map.append(False)
#         elif hasWildcard(current_bond):
#             bond_map.append(False)
#         elif idx in no_break_bonds:
#             bond_map.append(False)
#         else:
#             bond_map.append(True)
#         idx += 1
#     return bond_map
#
#
#     # elif heteroatoms.issuperset(atoms):  # EXCLUDE BONDS BETWEEN HETEROATOMS
#     #     # print("heteroatom bond")
#     #     bond_map.append(False)
#
#     # def in_functional_group(atoms):
#     #     if heteroatoms.intersection(atoms) and 'C' in atoms:  # bond heteroatom - C
#     #         carbon = atoms.pop('C')
#     #         current_bond_type = current_bond.GetBondTypeAsDouble()
#     #         for other_bond in carbon.GetBonds():
#     #             if other_bond.GetIdx() == current_bond.GetIdx():
#     #                 continue
#     #             other_atom_idx = other_bond.GetOtherAtomIdx(carbon.GetIdx())
#     #             other_atom = mol.GetAtomWithIdx(other_atom_idx)
#     #             # C double bonded with NOS has also bond with hetero
#     #             if current_bond_type == 2:
#     #                 if other_atom.GetSymbol() in heteroatoms:
#     #                     heteroatom = atoms.popitem()[1]
#     #                     # print("functional group:", heteroatom.GetSymbol(), carbon.GetSymbol(),
#     #                     #       'in double plus bond with', other_atom.GetSymbol())
#     #                     return True
#     #             # C has also double bond with NOS
#     #             elif other_bond.GetBondTypeAsDouble() == 2:
#     #                 if other_atom.GetSymbol() in "NOS":
#     #                     heteroatom = atoms.popitem()[1]
#     #                     # print("functional group:", heteroatom.GetSymbol(), carbon.GetSymbol(), 'in double with',
#     #                     #       other_atom.GetSymbol())
#     #                     return True
