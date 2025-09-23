#  DeltaQSAR, a knowledge discovery toolkit for automatic SAR/QSAR rules induction.
#  Copyright (C) 2021 - 2023  Thomas Ferrari

import json
from pathlib import Path
from time import time

import chem
import datasets
from datasets import QSARTrainingset, Testset, SARTrainingset
from QSAR_model import QSARModel
from QSAR_cluster_model import ClusterModel
from SAR_model import SARModel
import output
import delta_DB
from parameters import *  # DECIMAL, Progress  # only if GUI usage


def loadDataset(path, targetKey=None, idKey=None, smilesKey=None):  # , excludeOutliers=False, isDiscrete=False):
    print("\n\n Loading dataset...\n")
    dataset = datasets.loadDataset(path, idKey, smilesKey, targetKey)
    print("  Read %s molecular structures" % len(dataset))
    return dataset


def makeQSARTrainingset(dataset):  # , excludeOutliers=False):
    try:
        trainingset = QSARTrainingset(dataset)
    except ValueError:  # not numeric
        raise
    if trainingset.iqr == 0:
        print("\n ERROR: endpoint values not valid")
        return
    trainingset.clearMols()  # in QSAR is not required after fingerprints are calculated
    if 0 < len(dataset) < 500: print("  (probably too small as a training set)")
    # print("\n Numeric endpoint: QSAR numeric prediction")
    return trainingset


def makeSARTrainingset(dataset, target=None):
    trainingset = SARTrainingset(dataset)
    if 0 < len(dataset) < 500: print("  (probably too small as a training set)")
    trainingset.target = target
    if len(trainingset.labels) < 2:
        print("\n ERROR: endpoint values not valid")
        return
    print("\n   > CLASSES: ", trainingset)
    return trainingset


def fragmentize(dataset, depth=None, speed_up=False, pop_outliers=False, clear=True):
    print("\n\n Fragmenting...\n")
    start = time()
    if clear:
        dataset.clear_fragmentation()
    if isinstance(dataset, QSARTrainingset):
        if pop_outliers:
            outliers = dataset.popOutliers()
            print(f"  {len(outliers)} OUTLIERS excluded from training\n")
        else:
            dataset.reset()
    progress = Progress(len(dataset.structures()))
    for i in dataset.grindProgress(depth, speed_up):
        progress.printProgress(i)
    print("\n  Found %s molecular substructures" % len(dataset.fragments))
    print("\n     -> fragmenting time: %.0f seconds" % (time() - start))


def _makeModel(ModelConstructor, trainingset, verbose=True, **kwargs):
    if verbose: print("\n\n Extracting...\n")
    # print(f" [max_error = {kwargs['max_error']}, beta = {kwargs['beta']}]\n")
    start = time()
    model = ModelConstructor(trainingset, **kwargs)
    model['chem_library'] = chem.LIBRARY
    model.extract()
    if verbose:
        print(model.summary(verbose=True))
        print("\n     -> extracting time: %.0f seconds" % (time() - start))
    return model


def makeSARModel(trainingset: SARTrainingset, **kwargs):
    """
    Create a SAR model

    :param trainingset: the training set
    :type trainingset: SARTrainingset
    :key min_samples: minimum number of correct observation for each SA (int)
    :key max_atoms: max number of atoms for each SA (int)
    :key max_error: max admitted error for each SA (float)
    :key beta: F-score factor: recall is considered beta times as important as precision (float)
    :key explicit_Hs: allow explicit hydrogens (bool)
    :key exclusions: allow eclusions (bool)
    :key mono_target: prediction target label (str)
    :return: SAR classification model
    :rtype: SARModel
    """
    return _makeModel(SARModel, trainingset, **kwargs)


def analyzeDecomposition(model: SARModel):
    assert len(model.trainingset.labels) == 2  # BINARY DATASET
    depth = min(model.depth, 2)  # OVERWRITE depth > deeper than 2, accuracy decreases
    print("\n\n Analyzing decompositions...")
    model.extractFragments(depth)
    print(f"\n  Detected {len(model.fragments)} modulator fragments")


def makeQSARModel(trainingset: QSARTrainingset, min_samples, max_atoms, max_error, deep_learning=False, verbose=True):
    return _makeModel(QSARModel, trainingset, min_samples=min_samples, max_atoms=max_atoms,
                      max_error=max_error, deep_learning=deep_learning, verbose=verbose)


def getClustered(model: QSARModel, tolerance, verbose=True):
    print(f"\n\n Clustering modulators...")  # \t\t\t[tolerance = {tolerance}]")
    start = time()
    model = ClusterModel(model)
    model.clusterize(tolerance)
    if verbose:
        print(model.summary(verbose))
    print("\n     -> clustering time: %.0f seconds" % (time() - start))
    return model


def validateOnTrain(model):
    print("\n\n Validating on training set...")
    model.validateOnTrain()


def predict(testset, model):
    print("\n\n Predicting...")
    # start = time()
    if identities := testset.predictIdentities(model.trainingset):
        print("\n  * Detected %s IDENTITIES" % len(identities))
    if model.trainingset.excludeOutliers:
        if outliers := testset.predictOutliers(model.trainingset):
            print("\n  * Detected %s structures SIMILAR to OUTLIERS" % len(outliers))
    console_output = model.predict(testset)
    print(f"\n  {console_output}")
    # print("\n     -> predicting time: %.0f seconds" % (time() - start))


def predict_further(testset, model):
    if model['type'] == 'QSAR':
        depth = testset.fragmentation_level + 1
    else:
        depth = min(model.depth, 3)  # OVERWRITE depth, too deep and accuracy decreases
    fragmentize(testset, depth=depth, clear=False)
    # print("\n\n Fragmenting...\t[depth = {}]\n".format(depth))
    # quartiles = [100, 75, 50, 25]
    # for fragmented in testset.grindProgress(depth, progress=True):
    #     if fragmented >= quartiles[-1]:
    #         print('   {}%\tDone...'.format(quartiles.pop()))
    print("\n\n Predicting...")
    console_output = model.predictFurther(testset, depth)
    print(f"\n  {console_output}")


def saveModel(model, path):
    dumped_dataset = datasets._dumpDataset(model.trainingset)
    ruleset = {str(rule): str(rule.prediction) for rule in model.ruleset}
    parameters = model.getParameters()  # max_error is a float, all the other integers or strings
    model_stamp = {'PARAMETERS': parameters, 'RULESET': ruleset, 'TRAININGSET': dumped_dataset}
    with open(path, 'w') as f:
        try:
            json.dump(model_stamp, f)
        except Exception as e:
            raise e
        # else:
        #     model.path = path


def loadModel(path):
    with open(path, 'rb') as f:
        try:
            model_stamp = json.load(f)
        except Exception as e:
            raise e
    return model_stamp


def re_extractModel(model_stamp):
    try:
        parameters = model_stamp['PARAMETERS']
        chem.setLibrary(parameters['chem_library'])
        model_type = parameters['type']
        print("\n\n Loading dataset...\n")
        dataset = datasets._restoreDataset(model_stamp['TRAININGSET'])
        print("  Read %s molecular structures" % len(dataset))
        if model_type == 'SAR':
            trainingset = makeSARTrainingset(dataset)
        else:  # QSAR
            trainingset = makeQSARTrainingset(dataset)
        trainingset.path = parameters['train_name']
        fragmentize(trainingset, depth=parameters['depth'], speed_up=parameters['speed_up'],
                    pop_outliers=parameters.get('pop_outliers'))
        if model_type == 'SAR':
            model = makeSARModel(trainingset, **parameters)
            assert {str(rule): str(rule.prediction) for rule in model.ruleset} == model_stamp['RULESET']
            if parameters['decomposed']:
                analyzeDecomposition(model)
        else:  # QSAR
            tolerance = parameters.get('tolerance')
            model = makeQSARModel(trainingset, parameters['min_samples'], parameters['max_atoms'],
                                  parameters['max_error'], parameters['deep_learning'])
            if tolerance is not None:
                model = getClustered(model, tolerance)
            assert {str(rule): str(rule.prediction) for rule in model.ruleset} == model_stamp['RULESET']
    # except AssertionError:  # DEBUG
    #     print({mod.cansmiles: mod.bias for mod in model.modulators.values()})
    #     print(model_stamp['MODULATORS'])
    #     return None
    except Exception as e:
        raise e
    # model.path = path
    return model


def exportDB(db_path, model, info_dict):
    print("\n\n Exporting DATABASE...\n")
    start = time()
    progress = Progress(len(model.neuron_log))
    for i in delta_DB.exportDB(model, db_path, info_dict):
        progress.printProgress(i)
    print("\n {} FRAGMENTS INSERTED INTO DATABASE:\n\n '{}'".format(i, db_path))
    if delta_DB.invalid_smiles:
        print("\n WARNING: some entries have been skipped"
              "\n ({} invalid SMILES, see terminal for error log)".format(len(delta_DB.invalid_smiles)))
    print("\n     -> exporting time: %.0f seconds" % (time() - start))


# MAIN ########################################################################

# @profile
def mainLoop(trainingset, testset=None):
    start = time()

    # FRAGMENTING TRAININGSET
    if DISCRETE:
        trainingset = makeSARTrainingset(trainingset, TARGET)
        fragmentize(trainingset, depth=DEPTH, speed_up=SPEED_UP)
    else:
        trainingset = makeQSARTrainingset(trainingset)
        fragmentize(trainingset, depth=DEPTH, speed_up=SPEED_UP, pop_outliers=OUTLIERS_CHK)
        # if OUTLIERS_CHK:
        #     trainingset.popOutliers()

    # FRAGMENT TRAININGSET
    # fragmentize(trainingset, depth=DEPTH, speed_up=SPEED_UP, pop_outliers=OUTLIERS_CHK)

    # EXTRACT MODEL
    if trainingset.type == 'SAR':
        model = makeSARModel(trainingset, min_samples=MIN_SAMPLES, max_atoms=MAX_ATOMS, max_error=MAX_ERR, beta=BETA,
                             explicit_Hs=EXPLICIT_Hs, exclusions=EXCLUSIONS, mono_target=TARGET)
    elif trainingset.type == 'QSAR':
        model = makeQSARModel(trainingset, MIN_SAMPLES, MAX_ATOMS, MAX_ERR or round(trainingset.iqr / 2, 1),  # trainingset.iqr / 2  trainingset.abs_dev
                              DEEP_LEARNING)
    # OUTPUT MODEL
    Path('../EXPORTS').mkdir(exist_ok=True)
    # output.output_model(model, Path('../EXPORTS/model.csv'))  # WARNING: if SAR binary ONLY!!!
    # output.output_modAnalysis(model, Path('../EXPORTS/analysis.csv'))  # serve ancora ??

    # CLUSTER MODEL (optional)
    if not DISCRETE and GENERALIZE:
        # I could make it use the same grinder of the training set
        model = getClustered(model, BIAS_DIFF_THRES)
        # output.output_cluster_model(model, Path('../EXPORTS/cluster_model.csv'))

    # # # LOAD / SAVE MODEL
    # Path('../MODELS').mkdir(exist_ok=True)
    # saveModel(model, Path('../MODELS/modelSAR.json'))
    # model_stamp = loadModel(Path('../MODELS/modelSAR.json'))
    # model = re_extractModel(model_stamp)

    # VALIDATING
    validateOnTrain(model)

    # ANALYZE DECOMPOSITION
    if model.type == 'SAR' and len(model.trainingset.labels) == 2:
        analyzeDecomposition(model)

    # TESTING
    print("\n\n\n ***** TESTING *****")
    if testset is None:  # from file
        dataset = loadDataset(Path(TE_path), TE_targetKey, TE_IDKey, TE_smilesKey)
        testset = Testset(dataset)
    else:  # from split
        testset = Testset(testset)
    if not DISCRETE:
        # Here I fragment also potential identities...
        fragmentize(testset, depth=model.depth)  #, SPEED_UP)  # NO speed up for test set
    predict(testset, model)
    if TE_targetKey:
        model.printStatistics(testset, verbose=True)
    print("\n[ TOTAL TIME: %.0f seconds ]" % (time() - start))

    if DISCRETE and testset.unpredicted and model.fragments:
        print("\nDECOMPOSITION:")
        predict_further(testset, model)
        model.printStatistics(testset, verbose=True)

    # if testset.predictions:
    #     output.output_predictions(testset, Path('../EXPORTS/predictionsSAR.csv'), model.type)
    #     output.output_predictions(trainingset, Path('../EXPORTS/predictionsSARTrain.csv'), model.type)

    # return model


def crossvalidation(path, targetKey, IDKey, smilesKey):
    for fold in range(5):
        print(f"\n\n\n ************ FOLD {fold} ***********")
        dataset = loadDataset(Path(path), targetKey, IDKey, smilesKey)
        training, test = dataset.getSplit(fold)
        # print("split", len(training), len(test))
        yield training, test


if __name__ == "__main__":
    chem.setLibrary(LIBRARY)
    import datetime
    print(datetime.datetime.now())
    print("\n\nTRAININGSET", TR_path)
    print("CHEM LIBRARY", LIBRARY)
    print("depth", DEPTH)
    print("min_samples", MIN_SAMPLES)
    print("max_atoms", MAX_ATOMS)
    print("max_error", MAX_ERR)
    # print('SPEED_UP: %s\nOUTLIERS_CHK: %s\nGENERALIZE: %s' % (SPEED_UP, OUTLIERS_CHK, GENERALIZE))

    # # EXTERNAL TEST SET
    # trainingset = loadDataset(Path(TR_path), TR_targetKey, TR_IDKey, TR_smilesKey)
    # mainLoop(trainingset)

    # SPLIT DATASET
    # FOLD = 0
    # print("SPLIT FOLD", FOLD)
    # dataset = loadDataset(Path(TR_path), TR_targetKey, TR_IDKey, TR_smilesKey)
    # trainingset, testset = dataset.getSplit(FOLD)
    # mainLoop(trainingset, testset)

    # # CROSS-VALIDATION
    for trainingset, testset in crossvalidation(TR_path, TR_targetKey, TR_IDKey, TR_smilesKey):
        mainLoop(trainingset, testset)







    # if TE_path is None:
    #     # SPLIT 80/20
    #     TE_targetKey = TR_targetKey
    #     trainingset, testset = dataset.getSplit(2)
    #     print("split", len(trainingset), len(testset))
    # else:
    #     trainingset = dataset
    #     testset = Testset(loadDataset(Path(TE_path), TE_targetKey, TE_IDKey, TE_smilesKey))
    # # del dataset
    #
    # # REMOVE IDENTITIES
    # count = 0
    # for structure in testset.structures():
    #     if trainingset.pop(structure.cansmiles, None):
    #         count += 1
    # for structure in trainingset.structures():
    #     if testset.pop(structure.cansmiles, None):
    #         count += 1
    # print(f"REMOVED {count} identities from testset")


    # # SAVE TRAIN AND TEST
    # with open(Path('testset4_1000.csv'), 'w', encoding='utf8') as f:
    #     sep = str(',')
    #     header = sep.join(["ID", "SMILES", 'target'])
    #     f.write(header + '\n')
    #     for structure in testset.structures():
    #         row = sep.join([structure.S_id, structure.smiles, structure.target])
    #         f.write(row + '\n')
    # with open(Path('trainingset4_1000.csv'), 'w', encoding='utf8') as f:
    #     sep = str(',')
    #     header = sep.join(["ID", "SMILES", 'target'])
    #     f.write(header + '\n')
    #     for structure in trainingset.structures():
    #         row = sep.join([structure.S_id, structure.smiles, structure.target])
    #         f.write(row + '\n')

# from statistics import mean
#
# def GetRingSystems(mol, includeSpiro=False):
#     ri = mol.GetRingInfo()
#     systems = []
#     for ring in ri.AtomRings():
#         ringAts = set(ring)
#         nSystems = []
#         for system in systems:
#             nInCommon = len(ringAts.intersection(system))
#             if nInCommon and (includeSpiro or nInCommon > 1):
#                 ringAts = ringAts.union(system)
#             else:
#                 nSystems.append(system)
#         nSystems.append(ringAts)
#         systems = nSystems
#     return systems
#
# predicted = []
# unpredicted = []
#
# for predicted_smiles in testset.predictions:
#     mol = testset[predicted_smiles].mol
#     rings_num = len(GetRingSystems(mol))
#     predicted.append(rings_num)
#
# for struct in testset.unpredicted.values():
#     rings_num = len(GetRingSystems(struct.mol))
#     unpredicted.append(rings_num)
#
# print("PRED mean num of rings", mean(predicted))
# print("UNPRED mean num of rings", mean(unpredicted))


# predicted = [testset[s].mol for s in testset.predictions]
# smiles = [s for s in testset.predictions]
# pred_grid = rdkit.Draw.MolsToGridImage(predicted, molsPerRow=4, legends=smiles)
# pred_grid.save('./predictions.png')
# unpredicted = [s.mol for s in testset.unpredicted.values()]
# smiles = [s for s in testset.unpredicted]
# unpred_grid = rdkit.Draw.MolsToGridImage(unpredicted, molsPerRow=4, legends=smiles)
# unpred_grid.save('./unpredicted.png')
