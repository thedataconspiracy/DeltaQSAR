#  DeltaQSAR, a knowledge discovery toolkit for automatic SAR/QSAR rules induction.
#  Copyright (C) 2021  Thomas Ferrari

from parameters import DECIMAL

SGN_FRMT = '{' + ':+.{}f'.format(DECIMAL) + '}'
FLT_FRMT = '{' + ':.{}f'.format(DECIMAL) + '}'
STR_FRMT = '"{}"'


def string(s):
    return STR_FRMT.format(s)


def output_predictions(dataset, path, type='SAR'):
    predictions = dataset.predictions

    with open(path, 'w', encoding='utf8') as f:
        sep = str(',')
        header = sep.join(["ID", "SMILES", "PREDICTION", "TARGET", "ERROR", "Predicted by (means of)",
                           "Predicted by (SMILES)", "Predicted by (ID)", "Prediction sequence"])

        f.write(header + '\n')
        for structure in dataset.values():
            pred_by = smiles_cell = id_cell = seq = error = prediction = ''
            target = str(structure.target or '')
            if structure.cansmiles in predictions:
                record = predictions[structure.cansmiles]
                predictor = record.get('predicted by')
                prediction = record.get('prediction')
                if structure.target is not None:
                    error = FLT_FRMT.format(record.get('error'))
                if record['mode'] == 'identity':
                    pred_by = "Identity"
                    smiles_cell = str(predictor)
                    id_cell = predictor.getID()
                elif record['mode'] == 'similar to outliers':  # QSAR only
                    pred_by = "Similarity"
                    outliers = sorted(predictor, key=lambda s: s.getID())
                    smiles_cell = '  '.join((repr(s) for s in outliers))
                    id_cell = '  '.join((s.getID() for s in outliers))
                    seq = '  '.join((str(s.target) for s in outliers))
                elif record['mode'] == 'Same fragments':
                    pred_by = "Permutation"
                    same_frags_structs = sorted(predictor, key=lambda s: s.getID())
                    smiles_cell = '  '.join((repr(s) for s in same_frags_structs))
                    id_cell = '  '.join((s.getID() for s in same_frags_structs))
                    # if type == 'QSAR':
                    seq = '  '.join((str(s.target) for s in same_frags_structs))
                    # smiles_cell = str(predictor)
                    # id_cell = predictor.S_id
                elif record['mode'] == 'SA':
                    pred_by = 'Structural Alert'
                    smiles_cell = str(predictor)
                    id_cell = predictor.getID()
                else:
                    pred_by = "Decomposition"  # predictor.type()  # "Factorization"
                    network = predictor
                    smiles_cell = str(network.base)
                    id_cell = network.base.getID()
                    seq = str(network.base.target)
                    i = 1
                    for mod in network.neurons():
                        if mod in network.neuron_dict['PLUS']:
                            op = "  + "
                        else:
                            op = "  - "
                        smiles_cell += op + str(mod)
                        id_cell += op + mod.getID()
                        if type == "QSAR":
                            seq += " > " + FLT_FRMT.format(network.output_log[i])
                        i += 1
                    if type == "SAR":
                        seq += ' > ' + str(network.output_log) + ' > ' + network.prediction

            row = sep.join([string(structure.getID()), string(structure), str(prediction), target, error,
                            pred_by, string(smiles_cell), string(id_cell), string(seq)])
            f.write(row + '\n')

            ## GHOST MODULATORS
            ##if neuron.aggregates:
            ##last_activation = structure.network.outputs[i-1]
            ##for agg in neuron.aggregates: #.sort(key=attrgetter('cansmiles')):
            ##neurons += op + str(agg)
            ### trick temporaneo per rendere invisibili i ghost modulators nell'output
            ##activation = last_activation + agg.bias
            ##if rev: activation = last_activation - agg.bias
            ##actseq += " > "+FORMAT % activation
            ##last_activation = activation
            ###if "%.2f" % activation != "%.2f" % structure.network.outputs[i]:  # DEBUG
            ###pass


def output_model(model, path):
    sep = str(',')
    with open(path, 'w', encoding='utf8') as f:
        if model.type == "SAR":
            # STRUCTURAL ALERTS
            f.write(sep.join(['ID', 'SMARTS', 'Prediction', 'precision', f'F({model.beta})-score']) + '\n')
            for rule in model.ruleset:
                if model.beta == 0:
                    F_score = FLT_FRMT.format(rule.precision)
                else:
                    F_score = FLT_FRMT.format(rule.F_score)
                f.write(sep.join([string(rule.getID()), string(rule), rule.prediction, FLT_FRMT.format(rule.precision), F_score]) + '\n')
            if model.fragments:  # DECOMPOSITION MODULATORS
                l1, l2 = model.trainingset.labels
                f.write('\n\nBinary dataset: decomposition modulators\n\n')
                f.write(sep.join(['ID', 'SMILES', l1 + ' dir', l1 + ' rev', l2 + ' dir', l2 + ' rev']) + '\n')
                for rule in model.fragments:
                    f.write(sep.join([string(rule.getID()), string(rule),
                                      string(rule.bias_matrix[l1]['direct']),
                                      string(rule.bias_matrix[l1]['reverse']),
                                      string(rule.bias_matrix[l2]['direct']),
                                      string(rule.bias_matrix[l2]['reverse'])]) + '\n')
        else:
            f.write(sep.join(['ID', 'SMILES', 'endpoint SHIFT', 'abs ERROR']) + '\n')
            for rule in model.ruleset:
                f.write(sep.join([string(rule.getID()), string(rule), SGN_FRMT.format(rule.bias), str(rule.error)]) + '\n')
        # f.write('\n\n')
        # _appendTrainingset(model, f)
        # MODEL PARAMETERS
        f.write('\n\nMODEL PARAMETERS:\n\n')
        for key, value in model.getParameters().items():
            f.write(sep.join([key, STR_FRMT.format(value)]) + '\n')


def output_cluster_model(model, path):
    # CLUSTER MODULATOR HAS NO "GROUP ERROR"..............................
    dict_view = model.clustered_view
    sep = str(',')
    headers = ['Group', 'variant', 'ID', 'SMILES', 'endpoint SHIFT', 'abs ERROR']
    with open(path, 'w', encoding='utf8') as f:
        f.write(sep.join(headers) + '\n')
        for label in dict_view.labels(sort=True):
            cluster = dict_view[label]
            # first = True
            for i, variant in enumerate(cluster):
                if len(cluster) == 1:
                    GL = ''
                else:
                    # GL = string(chr(945 + i))  # variants denoted by greek letters
                    GL = string(i + 1)  # variants denoted by integers
                shift = variant[0].bias
                error = variant[0].error  # GENERALIZED ERROR !!!!!!!!!!!!!!!!!!!!!!!!
                smiles = '  |  '.join((repr(mod) for mod in variant))
                id = '  |  '.join((mod.getID() for mod in variant))
                # if first:
                #     first = False
                # else:
                if i > 0:
                    label = ''
                f.write(sep.join((string(label), GL, string(id), string(smiles), SGN_FRMT.format(shift), str(error))) + '\n')
                # for mod in variant:
                #     assert mod.bias == shift  # DEBUG!
                #     assert mod.error == error # DEBUG!
        if dict_view['PLACEHOLDERS']:
            label = 'DUMMY'
            desc = '(shift = 0)'
            for ph in dict_view['PLACEHOLDERS']:
                f.write(sep.join((label, desc, string(ph.getID()), string(ph), SGN_FRMT.format(ph.bias), str(ph.error))) + '\n')
                desc = ''
        # if dict_view['COMBO']:
        #     f.write('\nREDUNDANCIES:' + '\n\n')
        #     f.write(sep.join(headers + ['', "Combination of:"]) + '\n\n')
        #     for combo in dict_view['COMBO']:
        #         label = prefix = ''
        #         for i, agg in enumerate(combo.aggregates):
        #             label += SGN_FORMAT.format(agg.bias) + ' [ %s ]' % agg
        #             # if i != 0:
        #             #     prefix = ' + '
        #             # label += prefix + FORMAT % agg.bias + ' [ ' + str(agg) + ' ]'
        #         f.write(sep.join([combo.getID(), str(combo), SGN_FORMAT.format(combo.original_bias),
        #                           SGN_FORMAT.format(combo.bias), '', label]) + '\n')
        # f.write('\n\n')
        # _appendTrainingset(model, f)
        # MODEL PARAMETERS
        f.write('\n\nMODEL PARAMETERS:\n\n')
        for key, value in model.getParameters().items():
            f.write(sep.join([key, STR_FRMT.format(value)]) + '\n')


# def _appendTrainingset(model, f):
#     sep = str(',')
#     # headers = ["ID", "endpoint value", '', '', '', "SMILES"]
#     headers = ["ID", "SMILES", "endpoint value"]
#     f.write('BASE STRUCTURES:\n\n')
#     f.write(sep.join(headers) + '\n\n')
#     for base in model.trainingset.structures():
#         f.write(sep.join([base.S_id, str(base), str(base.target)]) + '\n')
#     outliers = model.trainingset._outliers
#     if outliers:
#         f.write('\n')
#         f.write('OUTLIERS:\n\n')
#         # f.write(sep.join(headers) + '\n\n')
#         # for out in sorted(outliers.values(), key=attrgetter('position')):
#         for out in outliers.structures():
#             f.write(sep.join([out.S_id, str(out), str(out.target)]) + '\n')


# # MODULATORS ANALYSIS
# def output_modAnalysis(model, path):
#     if model.type == 'SAR':
#         rules = model.fragments
#     else:
#         rules = model.ruleset
#     sep = str(',')
#     header = sep.join(["ID", "SMILES", "TRAINING ref. (ID)", "Structural factorization (IDs)",
#                        "Structural factorization (SMILES)", "Equation", "INPUT", "OUTPUT",
#                        "X (OUTPUT - INPUT)", "FINAL endpoint SHIFT"])
#     with open(path, 'w', encoding='utf8') as f:
#         f.write(header + '\n')
#         for rule in rules:
#             for i in range(rule.sample_size):
#                 network = rule.learningNetworks[i]
#                 networkID = network.base.getID()
#                 networkSmiles = str(network.base)
#                 equation = str(network.base.target)
#                 for neuron in network.neurons():
#                     # if neuron.aggregates:
#                     # neuronID = " + ".join([rule.reprID() for rule in neuron.aggregates])
#                     # neuronSMI = " + ".join([str(agg) for agg in neuron.aggregates])
#                     # neuronBIAS = " ".join([SGN_FORMAT % agg.old_bias for agg in neuron.aggregates])
#                     # else:
#                     # neuronID = neuron.reprID()
#                     # neuronSMI = str(neuron)
#                     # neuronBIAS = SGN_FORMAT % neuron.old_bias
#                     networkID += " + %s" % neuron.getID()
#                     networkSmiles += " + %s" % neuron
#                     equation += " " + SGN_FRMT.format(neuron.bias)
#                 networkID += " + [ITSELF] = %s" % network.structure.getID()
#                 networkSmiles += " + [ITSELF] = %s" % network.structure
#                 equation += " + X = " + str(network.structure.target)
#                 row = sep.join((string(rule.getID()), string(rule), string(network.structure.getID()),
#                                string(networkID), string(networkSmiles), string(equation),
#                                str(rule.inputs[i]), str(rule.outputs[i]),
#                                SGN_FRMT.format(rule.outputs[i] - rule.inputs[i]), str(rule.bias)))
#                 f.write(row + '\n')


# def output_pretty_modulators(model, filename):
# dict_view = model.generalized
# with open(filename, 'w') as f:
# sep = ','
# f.write('MODULATORS:\n\n')
# f.write(sep.join(["ID", "LABEL", "SMILES", "endpoint SHIFT\n(before generalization)",\
# "endpoint SHIFT"]) + '\n')
# for label in dict_view.labels(sort=True):
# f.write(label + '\n')
# for cluster in dict_view[label]:
# print cluster
# many_smiles = getCleanString([str(mod) for mod in cluster])
# shift = cluster[0].bias
# f.write(sep.join(['', many_smiles, SGN_FORMAT % shift]) + '\n')
# f.write('\n')
# if dict_view['PLACEHOLDERS']:
# f.write('Placeholders' + '\n')
# for ph in dict_view['PLACEHOLDERS']:
# f.write(sep.join(['', str(ph), SGN_FORMAT % ph.bias]) + '\n')
# f.write('\n\n')
##f.write('MODULATORS (operational):\n\n')
##f.write(sep.join(["ID", "SMILES", "endpoint SHIFT"]) + '\n')
##for mod in sorted(modulators, key=attrgetter('id')):
##f.write(sep.join([mod.reprID(), mod.smiles, SGN_FORMAT % mod.bias]) + '\n')
# f.write('\n\n')
# f.write('BASE STRUCTURES:\n\n')
# f.write(sep.join(["ID", "SMILES", "endpoint value"]) + '\n')
# for base in model.trainingset.bases.getSorted():
# f.write(sep.join([str(base.id), base.smiles, str(base.target)]) + '\n')
# f.write('\n\n')
# f.write('OUTLIERS:\n\n')
# f.write(sep.join(["ID", "SMILES", "endpoint value"]) + '\n')
# for out in sorted(model.trainingset.outliers.values(), key=attrgetter('position')):
# f.write(sep.join([str(out.id), out.smiles, str(out.target)]) + '\n')