import networks as net
import SAR_networks


def networksForSAR(label_0, label_1):
    net.Network = SAR_networks.SARNetwork
    net.LearningNetwork = SAR_networks.SARLearningNetwork
    SAR_networks.SARNetwork.setBinaryLabels(label_0, label_1)


def extractFragments(model, depth):
    model.neuron_log = {}  # populated by _getNeurons (to keep neurons unique)

    modes = ('direct', 'reverse')  # SAR_networks.SARNeuron.learning_modes
    labels = model.trainingset.labels
    model.neurons = {label: {mode: {} for mode in modes} for label in labels}
    # self.neurons is a dict of dicts, e` ridondante: (se un neurone e stato attivato per diverse label/mode, compare piu` volte)

    verbose = False
    for mode in modes:
        # if verbose: print("\n\nMODE: ", mode)
        for label in labels:  # learning from this label
            # if verbose: print("\n", label)
            knowledge = model.neurons[label][mode]
            for level in range(1, depth + 1):
                # if verbose: print("\nLivello", level)
                trained = _getNeurons(model, level, knowledge, label, mode, verbose)
                selected = _selectModulators(trained, label, mode, model.min_samples)
                knowledge.update(selected)
                _interference(model.trainingset, knowledge, verbose=False)  # OKKIO: NON FA NULLA!
    # 'REVERSE' mode: il livello 1 puo essere dedotto dalle dirette, e ai livelli successivi si puo` dedurre i samples concordi
    # print("neuron_log", len(model.neuron_log))

    modulator_set = set((mod for mode in modes for label in labels for mod in model.neurons[label][mode].values()))
    # mod.sample_size e` implementato come somma di tutti i samples, lo uso solo qui....................................
    fragments = sorted(modulator_set, key=lambda mod: (mod.sample_size, mod.cansmiles), reverse=True)
    if verbose: print("\n\n RULES\n")
    for i, mod in enumerate(fragments):
        mod.index = i + 1
    # outputNeurons(self.ruleset)
    return fragments


def _getNeurons(model, level, knowledge, label, mode, verbose):
    if mode == 'reverse':
        sub_dataset = model.trainingset
        super_dataset = model.trainingset.subset[label]  # .copy()
    else:
        sub_dataset = model.trainingset.subset[label]  # .copy()
        super_dataset = model.trainingset
    # Modulatori INERTI rispetto alla stesso target (label)
    # Qui scarto il learning da modulatori con precisione sotto la soglia, ma ci sarebbero info da estrarre...
    # (se tutti i campioni sono concordi, che il modulatore precedente sia preciso o meno mi dice poco)
    modulators = {k: v for k, v in knowledge.items()
                  if (v.bias_matrix[label][mode] is not None) and
                  (abs(v.bias_matrix[label][mode]) <= model.max_error)}

    trained = set()
    # if verbose:
    #     print("super", len(super_dataset), '\n', "sub", len(sub_dataset), '\n', "knowledge", knowledge, '\n', modulators", modulators)
    #     print("BASES:", label)
    for structure in super_dataset.structures():
        # propagating max_atoms to the generators saves time, but it also excludes some obs from self.neuron_log...
        # AT PRESENT: the check is implemented for deep networks only
        for network in net.learningNetworks(structure, sub_dataset, modulators, level, reverse=(mode == 'reverse')):
            # if verbose: print("\nNetwork", network)
            neuron = network.learner
            if str(neuron) in knowledge:
                # NOTA: anche se in learningNetworks() sono esclusi i learner gia` approvati, questo e` vero solo
                # per gli INERTI (solo loro sono dentro modulators), qua escludo anche quelli in knowledge
                # if verbose: print(" REMOVED - Learner from previous level", str(neuron), "\n\n")
                continue
            if network.switch_probability > model.max_error:  # INERT network ONLY
                # if verbose: print(" REMOVED - Bad network", network.error, "\n\n")
                continue
            network.learner = model.neuron_log.setdefault(str(neuron), neuron)  # keep neurons UNIQUE (check and swap)
            network.trainLearner()
            if network.learner.fragment.atoms <= model.max_atoms:
                trained.add(network.learner)
            # print(f"LEVEL {level}, LABEL {label}, MODE {mode}: {len(trained)} trained")
    # if verbose: print(f"\nSELECTING for label {label}\ntotal neurons", len(trained))
    return trained


def _selectModulators(neurons, label, mode, min_samples):
    modulators = {}
    for neuron in neurons:
        if len(neuron.sample_matrix[label][mode]) >= min_samples:
            neuron.activate(label, mode)
            modulators[neuron.cansmiles] = neuron
    return modulators


def _interference(trainingset, neurons, verbose):
    # NON FA NULLA!
    if verbose: print("\nINTERFERENCE check:")
    l = sorted(neurons.values(), key=lambda x: x.fragment.atoms, reverse=True)
    for i, n1 in enumerate(l):
        for n2 in l[i + 1:]:
            n1_hits = set(s.cansmiles for s in trainingset.getBySub(n1.cansmiles))
            n2_hits = set(s.cansmiles for s in trainingset.getBySub(n2.cansmiles))
            # if str(n1) == '*[N+](=O)[O-]' and str(n2) == '*[N+](=*)[O-]':
            #     for h in n1_hits:
            #         if h not in n2_hits:
            #             print("...", h)
            if n1_hits <= n2_hits:  # AND n1.samples == n2.samples.....
                if verbose: print("INTERFERENCE", n1, n2)
                # print("removing", n2)
                # try:
                #     del neurons[n2.cansmiles]
                # except:
                #     print("already removed?", n2)
        if verbose: print(n1)


def outputNeurons(ruleset):
    nmd = md = nmr = mr = 0
    for i, mod in enumerate(ruleset):
        print(mod, mod.bias_matrix)
        print(" ", mod.sample_matrix)
        if mod.bias_matrix['Mutagenic']['direct'] is not None: md += 1
        if mod.bias_matrix['Mutagenic']['reverse'] is not None: mr += 1
        if mod.bias_matrix['NON-Mutagenic']['direct'] is not None: nmd += 1
        if mod.bias_matrix['NON-Mutagenic']['reverse'] is not None: nmr += 1
        # rev = len(mod.sample_matrix['Mutagenic']['reverse']) + len(mod.sample_matrix['NON-Mutagenic']['reverse'])
        # dir = len(mod.sample_matrix['Mutagenic']['direct']) + len(mod.sample_matrix['NON-Mutagenic']['direct'])
        # assert dir == rev
            # print("INCONGRU dir", dir, "rev", rev)
    print("NON direct", nmd)
    print("NON reverse", nmr)
    print("mut direct", md)
    print("mut reverse", mr)


#  ############ PREDICTION ##################

def predictSameFragmets(trainingset, testset, depth, max_error, mono_target=None):
    for test_structure in list(testset.structures()):
        matches = {label: [] for label in trainingset.labels}
        for train_structure in trainingset.structures():
            # if not (test_structure is train_structure):  # if validate on training set
            if mono_target is not None and train_structure.target != mono_target:
                continue
            if test_structure.checkSameFragments(train_structure, depth):
                matches[train_structure.target].append(train_structure)
        match_count = len([match for subset in matches.values() for match in subset])
        if match_count > 0:
            for target in matches.keys():
                if 1 - len(matches[target]) / match_count <= max_error:  # check in case of discordant predictions
                    testset._setPrediction(test_structure.cansmiles, "Same fragments", matches[target], target)
                    break


def predictFragments(model, testset, depth):
    labels = model.trainingset.labels
    bases = model.trainingset.subset
    unpredicted = testset.unpredicted
    # onTrain = bool(testset.path == self.trainingset.path)
    level = 1

    while unpredicted and level <= depth:  # x-networks can predict deeper

        network_dict = {}  # { cansmiles: [networks] }
        for label in labels:
            # Direct Networks
            for structure in unpredicted.structures():
                network_dict.setdefault(structure.cansmiles, []) \
                    .extend(net.networkGen(structure, bases[label], model.neurons[label]['direct'], level))
            # Reverse Networks
            if depth <= model.depth:  # disable reverse prediction if depth > trainingset.fragmentation_level
                for train_structure in model.trainingset.subset[label].structures():
                    for network in net.networkGen(train_structure, unpredicted, model.neurons[label]['reverse'], level,
                                                  reverse=True):
                        network_dict.setdefault(network.structure.cansmiles, []).append(network)

        # Prediction
        for cansmiles, networks in network_dict.items():
            if networks:
                best_network = min(networks)
                if model.isMonoTarget() and best_network.prediction != model.mono_target:
                    continue
                if best_network.error <= model.max_error:
                    testset.setPrediction(cansmiles, best_network)

        # if level > 1:  # IF UNPREDICTED: predict by substitution
        #     for label in labels:
        #         neurons = dict(model.neurons[label]['reverse'])
        #         neurons.update(model.neurons[label]['direct'])
        #         for structure in list(unpredicted.structures()):
        #             networks = []
        #             for te_depth, tr_depth in net.weak_compositions(level):
        #                 networks.extend(net.deepNetworks(structure, bases[label], neurons, (te_depth, tr_depth), onTrain))
        #             if networks:
        #                 best_network = min(networks)
        #                 if model.isMonoTarget() and best_network.prediction != model.mono_target:
        #                     continue
        #                 if best_network.error <= model.max_error:
        #                     testset.setPrediction(structure.cansmiles, best_network)

        level += 1


# TEST DIRECT AND REVERSE

# PRINT TO SCREEN
# pred = best_network.prediction()
# assert min(pred, 1 - pred) <= 1 - self.min_precision  # DEBUG
# print()
# print(best_network.type(), 'level', level)
# print(testset.get(cansmiles).target, cansmiles)
# print("PRED:", best_network.prediction(), best_network)
# label = best_network.base.target
# mode = best_network.type()
# print("BASE:", label)
# print([n.bias_matrix[label][mode] for n in best_network.neurons()])
# print("BIAS", best_network.error)
# if len(networks) > 1:
#     print("\n\tDISCARDED NETWORKS")
# for net in networks[1:]:
#     print()
#     print('\t', net.type(), 'level', level)
#     print('\t', testset.get(cansmiles).target, cansmiles)
#     print('\t', "PRED:", net.prediction(), net)
#     label = net.base.target
#     mode = net.type()
#     print('\t', [n.bias_matrix[label][mode] for n in net.neurons()])
#     print('\t', "UPDATE", net.error)


# TEST SUBSTITUTION

# pred = best_network.prediction()
# print()
# print("***DEEP level", level)
# print(best_network.type())
# print(f"{testset.get(structure.cansmiles).target} predicted {best_network.prediction()}", structure.cansmiles)
# print("NETWORK:", best_network)
# label = best_network.base.target
# mode = best_network.type()
# print("BASE:", label)
# if mode in ('reverse', 'direct'):
#     print("\tNetwork sequence:", [n.bias_matrix[label][mode] for n in best_network.neurons()])
# else:
#     print("\tNetwork sequence:", [n.bias_matrix[label]['reverse'] for n in best_network.neurons('MINUS')],
#           [n.bias_matrix[label]['direct'] for n in best_network.neurons('PLUS')])
# print("BIAS", best_network.error)
# if len(networks) > 1:
#     print("\n\tDISCARDED NETWORKS")
#     for net in networks[1:]:
#         print()
#         print('\t', net.type(), 'level', level)
#         print('\t', testset.get(cansmiles).target, cansmiles)
#         print('\t', "PRED:", net.prediction(), net)
#         label = net.base.target
#         mode = net.type()
#         if mode in ('reverse', 'direct'):
#             print("\tNetwork sequence:",
#                   [n.bias_matrix[label][mode] for n in best_network.neurons()])
#         else:
#             print("\tNetwork sequence:", [n.bias_matrix[label]['reverse'] for n in
#                                         best_network.neurons('MINUS')],
#                   [n.bias_matrix[label]['direct'] for n in
#                    best_network.neurons('PLUS')])
#         print('\t', "UPDATE", net.error)