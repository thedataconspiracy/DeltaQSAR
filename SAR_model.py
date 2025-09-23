#  DeltaQSAR, a knowledge discovery toolkit for automatic SAR/QSAR rules induction.
#  Copyright (C) 2021  Thomas Ferrari


from StructuralAlert import StructuralAlert
from datasets import Testset
from SAR_model_decomposer import extractFragments, predictFragments, predictSameFragmets, networksForSAR
from parameters import DECIMAL, Progress


class SARModel(dict):

    def __init__(self, trainingset, **kwargs):
        """ kwargs: min_samples, max_error, max_atoms, mono_target, beta, explicit_Hs, exclusions """
        super().__init__(kwargs)
        self['train_name'] = trainingset.path
        self['speed_up'] = trainingset.fast_fragmentation
        self['depth'] = trainingset.fragmentation_level
        self['decomposed'] = False
        self['type'] = 'SAR'
        self.__dict__ = self.copy()  # SET keys AS ATTRIBUTES
        self.trainingset = trainingset
        self.ruleset = None
        self.min_precision = 1 - self.max_error
        self.fragments = None  # QSARpy decomposition

        # MONKEY-PATCH to add QSARpy decomposition technology to SARModel
        if len(trainingset.labels) == 2:  # BINARY CLASSIFICATION
            networksForSAR(trainingset.labels[0], trainingset.labels[1])

    def __repr__(self):
        if self.isMonoTarget():
            summary = f"{len(self.ruleset)} Structural Alerts ({self['mono_target']})"
        else:
            summary = f"{len(self.ruleset)} Structural Alerts:"
            counter = dict.fromkeys(self.targets(), 0)
            for SA in self.ruleset:
                counter[SA.prediction] += 1
            details = " {} ({}) /" * len(counter)
            details = details.rstrip('/')
            args = [elem for pair in ((k, v) for k, v in counter.items()) for elem in pair]
            summary += details.format(*args)
        if self.fragments:
            summary += f" + {len(self.fragments)} Modulators"
        return summary

    def summary(self, verbose=False):
        return f"\n  MODEL: {self}"

    def getParameters(self):
        return dict(self)

    def isMonoTarget(self):
        return self.mono_target is not None

    def targets(self):  # FIX THIS
        if self.isMonoTarget():
            return {self.mono_target}
        return self.trainingset.labels

    def validateOnTrain(self):
        testset = Testset(self.trainingset)
        testset.labels = self.trainingset.labels
        self.predict(testset)
        # if testset.unpredicted and self.isBinary():
        #     self.predictFurther(testset)
        self.printStatistics(testset)
        self.trainingset.predictions = testset.predictions  # to output training predictions
        assert testset.path == self.trainingset.path  # DEBUG

    def extract(self, verbose=True):
        if self.isMonoTarget():
            targets = len(self.trainingset.subset[self.mono_target])
            try:
                LR = (self.min_precision / self.max_error) * (targets / (len(self.trainingset) - targets))  # likelihood ratio
            except ZeroDivisionError:
                precision_thresholds = {self.mono_target: 1}
            else:
                precision_thresholds = {self.mono_target: LR / (1 + LR)}  # adjusted precision to consider prevalence
        else:
            precision_thresholds = {target: len(self.trainingset.subset[target]) / len(self.trainingset)\
                                    for target in self.trainingset.labels}

        alerts = self.getAlerts(precision_thresholds, verbose)
        low_fence = min(precision_thresholds.values())
        # print("POTENTIAL ALERTS", len(alerts))  # DEBUG
        self.ruleset = []
        for idx, rule in enumerate(self.RuleGenerator(alerts, low_fence), 1):
            rule.index = idx
            self.ruleset.append(rule)
            if verbose:
                if idx == 1:
                    print("\n\n Structural Alerts (SA):", end='')
                    if self.isMonoTarget():
                        print(f"\t[{self.mono_target} ONLY]", end='')
                    print('\n\n  ' + StructuralAlert.SUMMARY_HEADERS)
                print(rule.summary())
        # rule.interferenceCheck(ruleset[i + 1:])  # INTERFERENCE CHECK !

        # atom_distr = dict.fromkeys(range(1, self.max_atoms + 1), 0)
        # depth_distr = dict.fromkeys(range(1, self.depth + 1), 0)
        # for r in self.ruleset:
        #     atom_distr[r.fragment.atoms - r.fragment.wildcards] += 1
        #     depth_distr[r.fragment.wildcards] += 1
        # print("ATOMS:", atom_distr)
        # print("DEPTH:", depth_distr)

    def extractFragments(self, depth):
        self.fragments = extractFragments(self, depth)
        if self.fragments:
            self['decomposed'] = True

    def predict(self, dataset):
        c = 0
        for SA in self.ruleset:
            for structure in SA.getMatches(dataset.structures()):
                dataset.setPrediction(structure.cansmiles, SA)
                c += 1
        return f"Found {c} molecular structures with Structural Alerts"

    def predictFurther(self, dataset, depth):
        previous_predictions = len(dataset.predictions)
        predictSameFragmets(self.trainingset, dataset, depth, self.max_error, self.mono_target)
        predictFragments(self, dataset, depth)
        return f"Predicted {len(dataset.predictions) - previous_predictions} molecular structures by structural decomposition"

    def getAlerts(self, precision_thresholds, verbose):
        alerts = {}
        # c = 0
        if verbose: progress = Progress(len(self.trainingset.fragments))
        for i, fragment in enumerate(self.trainingset.fragments):
            if verbose: progress.printProgress(i)
            if fragment.atoms - fragment.wildcards > self.max_atoms:
                continue

            if fragment.cansmiles_unstarred not in alerts:  # prova a toglierlo e controlla sia stabile
                SA = StructuralAlert(fragment, self.trainingset.labels, nostar=True)
                SA.match(self.trainingset.structures(), prescreen_FP=self.speed_up)
                if self._activateAlert(SA, precision_thresholds):
                    alerts[SA.smarts_string] = SA

    # DEBUG
    #         if not SA.train_hits:  # or 'H' in SA.smarts.string:
    #             c += 1
    #             print(SA.smarts.string, SA.fragment)
    #             for struct in self.trainingset.getBySub(SA.fragment.cansmiles):
    #                 print(struct.smiles, struct.cansmiles)
    # print(c, "BAD SMARTS")

            if self.explicit_Hs:
                base_SA = alerts.get(fragment.cansmiles_unstarred, SA)
                SA_H = StructuralAlert(fragment, self.trainingset.labels, addHs=True)
                SA_H.match(base_SA.train_hits, prescreen_FP=self.speed_up)
                if self._activateAlert(SA_H, precision_thresholds):
                    alerts[SA_H.smarts_string] = SA_H

        return list(alerts.values())

    def _activateAlert(self, alert, precision_thresholds):
        for target in self.targets():
            alert.calcMetrics(self.trainingset, self.beta, target)
            if self._checkAlert(alert, precision_thresholds[target], target):
                alert.activate(target)
                return True
            elif self.isMonoTarget() and self.exclusions:  # MONKEY-PATCH EXCEPTIONS MONO
                if alert.countFalse(target) >= self.min_samples:
                    alert.activate('exception')
                    return True

    def _checkAlert(self, alert, precision_threshold, target=None):
        return alert.countTrue(target) >= self.min_samples and alert.precision > precision_threshold  # NO >= (multiclass)

    def RuleGenerator(self, alerts, threshold):
        if self.exclusions: self.initExclusions(alerts)
        current_dataset = self.trainingset.copy()

        while alerts:  # and current_dataset:  no need if check other_alert.countTrue()
            alerts.sort()
            if self.isMonoTarget() and self.exclusions:
                # Even if mono-target, I get alerts for all classes (as exclusions), but I make them unreachable
                alerts.sort(key=lambda alert: alert.prediction == self.mono_target)

            SA = alerts.pop()
            if SA.precision < self.min_precision:
                continue
            if self.isMonoTarget() and SA.prediction != self.mono_target:  # No more alerts for class mono_target
                break

            # Apply the selected SA to the dataset
            for cansmiles, target in SA.yieldHits():
                current_dataset.pop(cansmiles)
            # Update and prune other alerts
            popped = {SA}
            for i in range(len(alerts) - 1, -1, -1):
                other_alert = alerts[i]
                other_alert.removeHits(SA)
                if other_alert.prediction == 'exception':  # MONKEY-PATCH EXCEPTIONS MONO
                    if other_alert.countFalse() < self.min_samples:
                        popped.add(alerts.pop(i))
                    continue
                other_alert.resetExclusions()
                other_alert.calcMetrics(current_dataset, self.beta)
                if not self._checkAlert(other_alert, threshold):
                    popped.add(alerts.pop(i))
                    continue
            # Here discarded alerts are popped also from potential_exclusions: but a bad alert is not a bad exclusion...
            if self.exclusions: self.updateExclusions(alerts, current_dataset, popped)
            yield SA

    def initExclusions(self, alerts):
        alerts = sorted(alerts, key=lambda a: a.size())
        for i, alert in enumerate(alerts):
            if alert.prediction == 'exception':  # MONKEY-PATCH EXCEPTIONS MONO
                continue
            alert.setPotentialExcusions(alerts[i + 1:])  # assign superstructures
            if True:
            # if alert.precision >= self.min_precision:
                alert.assignExclusions(self.min_samples, self.beta)
                alert.calcMetrics(self.trainingset, self.beta)
                # TEST: ASSERT CHE E` MIGLIORATO

    def updateExclusions(self, alerts, dataset, discard_set=None):
        for i, alert in enumerate(alerts):
            if alert.prediction == 'exception':  # MONKEY-PATCH EXCEPTIONS MONO
                continue
            alert.potential_exclusions -= discard_set
            if True:
            # if alert.precision >= self.min_precision:
                alert.assignExclusions(self.min_samples, self.beta)  # after updating ALL alerts (so that superSA updated also), compute metrics
                alert.calcMetrics(dataset, self.beta)

    def printStatistics(self, dataset, verbose=False):
        if not dataset.predictions:
            print("\n *** Unpredicted DATASET")
            return
        dataset.validate()
        labels = {structure.target for structure in dataset.structures(all=True)}
        if not labels.issuperset(self.targets()):  # self.targets() <= labels:
            print("\n *** Can't compute metrics: MODEL and DATASET endpoint labels are incompatible!")
            # print(" There is a '%s' prediction, but the DATASET doesn't have such key..." % targets - dataset.labels)
            return

        def printMAE(description, errors, total=None, MAE=True):
            if errors:
                samples = len(errors)
                sum_abs_errors = sum((abs(error) for error in errors))
                if total:
                    err_desc = "  MAE (Mean Absolute Error) = "
                    perc = " ({0:.1f}%)".format(samples * 100 / total)
                else:
                    err_desc = "   MAE = "
                    perc = ""
                print(description + str(samples) + perc)
                if MAE:
                    print(err_desc, round(sum_abs_errors / samples, DECIMAL))

        for level in range(1):  # ,4):
            # print("\nLEVEL:", level)
            # print("GENERIC")
            error_dict = {}
            for record in dataset.predictions.values():
                # if record['mode'] != 'Same fragments' and len(record['predicted by'].neuron_dict['PLUS']) + len(record['predicted by'].neuron_dict['MINUS']) == level:
                error_dict.setdefault(record['mode'], []).append(record['error'])
            if verbose:  # DEBUG
                printMAE("\n  Predicted by Structural Alerts: ", error_dict.get('SA'))
                printMAE("\n  Predicted by same fragments: ", error_dict.get('Same fragments'))
                printMAE("\n  Predicted by direct factorization: ", error_dict.get('direct'))
                printMAE("\n  Predicted by reverse factorization: ", error_dict.get('reverse'))
                printMAE("\n  Predicted by substitution: ", error_dict.get('generic'))
            else:
                printMAE("\n  Predicted by factorization: ", error_dict.get('direct', [])
                         + error_dict.get('reverse', []) + error_dict.get('generic', []))

        all_errors = [error for sublist in error_dict.values() for error in sublist]

        if self.isMonoTarget():
            print(f"\n  Binary classification:  {{'{self['mono_target']}': POSITIVE, otherwise: NEGATIVE}}")
            binaryOutput(dataset, self['mono_target'])
        else:
            printMAE("\n  TOTAL predicted: ", all_errors, len(dataset), MAE=False)
            multiclassOutput(dataset, labels)


def binaryOutput(dataset, target):
    TP = TN = FP = FN = 0
    for record in dataset.predictions.values():
        if record['target'] == target and record['prediction'] == target:
            TP += 1
        elif record['target'] != target and record['prediction'] != target:
            TN += 1
        elif record['target'] == target and record['prediction'] != target:
            FN += 1
        elif record['target'] != target and record['prediction'] == target:
            FP += 1
    for structure in dataset.unpredicted.structures():
        if structure.target != target:
            TN += 1
        else:
            FN += 1
    accuracy = (TP + TN) / (TP + FN + TN + FP)
    targets = TP + FN
    nontargets = TN + FP
    if targets == 0:
        sensitivity = 0
    else:
        sensitivity = TP / targets
    if nontargets == 0:
        specificity = 0
    else:
        specificity = TN / nontargets
    print("\n  ACCURACY:\t{:.{d}f}".format(accuracy, d=DECIMAL))
    print("  sensitivity:\t{:.{d}f}".format(sensitivity, d=DECIMAL))
    print("  specificity:\t{:.{d}f}".format(specificity, d=DECIMAL))
    tab = 9
    print(f"\n\n  CONFUSION MATRIX:\n")
    print("\t".expandtabs(4), end='')
    print("Positive\tNegative\t \u2190 PREDICTIONS".expandtabs(tab))
    print("\t".expandtabs(4), end='')
    print(f"  {TP}\t  {FN}\t   Positive".expandtabs(tab))
    print("\t".expandtabs(4), end='')
    print(f"  {FP}\t  {TN}\t   Negative".expandtabs(tab))


def multiclassOutput(dataset, labels):
    labels = sorted(labels)  # for reproducibility
    outcomes = labels + ['unpredicted']
    matrix = {label: {outcome: 0 for outcome in outcomes} for label in labels}
    # populate confusion matrix
    for cansmiles in dataset:
        if cansmiles in dataset.predictions:
            record = dataset.predictions[cansmiles]
            matrix[record['target']][record['prediction']] += 1
        else:
            structure = dataset[cansmiles]
            matrix[structure.target]['unpredicted'] += 1

    total = len(dataset)
    hits = sum(matrix[label][label] for label in labels)
    unpredicted = total - len(dataset.predictions)
    errors = len(dataset.predictions) - hits
    assert unpredicted == sum(matrix[label]['unpredicted'] for label in labels)  # DEBUG
    assert hits + errors + unpredicted == total  # DEBUG
    print(f"  ACCURACY: {hits / len(dataset.predictions):.{DECIMAL}f}")
    # print("  ERRORS: {perc:.{d}f}%".format(perc=100 * errors / len(dataset.predictions), d=1))
    # print("  UNPREDICTED: {perc:.{d}f}%".format(perc=100 * unpredicted / total, d=1))
    print("\n\n  CONFUSION MATRIX:\n")
    print("\t".expandtabs(4), end='')
    tab = 9
    for outcome in outcomes:
        print(f"{outcome[:tab - 1]}\t".expandtabs(tab), end='')
    print(" \u2190 PREDICTIONS")
    for label in labels:
        print("\t".expandtabs(4), end='')
        for outcome in outcomes:
            print(f"  {matrix[label][outcome]}\t".expandtabs(tab), end='')
        print("   " + label)


        # # NON-Mutagenic
        # print("\nFROM NON MUTA")
        # nm_error_dict = {}
        # for record in dataset.predictions.values():
        #     if record['predicted by'].base.target == 'NON-Mutagenic':
        #         nm_error_dict.setdefault(record['mode'], []).append(record['error'])
        # if True:
        #     printMAE("\n  Predicted by direct factorization: ", nm_error_dict.get('direct'))
        #     printMAE("\n  Predicted by reverse factorization: ", nm_error_dict.get('reverse'))
        #     printMAE("\n  Predicted by substitution: ", nm_error_dict.get('generic'))
        #
        # print("\nFROM MUTA")
        # muta_error_dict = {}
        # for record in dataset.predictions.values():
        #     if record['predicted by'].base.target == 'Mutagenic':
        #         muta_error_dict.setdefault(record['mode'], []).append(record['error'])
        # if True:
        #     printMAE("\n  Predicted by direct factorization: ", muta_error_dict.get('direct'))
        #     printMAE("\n  Predicted by reverse factorization: ", muta_error_dict.get('reverse'))
        #     printMAE("\n  Predicted by substitution: ", muta_error_dict.get('generic'))
        #
        # # SWITCHING
        # print("\nby switching")
        # sw_error_dict = {}
        # for record in dataset.predictions.values():
        #     if record['predicted by'].base.target != record['prediction']:
        #         sw_error_dict.setdefault(record['mode'], []).append(record['error'])
        # if True:
        #     printMAE("\n  Predicted by direct factorization: ", sw_error_dict.get('direct'))
        #     printMAE("\n  Predicted by reverse factorization: ", sw_error_dict.get('reverse'))
        #     printMAE("\n  Predicted by substitution: ", sw_error_dict.get('generic'))
        #
        # print("\nby switching from MUTA")
        # swm_error_dict = {}
        # for record in dataset.predictions.values():
        #     if record['predicted by'].base.target != record['prediction'] and record['predicted by'].base.target == 'Mutagenic':
        #         swm_error_dict.setdefault(record['mode'], []).append(record['error'])
        # if True:
        #     printMAE("\n  Predicted by direct factorization: ", swm_error_dict.get('direct'))
        #     printMAE("\n  Predicted by reverse factorization: ", swm_error_dict.get('reverse'))
        #     printMAE("\n  Predicted by substitution: ", swm_error_dict.get('generic'))
        #
        # print("\nby switching from NON MUTA")
        # swnm_error_dict = {}
        # for record in dataset.predictions.values():
        #     if record['predicted by'].base.target != record['prediction'] and record['predicted by'].base.target == 'NON-Mutagenic':
        #         swnm_error_dict.setdefault(record['mode'], []).append(record['error'])
        # if True:
        #     printMAE("\n  Predicted by direct factorization: ", swnm_error_dict.get('direct'))
        #     printMAE("\n  Predicted by reverse factorization: ", swnm_error_dict.get('reverse'))
        #     printMAE("\n  Predicted by substitution: ", swnm_error_dict.get('generic'))
        #
        # # NOT SWITCHING
        # print("\nby NOT switching")
        # nsw_error_dict = {}
        # for record in dataset.predictions.values():
        #     if record['predicted by'].base.target == record['prediction']:
        #         nsw_error_dict.setdefault(record['mode'], []).append(record['error'])
        # if True:
        #     printMAE("\n  Predicted by direct factorization: ", nsw_error_dict.get('direct'))
        #     printMAE("\n  Predicted by reverse factorization: ", nsw_error_dict.get('reverse'))
        #     printMAE("\n  Predicted by substitution: ", nsw_error_dict.get('generic'))