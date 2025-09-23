# LIBRARY = 'OB'
LIBRARY = 'RDKit'
MAX_DEC_FAST = 500
MAX_DEC = 5000
# print("WARNING: MAX_DEC = 50")

DECIMAL = 2  # round

# DISCRETE = False
DISCRETE = True   # <- comment out if QSAR


# FRAGMENTING
SPEED_UP = True
# SPEED_UP = False
DEPTH = 4

# EXTRACTION
MAX_ATOMS = 18
MIN_SAMPLES = 4


if DISCRETE:  # SAR

    TARGET = None  # If TARGET == None: multiclass classification
    MAX_ERR = 0.25
    BETA = 0.05  # F-score factor (recall is considered BETA times as important as precision)
    EXPLICIT_Hs = False
    EXPLICIT_Hs = True
    EXCLUSIONS = False
    # EXCLUSIONS = True

    # # TRAINING SET
    # TR_path = "../DATASETS/Muta_500.sdf"
    # TR_IDKey = 'CAS_RN'
    # TR_targetKey = 'Ames_test_categorisation'
    # TR_smilesKey = None
    # # TARGET = 'mutagen'

    # TR_path = "../DATASETS/dataset_MUTA_KNN.csv"
    TR_path = "../DATASETS/dataset_MUTA_KNN_1000.csv"
    # TR_path = "../DATASETS/dataset_MUTA_KNN_200.csv"
    TR_IDKey = 'CAS'
    TR_targetKey = 'Experimental value'
    TR_smilesKey = 'SMILES'
    # TARGET = 'Mutagenic'

    # TR_path = "../DATASETS/18K_debug.csv"
    # TR_IDKey = None
    # TR_targetKey = 'exp'
    # # TARGET = 'mutagen'

    # TEST SET
    TE_path = "../DATASETS/765_KNN_2.csv"
    TE_targetKey = 'experimental'
    TE_IDKey = 'Id'
    TE_smilesKey = 'SMILES'

    # TE_path = "../DATASETS/Muta_500_testknn.sdf"  #None
    # TE_targetKey = 'Ames_test_categorisation'
    # TE_IDKey = 'CAS_RN'
    # TE_smilesKey = None  # 'SMILES'



else:  # QSAR

    MAX_ERR = None  # admitted error threshold: default (None) = absolute deviation
    OUTLIERS_CHK = False
    OUTLIERS_CHK = True
    DEEP_LEARNING = False
    # DEEP_LEARNING = True

    TR_path = "../DATASETS/QSAR_DEMO_TRAIN.csv"
    # TR_path = "../DATASETS/TRAIN_0.csv"
    TR_IDKey = 'CAS_mod'
    TR_targetKey = 'Experimental'
    TR_smilesKey = 'SMILES'

    TE_path = "../DATASETS/QSAR_DEMO_TEST.csv"
    # TE_path = "../DATASETS/TEST_0.csv"
    TE_IDKey = 'CAS_mod'
    TE_targetKey = 'Experimental'
    TE_smilesKey = 'SMILES'

    GENERALIZE = False  # Enable/disable the Generalization tool  (values: True/False)
    BIAS_DIFF_THRES = 0.3  # HIGHLY OPTIMIZED for VEGA dataset -> [0.3]
    # BIAS_DIFF_THRES -> value of the lower significant error (domain-dependent)

    # Related modulators can be grouped into a "cluster".
    # A cluster is a collection of items which "distance" is below a certain threshold.
    # "Distance" is defined as the BIAS difference.
    #
    # There are 3 ways of clusterizing dependig on the "relation" between modulators:
    #
    # REDUNDANCY: When a modulator structure is the union of the structures of other modulators
    #  (M = M1+M2), then generalize M with M1 and M2
    #  (i.e., remove modulator M and use modulators M1 and M2 instead)
    # CLUSTERIZATION: when more modulators share the same structure,
    #  then generalize them by a single meta-modulator that matches with all SMILES,
    #  the final bias is the weighted average
    #  (i.e., group C* and C=* in a meta-modulator "C" that matches C* and C=*)
    # PLACEHOLDERS (only for modulator with a very small bias):
    #  if |bias| < BIAS_DIFF/2  then  bias = 0
    #  all these modulators don't really affect the output so they are removed,
    #  they are just kept as placeholders without any effect
    #  (since QSARpy predict a structure only when all of its fragments are known)
    #
    # TIP:
    #  Keep it lower than the bias of C* if you don't want to ignore all C atoms!


class Progress:
    quartiles = [100, 75, 50, 25]

    def __init__(self, total):
        self.total = total

    def printProgress(self, i):
        for j, perc in enumerate(self.quartiles):
            if (i + 1) * 100 / self.total >= perc:
                print(f'   {perc}%\tDone...')
                self.quartiles = self.quartiles[:j]
                break
