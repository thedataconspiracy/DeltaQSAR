from sys import stderr
import sqlite3
from rdkit.Chem import MolFromSmiles, MolToSmiles
from parameters import DECIMAL


invalid_smiles = set()  # SMILES rising exceptions in RDKIT


def insert(table, row_dict):
    SQL_INSERT = '''INSERT INTO {table} ( {keys} ) VALUES ( {values} )'''
    for key, value in row_dict.items():
        if isinstance(value, bool):
            row_dict[key] = int(row_dict[key])
        elif isinstance(value, str):
            row_dict[key] = '"{}"'.format(row_dict[key])
        elif value is None:
            row_dict[key] = 'NULL'
    return SQL_INSERT.format(table=table,
                             keys=', '.join(row_dict.keys()),
                             values=', '.join((str(value) for value in row_dict.values()))
                             )


def convert2RDKit(smiles, chem_library):
    if chem_library == 'RDKit':
        return smiles
    try:
        return MolToSmiles(MolFromSmiles(smiles))
    except:
        if smiles not in invalid_smiles:
            invalid_smiles.add(smiles)
            stderr.write(' RDKit: {} SMILES not valid\n'.format(smiles))  # write to STDERR!


def exportDB(model, path, endpoint_info={}):
    invalid_smiles.clear()
    with sqlite3.connect(path) as con:  # using with looks faster....
        cur = con.cursor()

        # QSARpy parameters
        cur.execute('''CREATE TABLE parameters (
                        id INTEGER NOT NULL, 
                        endpoint TEXT, 
                        unit TEXT, 
                        colorbar TEXT, 
                        train_name TEXT, 
                        pop_outliers INTEGER, 
                        depth INTEGER, 
                        speed_up INTEGER, 
                        min_samples INTEGER, 
                        max_atoms INTEGER, 
                        max_error REAL, 
                        deep_learning INTEGER, 
                        type TEXT,
                        chem_library TEXT,
                        PRIMARY KEY (id)
                        )''')

        param = {}
        param.update(endpoint_info)
        param.update(model.getParameters())
        chem_library = param['chem_library']
        cur.execute(insert('parameters', param))

        # TABLE: dataset
        cur.execute('''CREATE TABLE dataset (
                        "ID" TEXT NOT NULL, 
                        cansmiles TEXT, 
                        value REAL, 
                        PRIMARY KEY ("ID")
                        )''')
        s_ids = set()
        for s in model.trainingset.structures():  # outliers already excluded
            smiles = convert2RDKit(s.cansmiles, chem_library)
            if smiles:
                cur.execute(insert('dataset', dict(ID=s.S_id, cansmiles=smiles, value=s.target)))

        # TABLE: modulators
        cur.execute('''CREATE TABLE modulators (
                        "ID" TEXT NOT NULL, 
                        cansmiles TEXT, 
                        value REAL, 
                        PRIMARY KEY ("ID")
                        )''')
        for mod in model.ruleset:
            smiles = convert2RDKit(mod.cansmiles, chem_library)
            if smiles:
                cur.execute(insert('modulators', dict(ID=mod.index, cansmiles=smiles, value=mod.bias)))

        # TABLE: observations
        cur.execute('''CREATE TABLE observations (
                        id INTEGER NOT NULL, 
                        cansmiles TEXT, 
                        "structure_IN" TEXT, 
                        "structure_OUT" TEXT, 
                        common TEXT, 
                        PRIMARY KEY (id)
                        )''')  # subtractive_mod_indexes TEXT,

        for level in range(1, model.depth):
            cur.execute('''ALTER TABLE observations ADD COLUMN mod_{} INTEGER'''.format(level))
            #  is that format safe?

        for progress, neuron in enumerate(model.neuron_log.values(), 1):
            smiles = convert2RDKit(neuron.cansmiles, chem_library)
            if not smiles:
                continue
            base_row = {'cansmiles': smiles}
            smiles_set = set()
            for network in neuron.learningNetworks:
                common_smiles = convert2RDKit(network.common_smiles, chem_library)
                if common_smiles is None:
                    continue
                row = base_row.copy()
                row['structure_IN'] = network.base.S_id
                row['structure_OUT'] = network.structure.S_id
                row['common'] = common_smiles

                smiles_set.add(network.base.cansmiles)
                smiles_set.add(network.structure.cansmiles)
                smiles_set.add(common_smiles)
                i = 1
                for n in network.neuron_dict['PLUS'].elements():
                    row['mod_' + str(i)] = n.index
                    i += 1
                    smiles_set.add(n.cansmiles)

                # DEEP LEARNING code
                # subtractive_indexes = ''
                # for n in network.neuron_dict['MINUS'].elements():
                #     row['mod_' + str(i)] = n.index
                #     subtractive_indexes += str(i)
                #     i += 1
                # row['subtractive_mod_indexes'] = subtractive_indexes

                # neuron_shift = sum(n.bias for n in network.neurons())
                # row['delta'] = round(network.structure.target - network.base.target - neuron_shift, DECIMAL)
                # # CONSISTENCY CHECK  (row['delta'] is redundant but helpful)  DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # structure_IN = db['dataset'].find_one(ID=row['structure_IN'])
                # structure_OUT = db['dataset'].find_one(ID=row['structure_OUT'])
                # delta = round(structure_OUT['value'] - structure_IN['value'] - neuron_shift, DECIMAL)
                # assert row['delta'] == round(delta, DECIMAL)

                # ids = set(row[key] for key in row if 'structure' in key or 'mod_' in key)
                if smiles_set.isdisjoint(invalid_smiles):
                    cur.execute(insert('observations', row))  # TO DO:   IGNORE if already present!
            yield progress
    con.close()


    # # BIAS ARRAY
    # deltas = np.array(frag.outputs) - np.array(frag.inputs)
    # # SAMPLE MEAN
    # mean = np.mean(deltas)
    # # STANDARD DEVIATION OF THE SAMPLE
    # std = np.std(deltas)
    # MODULATION RECORD

    # row['sample_size'] = neuron.sample_size
    #     row['mean'] = round(neuron.mean, 4)
    #     row['std'] = round(neuron.std, 4)
    #     row['record'] = json.dumps(record)
    #     try:
    #         db['fragments'].insert(row)
    #     except Exception as e:
    #         sys.exit("\n\nSQL INSERT ERROR in 'fragments' table!\n%s\n" % e.message)
    #
    # print("\n %s MOLECULAR FRAGMENTS INSERTED INTO DATABASE:  %s\n\n" % (len(fragments), dbname))
