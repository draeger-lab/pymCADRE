__author__ = "Nantia Leonidou"
__description__ = " Map high confidence genes to reactions. "

import numpy as np
# from rank.parse_gprs import *
# from test_Inputs.test_Inputs import *


def map_high_conf_to_rxns(model, GPR_file, GPR_rxns, C_H_genes=[]):
    """ Map high confidence genes to reactions.

        Input:
            - model --> cobra.io.core.model.Model
            - GPR_file --> list of lists, containing gene IDs as strings
            - GPR_rxns --> list of reaction IDs, obtained from GPRs, as strings
            - C_H_genes --> list filled with Entrez IDs (as strings) for any genes with particularly strong evidence for activity in the selected tissue
                            Within the rank_reactions module, the reactions corresponding to these genes will automatically become
                            assigned to the "core" set to avoid removal.

        Output:
            - is_C_H --> list of 0 or 1 indicating high confidence or not """

    if len(C_H_genes) != 0:
        # get all lengths of all sub-lists inside the GPR_file:
        l = [len(GPR_file[i]) for i in range(len(GPR_file))]
        # get the maximum and use it as number of columns, that C_H_GPR should maximally have
        columns = max(l)
        # check whether the input gene IDs from the C_H_genes file can be found in the GPR_file:
        # create a list of same length as GPR_file, containing 1 if the
        # ID from GPR_file can be found in C_H_genes list, 0 when it is not
        # found and NaN otherwise
        GPR_file = [[s.split('_')[0] for s in lst] for lst in GPR_file] # remove numbers after comma to match format of Entrez IDs
        # # --> split('_') for R3con3D, split('.') for humanModel
        new_GPR_file = [[r.strip() for r in lst] for lst in GPR_file] # remove unnecessary space before and after strings

        C_H_GPR = [[np.nan for x in range(columns)] for y in range(len(new_GPR_file))]
        for i in range(len(new_GPR_file)):
            for j in range(len(new_GPR_file[i])):
                if new_GPR_file[i][j] in C_H_genes:
                    C_H_GPR[i][j] = 1
                else:
                    C_H_GPR[i][j] = 0

        # get minimum value of each row (0 or 1)
        C_H_GPR_min = [min(lst) for lst in C_H_GPR]

        # Is confidence high?
        is_C_H = [0] * len(model.reactions)

        for idx,rxn in enumerate(model.reactions):
            # compare current reaction ID to each row of the GPR_rxns, looking for an exact match of the entire character vector
            rxn_GPRs = []
            for k in range(len(GPR_rxns)):
                if rxn.id == GPR_rxns[k]:
                    rxn_GPRs.append(k)  # if a match is found, store related indices from GPR_rxns file
            if len(rxn_GPRs) != 0:
                # based on reactions' indices, get corresponding values from the C_H_GPR_min
                to_find_max = [C_H_GPR_min[index] for index in rxn_GPRs]
                # store whether it is of high confidence
                is_C_H[idx] = max(to_find_max)
        # is_C_H = logical(is_C_H);

    else:
        # if C_H_genes file is not provided, fill the list with zeros
        is_C_H = [0] * len(model.reactions)

    return is_C_H

### test script
# model = io.mat.load_matlab_model('../../humanModel.mat')
# C_H_genes = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[2]
# G = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[0]
# U = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[1]
# GPR_file = parse_gprs(model)[1]
# GPR_rxns = parse_gprs(model)[0]
# print(map_high_conf_to_rxns(model, GPR_file, GPR_rxns, C_H_genes))
