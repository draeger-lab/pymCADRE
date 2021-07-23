__author__ = "Nantia Leonidou"
__description__ = " Calculate expression-based evidence "

# import numpy as np
# from cobra import *
# from rank.parse_gprs import *
# from rank.map_high_conf_to_rxns import *
# from rank.map_gene_scores_to_rxns import *
# from test_Inputs.test_Inputs import *
import math
import numpy as np


def calc_expr_evidence(model, GPR_rxns, U_GPR, is_C_H):
    """ Compute expression-bsed evidence: function of how often genes, associated to a reaction,
        is expressed in the selected context. E_X is ranged from 0 to 1 and its relation to
        the ubiquity scores is a composite of the Boolean gene-reaction rules defined in the generic
        model.

        Input:
            - generic model --> cobra.io.core.model.Model
            - GPR_rxns --> list of reaction IDs, obtained from GPRs, as strings
            - U_GPR --> list of lists, with ubiquity scores as float numbers
            - is_C_H --> list of 0 or 1 indicating high confidence or not

        Output:
            - E_X --> list of respective expression-based evidence. If 0, then this reaction
                      has strong evidence of not being active in the respective tissue context"""

    # create list with zeros
    E_X = [0.0] * len(model.reactions)

    # get minimum of each list in U_GPR(list of lists)
    U_GPR_min = [np.nanmin(l) for l in U_GPR]

    all_rxns = [r.id for r in model.reactions]

    # find existence of model.reactions in the GPR_rxns list
    for i in range(len(all_rxns)):
        rxn_GPRS = []  # store indices, where a model reaction is found
        if all_rxns[i] in GPR_rxns:
            indices = [index for index, value in enumerate(GPR_rxns) if value == all_rxns[i]]
            rxn_GPRS = indices

        # if rxns_GPRS not empty
        if len(rxn_GPRS) != 0:
            tmp_lst = [U_GPR_min[k] for k in rxn_GPRS]
            E_X[i] = np.nanmax(tmp_lst)

    # for reactions with no corresponding probe in expression data
    for i in range(len(is_C_H)):
        if is_C_H[i] == 1:
            E_X[i] = 1.0

    return E_X
