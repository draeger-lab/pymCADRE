__author__ = "Nantia Leonidou"
__description__ = " Create a consistent generic model by removing all inactive reactions"

# import numpy as np
# from cobra import *
# from rank.map_gene_scores_to_rxns import *
# from rank.parse_gprs import *
# from rank.calc_expr_evidence import *
# from cobra.flux_analysis import *
# from rank.map_high_conf_to_rxns import *
from prune.check_model_consistency import *
import math


def initialize_generic_model(model, C, E_X, confidence_scores, method=1):
    """ This function creates a consistent generic model by removing all inactive
        reactions. It will also return adjusted vectors for expression-based and
        literature-based evidence, corresponding to the subset of reactions in
        generic model GM.

        Input:
            - model --> cobra.io.core.model.Model
            - C --> list of Cobra core reactions
            - E_X --> list of respective expression-based evidence
            - confidence_scores --> experimental-based scores
            - method --> fastFVA (1) or FASTCC (0)

        Output:
            - GM --> generic model, includes only core reactions
            - C --> list of Cobra core reactions
            - E_X --> list of expression-based evidence
            - E_L --> list of confidence level-based evidence"""
    # define literature-based evidence from confidence_scores
    if len(confidence_scores):
        E_L = confidence_scores

        # reactions with no confidence information should not be ranked higher than
        # those with non-zero confidence
        for idx, elem in enumerate(E_L):
            if math.isnan(elem):
                E_L[idx] = 0
    else:
        # create an empty list with zeros
        E_L = [0.0] * len(model.reactions)

    # get inactive reactions
    inactive_rxns = check_model_consistency(model, method)[0]

    # update core set of reactions
    C = list(set(C) - set(inactive_rxns))

    # create generic model by removing inactive reactions
    GM = model.copy()
    GM.remove_reactions(inactive_rxns, remove_orphans=True) # remove_orphans--> to remove corresponding METABOLITES and GENES as well!
    #                                                                           (in MATLAB only metabolites are removed, but that is not a problem,
    #                                                                            as genes are not used further in the code )

    # update list of expression-based evidence
    GM_rxns_ids = [rxn.id for rxn in GM.reactions]  # get reaction IDs from GM model
    model_rxn_ids = [rxn.id for rxn in model.reactions]  # get reaction IDs from initial model

    GM_idx = []
    model_idx = []
    for idx, i in enumerate(GM_rxns_ids):
        for g_idx, j in enumerate(model_rxn_ids):
            if i == j and g_idx not in model_idx:
                model_idx.append(g_idx)
                GM_idx.append(idx)

    E_X_GM = [0] * len(GM.reactions)
    for i in range(len(GM_idx)):
        E_X_GM[i] = E_X[model_idx[i]]
    E_X = E_X_GM

    # update list of confidence level-based evidence
    E_L_GM = [0] * len(GM.reactions)
    for j in range(len(GM_idx)):
        E_L_GM[j] = E_L[model_idx[j]]
    E_L = E_L_GM

    return GM, C, E_X, E_L



### test script
# model = io.mat.load_matlab_model('../../humanModel.mat')
# C_H_genes = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[2]
# G = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[0]
# U = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[1]
# confidence_scores = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[3]
# GPR_rxns = parse_gprs(model)[0]
# GPR_file = parse_gprs(model)[1]
# is_C_H = map_high_conf_to_rxns(model, GPR_file, GPR_rxns, C_H_genes)
# U_GPR = map_gene_scores_to_rxns(model, G, U, GPR_file)
# E_X = calc_expr_evidence(model, GPR_rxns, U_GPR,is_C_H)
# C = []
# for i in range(len(E_X)):
#     if E_X[i] > 0.9:
#         C.append(model.reactions[i].id)
# GM, C, E_X, E_L = initialize_generic_model(model,C,E_X,confidence_scores,1)
# print(len(GM.reactions))
# print(len(C))
# print(len(E_X))
# print(len(E_L))