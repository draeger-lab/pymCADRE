__author__ = "Nantia Leonidou"
__description__ = " Rank reactions "

from rank.parse_gprs import *
from rank.map_high_conf_to_rxns import *
from rank.map_gene_scores_to_rxns import *
from rank.calc_expr_evidence import *
from rank.initialize_generic_model import *
from rank.calc_conn_evidence import *
import math
import pandas as pd


def rank_reactions(model, G, U, confidence_scores, C_H_genes, method=1):
    """
        Input:
            - model (cobra.core.model.Model)
            - G: gene IDs from expression data (lst)
            - U: gene ubiquity scores (i.e., from mas5callToExpression) (lst)
            - C_H_genes: high confidence genes (lst, optional)
            - method: FVA(1, default) or FASTCC (2) (int)

        Output:
            - GM: generic model with inactive reactions removed (cobra.core.model.Model)
            - C:  core reactions (lst)
            - NC: non-core reactions (lst)
            - P:  removal order of non-core reactions (lst) """

    ###############################################
    # Parse GPRs
    ###############################################
    GPRrxns, GPRmat = parse_gprs(model)

    ###############################################
    # Map high confidence genes to reactions
    ###############################################

    if len(C_H_genes) != 0:
        is_C_H = map_high_conf_to_rxns(model, GPRmat, GPRrxns, C_H_genes)
    else:
        is_C_H = [0] * len(model.reactions)

    ###############################################
    # Map gene ubiqiuty scores to reactions
    ###############################################
    U_GPR = map_gene_scores_to_rxns(model, G, U, GPRmat)

    ###############################################
    # Determine confidence level-based evidence
    ###############################################
    E_L = confidence_scores

    # Reactions with no confidence information should not be ranked higher than
    # those with non-zero confidence
    for idx, elem in enumerate(E_L):
        if math.isnan(elem):
            E_L[idx] = 0

    ###############################################
    # Calculate expression-based evidence
    ###############################################
    E_X = calc_expr_evidence(model, GPRrxns, U_GPR, is_C_H)
    C = []
    for i in range(len(E_X)):
        if E_X[i] > 0.9:
            C.append(model.reactions[i].id)
    model_C = C

    ######################################################################
    # Initialize the consistent generic model & update evidence vectors
    ######################################################################
    GM, C, E_X, E_L = initialize_generic_model(model, C, E_X, confidence_scores, method)
    R_G = GM.reactions
    # get IDs
    R_G_ids = [item.id for item in R_G]

    NC = []
    NC_idx = []
    for rg_idx, i in enumerate(R_G_ids):
        if i not in C and i not in NC:
            NC_idx.append(R_G.index(i))
            NC.append(i)

    ################################################
    # Calculate connectivity-based evidence
    ################################################
    E_C = calc_conn_evidence(GM, E_X)  # very slow due to matrix multiplication (2766x3742)*(3742x2766)

    ###############################################
    # Rank non-core reactions
    ###############################################
    E_X_NC = [E_X[i] for i in NC_idx]  # expression-based evidence for non-core reactions
    E_C_NC = [E_C[i] for i in NC_idx]  # connectivity-based evidence for non-core reactions
    E_L_NC = [E_L[i] for i in NC_idx]  # literature - based evidence for non - core reactions
    # put all evidences in a singe list
    data_tuples = list(zip(E_X_NC, E_C_NC, E_L_NC))
    # create data frame
    df_evidence = pd.DataFrame(data_tuples, columns=['E_X_NC', 'E_C_NC', 'E_L_NC'])
    # sort it
    E_NC = df_evidence.sort_values(by=['E_X_NC', 'E_C_NC', 'E_L_NC'])
    # get indices as a list
    NC_order = list(E_NC.index)

    # ordered (ranked) non-core reactions
    P = [NC[i] for i in NC_order]

    # convert data frame to numpy array, for easier analysis
    E_NC = E_NC.T.to_numpy()

    ###############################################
    # Identify zero-expression reactions
    ###############################################
    Z = []
    for i in range(len(E_NC[0])):
        # print(i)
        if E_NC[0][i] == -1e-6:
            # print(E_NC['E_X_NC'][i])
            Z.append(P[i])

    print('rank_reactions done ... ')
    return GM, C, NC, P, Z, model_C
