__author__ = "Nantia Leonidou"
__description__ = " "

from cobra import *
import numpy as np


def is_member(A, B):
    """ Check which elements from A are also in B """
    check = []
    for a in A:
        if a in B:
            check.append(1)
        else:
            check.append(0)
    return check


def check_core_dead_ends(model, core_rxns):
    """ Input:
            - model
            - core_rxns: list of core reactions

        Output:
            - dead_end_core_rxns: list of dead-end core reactions"""

    # create stoichiometric matrix
    # --> S[i,j] contains the quantity of metabolite i produced (negative for consumed) by reaction j
    s_matrix = util.array.create_stoichiometric_matrix(model)

    rev = []  # get reversibility of all reactions
    for react in model.reactions:
        r_id = react.id
        rev.append(model.reactions.get_by_id(r_id).reversibility)

    dead_end_core_rxns = []
    for met in range(len(model.metabolites)):
        de_mets = 0

        both_rxns = []
        for i in range(len(s_matrix[met, :])):  # iterate over selected row (based on met) in stoichiometric matrix
            if (s_matrix[met, :][i] != 0) and (
            rev[i]):  # find all non-zero values in current row and those reactions that are reversible
                both_rxns.append(i)  # store column index, when above criterion is fulfilled

        pos_rxns = []  # produced reactions
        for i in range(len(s_matrix[met, :])):  # iterate over selected row (based on met) in stoichiometric matrix
            if s_matrix[met, :][i] > 0:
                pos_rxns.append(i)
        # union of both lists
        produced_rxns = list(set(pos_rxns) | set(both_rxns))
        # print("produced",len(produced_rxns),produced_rxns)

        neg_rxns = []  # consumed reactions
        for i in range(len(s_matrix[met, :])):  # iterate over selected row (based on met) in stoichiometric matrix
            if s_matrix[met, :][i] < 0:
                neg_rxns.append(i)
        # union of both lists
        consumed_rxns = list(set(neg_rxns) | set(both_rxns))
        # print("consumed",len(consumed_rxns), consumed_rxns)

        # Check for produced-only metabolites
        if len(consumed_rxns) == 0:
            # dead-ends: 1 if dead-ends were found, otherwise 0
            de_mets = 1
        # Check for consumed-only metabolites
        elif len(produced_rxns) == 0:
            de_mets = 1
        # Check for metabolites both consumed and produced, but only in a single
        # reversible reaction
        else:
            if (len(produced_rxns) == 1) and (len(consumed_rxns) == 1) and (produced_rxns == consumed_rxns):
                de_mets = 1
        # print(de_mets)

        # If dead end found => check for overlap with core reactions
        # If any core reaction contains a dead-end metabolite => the reaction itself will be a dead end.
        # This check avoids subsequent optimizations, as the function terminates if any blocked core reactions are detected.
        if de_mets == 1:
            # store names of found dead-end reactions
            is_de = []
            # iterate over current matrix row (current metabolite)
            for i in range(len(s_matrix[met, :].astype(np.int64))):
                if s_matrix[met, :].astype(np.int64)[i] == 1 or s_matrix[met, :].astype(np.int64)[i] == -1:
                    is_de.append(model.reactions[i].id)
            # print(is_de,len(is_de))

            # check whether detected dead-ends are core reactions as well
            if any(is_member(is_de, core_rxns)):
                dead_end_core_rxns = is_de
                break

    return dead_end_core_rxns

### test script
# model = io.mat.load_matlab_model('../../humanModel.mat')
# ### create dummy C list with core reactions
# np.random.seed(0)
# perm = np.random.permutation(len(model.reactions))[0:500]  #take only first 500 elements
# C = []
# for i in perm:
#     C.append(model.reactions[i].id)
# ### or use the following list to better test it with matlab
# C=['P4504B1r','SO4HCOtex','A_MANASEly','SPMS','ILETA','AGPRim','ESTRIOLtr','GLU5Km','UAG4Ei','UDPGLCter','CITtbm','UGALGTg',
# 'THRt4','FRUt4','PROD2m','BILGLCURtr','ARTFR55','SPHMYLNtl','PI3P3Pn','THYMDtl','MESCOALm','ABUTt2r','RTOT_2','TDCHOLAte',
# 'MI4PP','NTD5m','NRPPHRtu','EX_rbt_e','MEPIVESSte','CYSt4','NABTNO','3HPPD','ENMAN2g','AASAD3m','AGMTm','PPDOy','NACHEX9ly',
# 'EX_acn13acngalgbside_hs_e','STRDNCCPT2','6DHFtm','EX_xolest2_hs_e','EX_lac__L_e','PCm','PI34P4Pn','DM_ethamp_r','FT',
# 'FAOXC2252053m','INOSTO','TMNDNCCPT1','NTD1m','ELAIDCRNt','MI1PP','KSII_CORE4t','CSPG_Bt','ANDRSTRNtr','FUCACGALFUCGALACGLCGALGLUSIDEtg',
# 'TTDCRNt','EX_dmhptcrn_e','GALACGLCGALGBSIDEtg','COAtp','LGNCCOAtx','DIGALSGALSIDEte','GASNASEly','NACHEX6ly','PA_HSter']
# deadends = check_core_dead_ends(model, C)
# print(len(deadends), deadends)
