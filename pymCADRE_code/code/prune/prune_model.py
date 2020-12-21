__author__ = "Nantia Leonidou"
__description__ = " Prune Model "

from check.check_model_function import *
from prune.check_salvage_path import *
from prune.check_model_consistency import *
import numpy as np


def prune_model(GM, P, C, Z, eta, precursor_mets, salvage_check=1, method=1, cutoff=np.inf):
    # default salvage_check value is 1 --> assume non-hepatic tissue, where nucleotides are not expected to be produced de novo

    # Initialize variables
    R_G = [r.id for r in GM.reactions]
    PM = GM
    R_P = R_G

    NC_removed = 0
    C_removed = 0
    cRes = [0] * 3000
    count = 1

    while len(P) != 0 and count < cutoff:
        print('Reaction no. ', count)
        r = P[0]
        print('Attempting to remove reaction: ', r, '...')
        print('Initial PM reactions:', len(PM.reactions))
        modelR = PM.copy()
        modelR.remove_reactions([r], remove_orphans=True)

        # First check precursor production;
        # if this test fails, no need to check model consistency with FVA (time-saving step)
        rStatus = check_model_function(modelR, 'required_mets', precursor_mets)[0]

        # If specified, check the salvage pathway as well
        if salvage_check:
            rSalvage = check_salvage_path(modelR)[0]
            rStatus = rStatus and rSalvage

        if rStatus:
            # Check for inactive reactions after removal of r using a copy of PM
            copyPM = PM.copy()
            inactive_G = sorted(check_model_consistency(copyPM, method, [r])[0]) # get inactive reactions, includes a remove command

            inactive_C = list(set(inactive_G) & set(C))
            inactive_NC = list(set(inactive_G) - set(inactive_C))
            print('Lengths all--> inactive_G: %s , inactive_C: %s, inactive_NC:%s' %(len(inactive_G), inactive_C, inactive_NC))

            # Remove reactions with zero expression (previously penalized in
            # rank_reactions) and corresponding inactive core reactions, only if
            # sufficiently more non-core reactions are removed
            if r in Z:
                print('Zero-expression evidence for reaction...')
                # Check model function with all inactive reactions removed
                modelTmp = PM.copy()
                modelTmp.remove_reactions(inactive_G,remove_orphans=True)
                tmpStatus = check_model_function(modelTmp, 'required_mets', precursor_mets)[0]

                # If specified, check the salvage pathway as well
                if salvage_check:
                    tmpSalvage = check_salvage_path(modelTmp)[0]
                    tmpStatus = tmpStatus and tmpSalvage

                if (len(inactive_C) / len(inactive_NC) <= eta) and tmpStatus:
                    R_P_ids = [item.id for item in GM.reactions]
                    R_P = sorted(list(set(R_P_ids) - set(inactive_G)))

                    PM.remove_reactions(inactive_G,remove_orphans=True)
                  
                    for i in inactive_G:
                        if i in P:
                            P.remove(i)
                    NC_removed = NC_removed + len(inactive_NC)
                    C_removed = C_removed + len(inactive_C)
                    num_removed = NC_removed + C_removed
                    print('Removed all inactive reactions')

                    # result = -1.x --> reaction r had zero expression evidence and was removed along with any
                    # consequently inactivated reactions;
                    # x indicates the number of core reactions removed
                    if len(inactive_C) > 100:
                        removed_C_indicator = len(inactive_C) / 100
                    else:
                        removed_C_indicator = len(inactive_C) / 10
                    result = -1 - removed_C_indicator

                else:
                    # no reactions (core or otherwise) are actually
                    # removed in this step, but it is necessary to update the
                    # total number of removed reactions to avoid errors below
                    num_removed = NC_removed + C_removed
                    P.pop(0)
                    print('No reactions removed')

                    # result = 1.x --> no reactions were removed because removal of r either led to a ratio
                    # of inactivated core vs. non-core reactions above the specified threshold
                    # eta (x = 1) or the removal of r and consequently inactivated reactions prevented
                    # production of required metabolites (x = 0)
                    result = 1 + tmpStatus / 10

            # If reaction has expression evidence, only attempt to remove
            # inactive non-core reactions
            else:
                print('Reaction %s is not in Z' %(r))
                # when r is not in Z
                # Check model function with non-core inactive reactions removed
                modelTmp = PM.copy()
                modelTmp.remove_reactions(inactive_NC, remove_orphans=True)
                tmpStatus = check_model_function(modelTmp, 'required_mets', precursor_mets)[0]

                # If specified, check the salvage pathway as well
                if salvage_check:
                    tmpSalvage = check_salvage_path(modelTmp)[0]
                    tmpStatus = tmpStatus and tmpSalvage

                if (len(inactive_C) == 0) and tmpStatus:
                    R_P_ids = [item.id for item in GM.reactions]
                    R_P = sorted(list(set(R_P_ids) - set(inactive_NC)))
                    PM.remove_reactions(inactive_NC, remove_orphans=True)
                    for i in inactive_G:
                        if i in P:
                            P.remove(i)
                    NC_removed = NC_removed + len(inactive_NC)
                    num_removed = NC_removed + C_removed
                    print('Removed non-core inactive reactions')

                    # result = -2.x --> reaction r had expression
                    # evidence and was removed along with (only) non-core inactivated reactions;
                    # x indicates the number of core reactions removed (should be zero!)
                    if len(inactive_C) > 100:
                        removed_C_indicator = len(inactive_C) / 100
                    else:
                        removed_C_indicator = len(inactive_C) / 10

                    result = -2 - removed_C_indicator

                else:
                    num_removed = NC_removed + C_removed
                    P.pop(0)
                    print('No reactions removed')

                    # result = 2.x --> no reactions were removed because removal of r either led to inactivated core
                    # reactions (x = 1) or the removal of r and consequently inactivated reactions prevented production
                    # of required metabolites (x = 0)
                    result = 2 + tmpStatus / 10

        else:
            num_removed = NC_removed + C_removed
            P.pop(0)

            # result = 3 --> no reactions were removed because
            # removal of r by itself prevented production of required metabolites
            result = 3

        cRes[count] = result
        count += 1
        print('Num. removed: %s (%s core, %s non-core); Num. remaining: %s' % (num_removed, C_removed, NC_removed, len(P)))

    cRes[count:] = []
    print('prune_model done .... ')
    return PM, cRes


