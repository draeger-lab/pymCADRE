__author__ = "Nantia Leonidou"
__description__ = " Detect inactive reactions using either fastFVA oder FASTCC "

from cobra.flux_analysis import flux_variability_analysis, fastcc
from cobra.flux_analysis import find_blocked_reactions
# from cobra import *
# from cobra.io.sbml import *


def find_inactive_rxns(model, method=1):
    """ ** Runnable only in jupyter notebook **
        Find all inactive (blocked) reactions in a cobra model using either fastFVA or FASTCC methods

        Input:
            - model: cobra.io.core.model.Model
            - method: fastFVA(1, default) or FASTCC (2)

        Output:
            - inactive_rxns: list of inactive reactions """

    # inactive_rxns = []
    if method == 1:  # fastFVA (by default, if no method number is given)

        print("Checking all reactions (fastFVA)...")
        # set all objective coefficients to 0
        for react in range(len(model.reactions)):
            model.reactions[react].objective_coefficient = 0

        inactive_rxns = find_blocked_reactions(model)
        #####################################################################################
        #      **       ALTERNATIVELY:          **
        # following code delivers the same inactive reactions as find_blocked_reactions
        #####################################################################################
        # conduct FVA analysis
        #fva = flux_variability_analysis(model) # add number of allowed tasks (threads) to run at the same time
        #print("FVA is done...")
        # store all maximum and all minimum flux-values after FVA
        #opt_max = []
        #opt_min = []
        #for i in range(len(fva)):
        #    opt_max.append(fva.iloc[i][1])
        #    opt_min.append(fva.iloc[i][0])
        #print("storing done...")
        # set constraints to consider a reaction "inactive"
        #abs_max = [i for i in range(len(opt_max)) if abs(opt_max[i]) < 1e-6]
        #abs_min = [i for i in range(len(opt_min)) if abs(opt_min[i]) < 1e-6]
        #print("Constraints done..")
        #is_inactive = [rxn_id for rxn_id in abs_max if rxn_id in abs_min]
        # get all (inactive) blocked reactions' id
        #inactive_rxns = []
        #for rxn_id in is_inactive:
        #    inactive_rxns.append(model.reactions[rxn_id].id)
        #########################################################################

    else:  # otherwise use FASTCC

        print("Checking all reactions (FASTCC)...")
        # get non-blocked reactions, * in MATLAB the fastcc function gives less is_active reactions *
        is_active = fastcc(model, zero_cutoff=1e-4).reactions
        # get all active reactions' id
        active_rxns = [rxn.id for rxn in is_active]
        # get all model's reactions
        all_model_rxns = [rxns.id for rxns in model.reactions]
        # get difference, which are then the inactive reactions
        inactive_rxns = list(set(all_model_rxns) - set(active_rxns))

    print('Model consists of %s blocked (inactive) reactions' % len(inactive_rxns))
    return inactive_rxns

