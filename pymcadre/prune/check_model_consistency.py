__author__ = "Nantia Leonidou"
__description__ = " Check a model's consistency by identifying inactive (dead-end) reactions"

from time import process_time
from cobra.flux_analysis import *
from prune.check_core_deadends import *
from prune.find_inactive_rxns import *


def check_model_consistency(model, method=1, rxn_to_remove=[], de_check=0, core_rxns=[]):
    """ Identification of dead-end reactions in a stoichiometric model.
        The function operates either independently to detect inactive reactions
        for a model, or within a pruning algorithm (e.g., MBA) to study the
        effect of removing reactions.
        It is mainly based on a heuristic speed-up version of the Flux Variability
        Analysis (FVA) proposed by (Jerby et al., 2010)

        Inputs:
            - model: COBRA model structure
            optional inputs:
            - method: either fastFVA (1) or fastcc (2)
            - rxn_to_remove: name of reaction to be removed (for model pruning in mCADRE or MBA)
            - de_check: check for core reactions containing dead end metabolites
                        (only for use with model pruning in MBA)
            - core_rxns: list of core reaction names (only for model pruning in MBA)

        Outputs:
            - inactive_rxns: list of IDs corresponding to reactions with 0 mininum and 0 maximum flux
            - time: CPU time required to complete function
            - result: summary indicator of dead-end effects on inactive reactions. If
                      result = 1: removal of rxn_to_remove * DID NOT *  create metabolite dead ends leading
                                  to inactivation of core reactions
                      result = 2: removal of rxn_to_remove created metabolite dead ends leading to
                                  inactivation of core reactions"""

    if len(rxn_to_remove):
        model.remove_reactions(rxn_to_remove, remove_orphans=True)

    # set all objective coefficients to 0
    for react in range(len(model.reactions)):
        model.reactions[react].objective_coefficient = 0

    inactive_rxns = rxn_to_remove

    t1_start = process_time()
    result = 1  # (until proven otherwise) assume removal or rxn_to_remove does not create any dead-end metabolite

    ###################################################################################
    # check whether any core reactions are blocked by the removal of rxn_to_remove
    ###################################################################################
    # A core reaction is considered to be dead-end , if it contains a dead-end metabolite.
    # Function exits if any blocked core reactions are detected.
    if de_check:
        dead_end_cores = check_core_dead_ends(model, core_rxns)
    else:
        dead_end_cores = []

    # If no core reactions are blocked based on metabolite dead
    # ends, maximize and minimize reactions to get those with zero flux
    # capacity.
    if len(dead_end_cores):
        # Setting inactive_rxns to include dead-end-containing reactions will
        # effectively cause the function to exit without checking non-core
        # reactions below; thus, the full list of inactive reactions will not
        # be enumerated
        inactive_rxns = list(set(dead_end_cores) | set(inactive_rxns))

        # indicates that dead-end-containing reactions were found in the core
        result = 2

    else:
        # Otherwise, fastFVA is used to quickly scan through all
        # reactions. **note: may want to include option to use fastFVA with GLPK
        fast_scan = find_inactive_rxns(model, method)
        inactive_rxns = list(set(inactive_rxns) | set(fast_scan))

    t1_stop = process_time()
    exec_time = t1_stop - t1_start
    print('Execution time of check_model_consistency: %s s' % exec_time)
    print('check_model_consistency done ... ')

    return inactive_rxns, exec_time, result
