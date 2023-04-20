__author__ = "Nantia Leonidou"
__description__ = " Modify metabolites' bounds "

from check.find_organic_ex_rxns import *


def set_organic_met_bounds(model, ex_rxns):
    """ Set lower bound of organic reactions to 0.0, aiming the turn-off of the uptake

        Input: a cobra model and a list of reaction IDs

        Output: a cobra model with the updated bounds  """

    # Identify all exchange reactions that include organic metabolites
    organic_ex_rxns = find_organic_ex_rxns(model, ex_rxns)

    # Reset all lower bounds for organic reactions to 0 to disable uptake
    for react in organic_ex_rxns:
        model.reactions.get_by_id(react).lower_bound = 0.0

    #############
    # check whether bounds were set to 0.0
    # for r in organic_ex_rxns:
    #     print(model.reactions.get_by_id(r).bounds)
    ############

    return model


##############################################################################
# not implemented in this version --> see MATLAB code
def set_media_ex_bounds(model, ex_rxns):
    # for rxns in ex_rxns:
    #
    #     # set all lower and upper bounds to 0
    #     model.reactions.get_by_id(rxns).lower_bound = 0.0
    #     model.reactions.get_by_id(rxns).upper_bound = 0.0
    print("Not implemented yet.")
    return model
###############################################################################
