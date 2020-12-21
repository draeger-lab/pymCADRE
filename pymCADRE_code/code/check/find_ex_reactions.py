__author__ = "Nantia Leonidou"
__description__ = " Identification of all exchange reactions in a given Cobra model"

########
# NOTE:
# Using the cobrapy functions to obtain the exchange, sink and demand reactions (model.sink, model.exchanges and model.demand), two reactions
# (DM_T_antigen(g) and sink_citr-L) of the mouse model could not be identified as exchange reactions.
# However, using the humanModel.mat, all reactions have been detected correctly.
# --> Thus, the identification of all exchanges reactions is done manually, using the length of the reactions' metabolites
#     as criterion, to ensure always the correct result.
########
from cobra import *


def find_ex_rxns(model):

    """ This function identifies reactions that exchange metabolites into and out of the extra- and intracellular space
        (sinks, demands and exchanges)

        Input:
            - model: cobra.io.core.model.Model

        Ouput:
            - ex_rxns_ids: list of all identified reactions' IDs """

    ex_rxns_ids = []

    for react in model.reactions:

        if len(react.metabolites) == 1:
            ex_rxns_ids.append(react.id)

    return ex_rxns_ids

#model = io.mat.load_matlab_model('../../humanModel.mat')
#print(find_ex_rxns(model))