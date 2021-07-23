__author__ = "Nantia Leonidou"
__description__ = " Identification of all exchange reactions in a given Cobra model"


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
