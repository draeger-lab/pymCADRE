__author__ = "Nantia Leonidou"
__description__ = " Detect required reactions and add them to the model"

from cobra import *


def add_metabolites(met_id_lst, coeffs_lst, compartment_lst):
    """ Create a dictionary with metabolites and their coefficients
        Input: lst with metabolites' IDs, list with repsective coefficients
               and lst with respective compartment letters
        Output: dictionary with metabolites"""

    if len(met_id_lst) != len(coeffs_lst) or len(met_id_lst) != len(compartment_lst) or len(coeffs_lst) != len(
            coeffs_lst):
        print("ERROR: Input lists should have the same size")
    else:
        met_dict = {}
        for i in range(len(met_id_lst)):
            met_dict[Metabolite(met_id_lst[i], compartment=compartment_lst[i])] = coeffs_lst[i]
        return met_dict


def add_single_reaction(model, id, lb, up, met_dict):
    """ Add reactions to a model.
        Input: model, reaction id, lower/upper bound and dictionary with the metabolites to add
        Output: model with a newly added reaction and the respective reaction."""

    react = Reaction()
    react.id = id
    react.lower_bound, react.upper_bound = lb, up
    # add metabolites by hand in the model, First checks whether the metabolites exist in the model
    react.add_metabolites(met_dict)

    return model, react


def find_required_rxns(model, metabolites):
    """ Input: model and list of metabollites in a .txt format (see example metabolites_humanModel.txt)
        Output: model and list of only required reactions' id """

    # define manually _coa metabolites
    coa_mets = ['accoa_m', 'succoa_m', 'pmtcoa_c']

    # read input file including diverse metabolites
    with open(metabolites, 'r') as file:
        # list of all input metabolites
        met_lst = [line.rstrip('\n') for line in file]

    # difference between coa_mets and input list of metabolites --> list of strings
    diff = (list(set(met_lst) - set(coa_mets)))
    print(diff)

    required_rxns = []
    # add metabolites from diff as demand reactions, if metabolite is in the model
    model_mets = [m.id for m in model.metabolites]
    for met in diff:
        r = model.add_boundary(model.metabolites.get_by_id(met),type="demand")  # add_boundary if metabolite is already in the model
        required_rxns.append(r.id)


    react1 = add_single_reaction(model, "DM_accoa_m",  0.0,1000.0, met_dict={Metabolite('accoa_m', compartment='m'): -1, Metabolite('coa_m', compartment='m'): 1})[1]
    required_rxns.append(react1.id)
    react2 = add_single_reaction(model, "DM_succoa_m", 0.0,1000.0, met_dict={Metabolite('succoa_m', compartment='m'): -1, Metabolite('coa_m', compartment='m'): 1})[1]
    required_rxns.append(react2.id)
    react3 = add_single_reaction(model, "DM_pmtcoa_c", 0.0,1000, met_dict={Metabolite('pmtcoa_c', compartment='c'): -1, Metabolite('coa_c', compartment='c'): 1})[1]
    required_rxns.append(react3.id)
    model.add_reactions([react1,react2,react3])

    return model, required_rxns

