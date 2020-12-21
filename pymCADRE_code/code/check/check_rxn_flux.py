__author__ = "Nantia Leonidou"
__description__ = "This function uses the heuristic speed-up proposed by Jerby et al. in the MBA paper for performing a pseudo-FVA calculation. "

# from cobra import *
# from check.find_required_rxns import *
import numpy as np


def check_rxn_flux(model, required_rxns):
    """
        Input:
            - model: cobra.io.core.model.Model
            - required_rxns: list of required reactions to check their flux

        Output:
            - inactive_required
    """
    rxn_lst = required_rxns
    inactive_required = []

    while len(rxn_lst):  # while lst is not empty
        # set number of reactions equal to number of input reactions
        num_rxn_lst = len(rxn_lst)

        # modify objective function
        model.objective = rxn_lst
        # optimize using glpk solver
        fba_solution = model.optimize(objective_sense="maximize")
        # get flux values of all reactions as an array
        opt_max = fba_solution.fluxes.values

        # If no solution was achieved while trying to maximize all reactions, skip the next
        # step of checking individual reactions
        if len(opt_max) == 0:
            inactive_required.append(1)
            break  # break the whole while-loop

        # check whether reactions from model are included in input required reactions
        # put 1 if reaction from model is in the list of required reactions
        is_member = [1 if react.id in required_rxns else 0 for react in model.reactions]

        # get fluxes for these specific reactions
        # only the order of numbers differs from MATLAB results
        required_flux = [opt_max[idx] for idx in range(len(is_member)) if is_member[idx]]

        # absolute values of all fluxes
        abs_fluxes = [abs(ele) for ele in required_flux]
        # filter out those with a specific value --> those reactions will be kept
        filtered = [1 if abs_fluxes[idx] >= 1e-6 else 0 for idx in range(len(abs_fluxes))]
        # get related reaction IDs
        active_required = [required_rxns[idx] for idx in range(len(filtered)) if filtered[idx]]

        # update rxn_lst to be the difference between reactions' list and required active reactions
        rxn_lst = (list(set(rxn_lst) - set(active_required)))
        # get how many reactions were removed from required_rxns list --> updated length of rxns_lst
        num_removed = num_rxn_lst - len(rxn_lst)

        if not num_removed:  # if this if-case is fulfilled, throws an out-of-bounds error, similar happens to MATLAB code as well
            # create a random series of numbers from 0 to len(rxn_lst)
            rand_ind = np.random.permutation(len(rxn_lst))

            # extract reaction from rxn_lst based on the position taken from the first number in the randomly generated series
            i = rxn_lst[rand_ind[0]]

            # change objective function
            model.objective = i
            # Maximize reaction i
            fba_solution = model.optimize(objective_sense="maximize")
            # get flux values of all reactions
            opt_max = fba_solution.fluxes
            if len(opt_max) == 0:
                inactive_required.append(i)
                break  # break the whole while-loop

            abs_opt = [abs(ele) for ele in opt_max]
            for n in abs_opt:
                if n < 1e-6:
                    inactive_required.append(i)
                    break

            rxn_lst.remove(i)

    return inactive_required


###### call function to test it
# print("Human model is loading...")
# model = io.mat.load_matlab_model('../../humanModel.mat')
# print("Metabolites.....")
# required_rxns = find_required_rxns(model, '../../metabolites_humanModel.txt')[1]
# check_rxn_flux(model, required_rxns)
