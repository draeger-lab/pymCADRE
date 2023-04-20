__author__ = "Nantia Leonidou"
__description__ = """ Check production of required metabolites  """

from time import process_time
from check.find_ex_reactions import *
from check.set_metabolite_bounds import *
from check.find_required_rxns import *
from check.check_rxn_flux import *
import os


def check_model_function(model, *args):
    """ Check production of required metabolites
        Inputs: - a cobra model
                - optional inputs: - 'required_mets', met_lst
                                   - 'biomass', biomass_rxn
                                   - 'media', media_def

        Output:  - required_met_status
                 - time """

    ################################################################################################
    #                               Parse input parameters                                         #
    ################################################################################################
    if len(args) != 0:
        if len(args) % 2 == 0:  # if division between len(args) and 2 gives rest 0

            options = ["required_mets", "biomass", "media"]

            for i in range(0, len(args), 2):  # i goes from 1 to len(args) incrementally in steps of 2
                arg_name = args[i]  # define arg_name, can be 'required_mets' or 'biomass' or 'media'
                arg_val = args[i + 1]  # define arg_val, can be met_lst or biomass_rxn or media_def

                # read input .txt file
                with open(arg_val, 'r') as file:
                    # list of all input metabolites
                    met_lst = [line.rstrip('\n') for line in file]

                # set 1 if one of the options is given as arguments, else set 0
                option = [1 if arg_name == opt else 0 for opt in options]

                # find indices of nonzero elements in the lst, called options
                # returns list with indices depending on whether an input argument from *args matches one of the strings in list 'options'
                # Example: if 'media' is given as an argument in the check_model_function function,
                #          then the find_idx_option = [2], because the 'media' string is located in the index 2 in 'options' list
                find_idx_option = [i for i, e in enumerate(option) if e != 0]

                # if indices of nonzero elements are found
                # biomass_rxn = []
                # media_def = []
                if len(option) != 0:
                    for elem in find_idx_option:
                        if elem == 0:
                            arg_val = met_lst
                        elif elem == 1:
                            raise Exception("%s option not yet implemented." % arg_name)
                            # arg_val = biomass_rxn
                        elif elem == 2:
                            raise Exception("%s option not yet implemented." % arg_name)
                            # arg_val = media_def
                # else-case
                else:
                    raise Exception("Unknown option %s." % arg_name)
        else:
            raise Exception("Incorrect number of input arguments to function.", os.path.basename(__file__))
    ################################################################################################

    # Start the stopwatch / counter
    # returns float value of time in seconds
    t1_start = process_time()

    # Identify exchange reactions in the model
    ex_rxns = find_ex_rxns(model)

    # Turn off uptake of organic metabolites
    if "media_def" in globals():  # check if variable named media_def is defined
        set_media_ex_bounds(model, ex_rxns)  # not implemented in this version
    else:
        set_organic_met_bounds(model, ex_rxns)

    # store reaction names
    names = []  # rxn name
    for r in model.reactions:
        names.append(r.name)

    ###########
    # Allow uptake of glucose and CO2, only if they are in the model
    ###########
    if 'D-Glucose exchange' in names:  # Reaction id: EX_glc__D_e
        idx = names.index('D-Glucose exchange')
        model.reactions[idx].lower_bound = -5.0
        
    if 'CO2 exchange' in names: # Reaction id: EX_co2_e
        idx = names.index('CO2 exchange')
        model.reactions[idx].lower_bound = -1000.0

    # Add demand reactions for required metabolites
    if "met_lst" in globals():
        model, required_rxns = find_required_rxns(model, met_lst)
    else:
        required_rxns = []

    # if "biomass_rxn" in globals():
    #     required_rxns = np.concatenate(required_rxns, biomass_rxn) # not implemented in this version

    inactive_required = check_rxn_flux(model, required_rxns)
    required_met_status = int(not len(inactive_required))  # return negation of number corresponding to length of inactive_required list

    # Stop the stopwatch / counter
    t1_stop = process_time()
    # compute elapsed time
    time = t1_stop - t1_start

    print('check_model_function done ....')
    return required_met_status, time
