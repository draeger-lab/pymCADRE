__author__ = "Nantia Leonidou"
__description__ = " Check the salvage pathway. " \
                  " Test the ability of all tissues to synthesize purine nucleotides from purines bases and PRPP, " \
                  " as de novo purine synthesis occurs primarily in the liver and other tissues uses salvage pathway."

from cobra import *
from check.find_ex_reactions import *
from check.set_metabolite_bounds import *
from time import process_time


def check_salvage_path(model):
    """  * not useful for tissues known to make purines de novo *
        When the model is allowed to use PRPP and guanine or hypoxanthine, test if it
        can make GMP or IMP. This is the salvage pathway that non-hepatic tissues use
        for purine synthesis.

        Input:
            - model: cobra.io.core.model.Model

        Output:
            - salvage_status and time of execution"""

    # Start the stopwatch / counter
    # returns float value of time in seconds
    t1_start = process_time()

    # find exchange reactions in the model
    ex_rxns = find_ex_rxns(model)

    # turn off uptake of organic metabolites
    if "media_def" in globals():  # check if variable named media_def is defined
        model = set_media_ex_bounds(model, ex_rxns)  # not implemented in this version
    else:
        model = set_organic_met_bounds(model, ex_rxns)

    # add PRPP sink reaction for subsequent checks
    try:
        model.metabolites.get_by_id("prpp_c")
        model.add_boundary(model.metabolites.get_by_id("prpp_c"), type="sink", lb=-5, ub=5)
    except:
        met_lst = [Metabolite("prpp_c", compartment='c')]
        model.add_metabolites(met_lst)
        model.add_boundary(model.metabolites.get_by_id("prpp_c"), type="sink", lb=-5, ub=5)

    # store reaction names
    names = []  # rxn name
    for r in model.reactions:
        names.append(r.name)

    ##############################
    # Check production of GMP:
    ##############################
    # allow uptake of guanine
    model_gmp = model.copy()
    if "Guanine exchange" in names:  # Reaction id: EX_gua_e
        idx = names.index("Guanine exchange")
        model_gmp.reactions[idx].lower_bound = -5.0
    # add demand reaction for GMP
    try:
        model_gmp.metabolites.get_by_id("gmp_c")
        gmp_react = model_gmp.add_boundary(model_gmp.metabolites.get_by_id("gmp_c"), type="demand")
    except:
        met_lst = [Metabolite("gmp_c", compartment='c')]
        model_gmp.add_metabolites(met_lst)
        gmp_react = model_gmp.add_boundary(model_gmp.metabolites.get_by_id("gmp_c"), type="demand")
    # modify objective function
    model_gmp.objective = {gmp_react: 1}
    # optimize using glpk solver
    sol = model_gmp.optimize(objective_sense="maximize")
    # True or False if objective value is greater than 1e-6
    status_gmp = sol.objective_value > 1e-6

    ##############################
    # Check production of IMP:
    ##############################
    # allow uptake of hypoxanthine
    model_imp = model.copy()
    if "Hypoxanthine exchange" in names:  # Reaction id: EX_hxan_e
        idx = names.index("Hypoxanthine exchange")
        model_imp.reactions[idx].lower_bound = -5.0
    # add demand reaction for IMP
    try:
        model_imp.metabolites.get_by_id("imp_c")
        imp_react = model_imp.add_boundary(model_imp.metabolites.get_by_id("imp_c"), type="demand")
    except:
        met_lst = [Metabolite("imp_c", compartment='c')]
        model_imp.add_metabolites(met_lst)
        imp_react = model_imp.add_boundary(model_imp.metabolites.get_by_id("imp_c"), type="demand")
    # modify objective function
    model_imp.objective = {imp_react: 1}
    # optimize using glpk solver
    sol = model_imp.optimize(objective_sense="maximize")
    # True or False if objective value is greater than 1e-6
    status_imp = sol.objective_value > 1e-6

    # get final salvage status
    salvage_status = status_gmp & status_imp

    # Stop the stopwatch / counter
    t1_stop = process_time()
    # compute elapsed time
    time = t1_stop - t1_start

    return salvage_status, time


### test script
# from cobra.io.sbml import *
# model = read_sbml_model('../../RECON1.xml')
# # model = io.mat.load_matlab_model('../../humanModel.mat')
# print(check_salvage_path(model))
