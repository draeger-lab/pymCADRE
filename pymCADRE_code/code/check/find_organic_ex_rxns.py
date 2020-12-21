__author__ = "Nantia Leonidou"
__description__ = """ Detection of all organic reactions given a cobra model and a list of exchange reactions 
                   (sink, demand and exchange reactions) """


# from cobra import *
# from cobra.io.sbml import *
# from check.find_ex_reactions import *

def find_organic_ex_rxns(model, ex_rxns):
    """ Detection of organic reactions included in all exchange reactions

        Input:
            - model: cobra.io.core.model.Model
            - ex_rxns: list of reaction IDs

        Output:
            - organic_ex_rxns: list of organic reactions' IDs"""

    # Organic metabolites contain carbon (C) and hydrogen (H) Thus, check molecular formulae
    organic_metabolites = []
    for met in model.metabolites:
        if "C" in met.formula and "H" in met.formula:
            organic_metabolites.append(met.id)
        else:
            # add metabolites, that are organic but don't contain only "C" and "H" in their chemical formulae OR include protein compounds
            if ("C" in met.formula and "R" in met.formula) or ("X" in met.formula and "H" in met.formula) \
                    or ("C" in met.formula and "X" in met.formula) or ("H" in met.formula and "R" in met.formula):
                organic_metabolites.append(met.id)

    # Find all organic reactions (i.e reactions evolving organic metabolites)
    all_organic_reactions = []
    for met in organic_metabolites:
        # all organic reactions associated with current organic metabolite
        org_rxns = model.metabolites.get_by_id(met).reactions
        for react in org_rxns:
            # put all of them in a list, after excluding reactions that may come up multiple times
            if react.id not in all_organic_reactions:
                all_organic_reactions.append(react.id)

    # extract all exchange reactions that involve organic metabolites
    # organic_ex_rxns = [react for react in ex_rxns if react in all_organic_reactions]
    organic_ex_rxns = list(set(all_organic_reactions) & set(ex_rxns))

    return organic_ex_rxns

### test script
# model = read_sbml_model('../../RECON1.xml')
# ex_rxns = find_ex_rxns(model)
# print(len(find_organic_ex_rxns(model, ex_rxns)))
# print(find_organic_ex_rxns(model, ex_rxns))

#  Difference to original MATLAB code when using humanModel.mat:  the following reactions were added into the organic_ex_rxns list,
#  as they were found to involve organic metabolites as well:
#   1)  'DM_Asn_X_Ser_Thr_l'
#   2)  'DM_Ser_Thr_l'
#   3)  'DM_Ser_Gly_Ala_X_Gly_l'
#   4)  'EX_peplys_e'
