__author__ = "Nantia Leonidou"
__description__ = " pymCADRE main"

from time import process_time
from rank.rank_reactions import *
from prune.prune_model import *
import pandas as pd
import numpy as np
from cobra import *
from cobra.io.sbml import *


# read model
model = io.read_sbml_model('../test_dataset/recon1_with_BOF_and_VBOF.xml')
# genes
G = pd.read_csv('../test_dataset/1_GPL570_GSE3397_entrez_ids.csv')
G = list(G['ENTREZ_GENE_ID'])
# ubiquity scores
U = pd.read_csv('../test_dataset/1_GPL570_GSE3397_ubiquity.csv', header=None)
U = U.rename(columns={0: "Scores"})
U = list(U['Scores'])
#confidence_scores
confidence_scores = pd.read_csv('../test_dataset/Recon1_confidence_scores_with_BOF_and_VBOF.csv')
confidence_scores = np.float64(list(confidence_scores['Confidence Score']))

##############################################
# Rank reactions
##############################################

print('Processing inputs and ranking reactions...')
GM, C, NC, P, Z, model_C = rank_reactions(model, G, U, confidence_scores, [], method=1)

##################################################
# Define inputs to the model pruning step
##################################################
# Define core vs. non-core ratio threshold for removing reactions
eta = 1 / 3
# Check functionality of generic model
genericStatus = check_model_function(deepcopy(GM), 'required_mets', '../test_dataset/key_metabolites_RECON1.txt')[0]

if genericStatus:
    print('Generic model passed precursor metabolites test')

    ##############################################################
    # If generic functionality test is passed, prune reactions
    ###############################################################
    print('Pruning reactions...')
    t0 = process_time()
    PM, cRes = prune_model(GM, P, C, Z, eta, '../test_dataset/key_metabolites_RECON1.txt', salvage_check=1, method=1)
    # Stop the stopwatch / counter
    t_stop = process_time()
    # compute elapsed time
    pruneTime = t_stop - t0

else:
    print('Generic model failed precursor metabolites test!!')

print('PM.reactions: ',len(PM.reactions))
print('PM.metabolites: ',len(PM.metabolites))
print('PM.genes: ', len(PM.genes))
print('prune Time: ', pruneTime)

#########################################
# * store pruned model in SBML format *
#########################################
io.write_sbml_model(PM, "python_recon1_BEC1_method1.xml")
