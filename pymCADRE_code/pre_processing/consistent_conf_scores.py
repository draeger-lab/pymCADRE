__author__ = "Nantia Leonidou"
__description__ = " Confidence Scores of consistent Recon3D"

from cobra import *
from cobra.io.sbml import *
import pandas as pd

#############################################################################
#  Confidence Scores for rxns in generic model ##############################
#############################################################################

# read model
model = read_sbml_model('dataset/Recon3_blood_consistent.xml')
# model.remove_reactions(['VBOF'])
# io.write_sbml_model(model, "Recon3_blood_consistent.xml")
print('Human Model is imported ... ')

# read initial csv
file = pd.read_csv('dataset/Recon3D_confidence_scores.csv', sep=',')

consistent_rxns = [r.id for r in model.reactions]
# print(consistent_rxns)

file = file.loc[file['Reaction ID'].isin(consistent_rxns)]
print(file)
print(len(model.reactions))
# store file
file.to_csv('Recon3D_consistent_conf_scores.csv', index=False)

