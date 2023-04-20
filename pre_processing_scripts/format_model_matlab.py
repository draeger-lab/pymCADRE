__author__ = "Nantia Leonidou"
__description__ = " Format Recon3D to be read from Matlab code"


import cobra
from cobra.io import *

# format Recon2.2
model_recon2 = read_sbml_model('dataset/Recon2.2.xml')
for idx,g in enumerate(model_recon2.genes):
    #print('before',model.genes[idx].id)
    new_id = g.id.replace(':','_')
    model_recon2.genes[idx].id.replace(model_recon2.genes[idx].id,new_id)
    model_recon2.genes[idx].name.replace(model_recon2.genes[idx].name, new_id)
    #print('after',model.genes[idx].id)
cobra.io.write_sbml_model(model_recon2, "Recon2.2_formatted.xml")

model= read_sbml_model('dataset/Recon3D.xml')
# format Recon3D
# format reactions
for idx in range(len(model.reactions)):
    react = model.reactions[idx]

    if react.id.startswith('EX_') or react.id.startswith('DM_'):
        new_id = react.id.replace(react.id[-2:], '('+react.id[-1]+')')
        model.reactions[idx].id = new_id

# format metabolites
for j in range(len(model.metabolites)):
    met = model.metabolites[j]
    new_met = met.id.replace(met.id[-2:], '[' + met.id[-1] + ']')
    model.metabolites[j].id = new_met

# format gene names
# for i in range(len(model.genes)):
#     gene = model.genes[i]
#     new_gene = gene.id.replace(gene.id[-4:],'.0')
#     #print(new_gene)
#     model.genes[i].id = new_gene
#     #print(model.genes[idx].id)


# save model in mat file
cobra.io.write_sbml_model(model, "Recon2.2_formatted.xml")

