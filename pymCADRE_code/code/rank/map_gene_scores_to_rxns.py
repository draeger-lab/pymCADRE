__author__ = "Nantia Leonidou"
__description__ = " Map gene ubiquity scores to reactions. "

# from cobra import *
# from rank.parse_gprs import *
# from test_Inputs.test_Inputs import *
import numpy as np


def map_gene_scores_to_rxns(model, G, U, GPR_file):
    """ Map gene ubiqiuty scores (i.e., how often a gene is expressed in tissue samples) to reactions

        Input:
            - model --> cobra.io.core.model.Model
            - G --> column vector, where each entry is the human- or mouse-specific Entrez ID for gene
            - U --> numeric column vector, where each entry represents the calculated "ubiquity score" for the corresponding gene
            - GPR_file

        Output: U_GPR --> list of lists, with ubiquity scores as float numbers"""

    # get all Entrez IDs from given model
    genes = []
    for g in model.genes:
        g = g.id
        g = g.split('_')[0]
        g = g.strip()
        genes.append(g)
    genes = [int(g) for g in genes]

    l = [len(GPR_file[i]) for i in range(len(GPR_file))]
    # get the maximum and use it as number of columns, that C_H_GPR should maximally have
    columns = max(l)
    # create list of lists
    U_GPR = [[np.nan for x in range(columns)] for y in range(len(GPR_file))]

    # get intersection of genes and G and their indices in G
    gene_inter = list(set(G) & (set(genes)))
    # which genes are in the intersection, sorted
    gene_id = sorted([str(g) for g in gene_inter])
    # indexes of those genes in G
    G_idx = [G.index(int(x)) for x in gene_id]

    # *  faster than the corresponding MATLAB part, with small modifications (gene_GPR doesn't need to be computed) *
    # U_GPR contains the ubiquity scores for those genes that are found to be as
    # well in the GPR file as in the gene_idx list
    for gpr_idx in range(len(GPR_file)):
        for gpr_idx_nested in range(len(GPR_file[gpr_idx])):
            if GPR_file[gpr_idx][gpr_idx_nested] in gene_id:
                get_idx = gene_id.index(GPR_file[gpr_idx][gpr_idx_nested])
                U_GPR[gpr_idx][gpr_idx_nested] = U[G_idx[get_idx]]
    # now penalize genes with zero expression, such that corresponding reactions
    # will be ranked lower than non-gene associated reactions.
    for x in range(len(U_GPR)):
        for y in range(len(U_GPR[x])):
            if U_GPR[x][y] == float(0):
                U_GPR[x][y] = -1e-6

    return U_GPR


### test script
# model = io.mat.load_matlab_model('../../humanModel.mat')
# C_H_genes = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[2]
# G = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[0]
# U = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[1]
# GPR_file = parse_gprs(model)[1]
# GPR_rxns = parse_gprs(model)[0]
# print('U_GPR:', map_gene_scores_to_rxns(model, G, U, GPR_file))
