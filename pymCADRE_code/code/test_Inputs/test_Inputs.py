__author__ = "Nantia Leonidou"
__description__ = " Create Test Inputs in the correct format "

import scipy.io
from cobra import *


def get_test_inputs(file, file2):
    """ Extract data for G, U and C_H_genes from given .mat file containing these info"""
    # read mat file
    mat = scipy.io.loadmat(file)
    model_mat = scipy.io.loadmat(file2)

    G = []
    for g in range(len(mat['G'])):
        G.append(mat['G'][g][0][0])

    U = []
    for u in range(len(mat['U'])):
        U.append(float("{:.4f}".format(mat['U'][u][0])))

    C_H_genes = []
    for ch_gene in range(len(mat['C_H_genes'])):
        C_H_genes.append(mat['C_H_genes'][ch_gene][0][0])

    confidence_Scores = []
    for idx in range(len(model_mat['confidenceScores'])):
        confidence_Scores.append(model_mat['confidenceScores'][idx][0])

    return G, U, C_H_genes, confidence_Scores


# C_H_genes = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[2]
# G = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[0]
# U = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[1]
# confidence_scores = get_test_inputs('../../testInputs.mat','../../humanModel.mat')[3]
#

def get_model_from_mat(file):

    model = io.mat.load_matlab_model(file)

    return model

