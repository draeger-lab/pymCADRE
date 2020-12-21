__author__ = "Nantia Leonidou"
__description__ = " Calculate connectivity evidence "

from cobra import *
import numpy as np


def calc_conn_evidence(model, E_X):
    """ Input:
            - model: cobra.io.core.model.Model --> GM model computed from initialize_generic_model
            - E_X: list of E_X computed from initialize_generic_model

        Output:
            - E_C: list of connectivity scores """

    # create stoichiometric matrix
    s_mtrx = util.array.create_stoichiometric_matrix(model)
    # S matrix should be binarized to indicate metabolite participation in each reaction
    S_bin = abs(s_mtrx.astype(np.int64))

    # Adjacency matrix (i.e., binary reaction connectivity); the connection between
    # a reaction and itself is ignored by subtracting the identity matrix
    S_bin_transpose = np.transpose(S_bin)  # transpose matrix S_bin
    mult = S_bin_transpose.dot(S_bin)  # multiply both matrices, takes a bit time (!)
    for i in range(len(mult)):
        for j in range(len(mult[i])):
            # put in matrix 1, if cell contains non-zero element
            if mult[i][j] != 0:
                mult[i][j] = 1
    # subtract the identity matrix
    id_mtrx = np.identity(S_bin.shape[1])
    # final Adjacency matrix
    A = mult - id_mtrx

    # Influence matrix --> describes the divided connectivity of reactions --
    #                      e.g., if RXN1 is connected to 4 reactions, its influence on each of those reactions would be 0.25
    row_sum = np.sum(A, axis=1)  # sum of elements in each row of matrix A
    rep = np.array([row_sum, ] * A.shape[
        0]).transpose()  # construct an array by repeating column vector row_sum as many times as shape of matrix A
    I = A / rep  # define Influence Matrix

    # Weighted influence matrix --> each reaction's influence on others is weighted by its expression score. Thus, reactions with with 0 in the
    #                               corresponding cell, have no weighted influence on adjacent reactions
    rep2 = np.array([E_X, ] * A.shape[0]).transpose()
    # element-wise multiplication
    WI = np.multiply(rep2, I)

    # Connectivity score --> sum the influence of all connected reactions
    E_C = np.nansum(WI,axis=0).tolist()

    return E_C
