<img src="pymCADRE_logo.png" alt="drawing" width="400"/>

## Overview

The **pymCADRE**  constructs context-specific models by leveraging gene  expression  
data and literature-based evidence, along with network topology information. The reactions withing the global generic 
model are being ranked and the ones with the lowest supporting evidence for the tissue
of interest are given the highest priority for removal:
```
GM, C, NC, P, Z, model_C = 
        rank_reactions(model, G, U, confidence_scores, [], method=1)
```
If the generic functionality test is passed, the model undergoes pruning, which results
to a context-specific reconstruction:
```
PM, cRes = 
    prune_model(GM, P, C, Z, eta, precursorMets, salvage_check=1, method=1)
```

## Prerequisites

This tool has the following dependencies:

python 3.8

Packages:
* pandas
* numpy
* cobra
* requests
* os

## Input data
+ `model`: generic global metabolic model of the organism of interest
+ `precursorMets`: list of precursor, key, metabolites
+ `confidence_scores`: literature/experimental-based confidence assigned to reactions in generic model
###### Tissue-specific expression evidence
+ `G`: list of Entrez IDs for all genes 
+ `U`: list of genes' ubiquity scores

#### Optional Inputs
+ `salvageCheck`: (1) test tissues on  whether  they  synthesize  purines from  purine  bases  
and  Phosphoribosyl  pyrophosphate  (PRPP)  through  the salvage pathway, (0) no test.
+ `C_H_genes`: list with Entrez IDs for genes with particularly high evidence of activity in the tissue of interest
+ `method`: method for internal optimization, (1) fastFVA or (2) fastcc

## Outputs
+ `PM`: pruned COBRA tissue-specific model
+ `GM`: COBRA model after removing blocked reactions from the input global model
+ `C`: core reactions in `GM`
+ `NC`: non-core reactions in `GM` 
+ `Z`: reactions with zero expression across all samples after binarization
+ `model_C`: core reactions in the original global model (including blocked reactions)
+ `pruneTime`: total reaction pruning time 
+ `cRes`: result of model checks (consistency/function) during pruning
    + **1**: Reaction has zero expression evidence.
    + **-1.x**: Reaction removed along with any inactivated reactions
      (x: removed core reactions)
    + **+1.x**: No reactions removed.
      + x = 1: ratio of inactivated core vs. non-core reactions above a pre-defined threshold
      + x = 0: production of precursor metabolites not possible 
    + **2**: Reactions has expression evidence
    + **-2.x**: Reaction removed along with only non-core inactivated reactions
      (x: removed core reactions, always zero!)
    + **+2.x**: No reactions removed.
      + x= 1: removal led to inactive cores 
      + x= 0: production of precursor metabolites not possible
    + **3**: Removal or reaction prevented the production of required metabolites. No reaction removed.
  

## Usage
To run pymCADRE simply execute the notebook named main_pymcadre.ipynb or or the python script 
named pymcadre.py. One can also modify one of the test_*.py scripts to adjust
the parameters as wanted and read the desired input files. Jupyter notebooks with
test-runs are also provided.

## Authors
**Nantia Leonidou**


