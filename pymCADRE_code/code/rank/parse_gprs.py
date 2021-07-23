__author__ = "Nantia Leonidou"
__description__ = " Parse Gene–protein–reaction relationships (GPRs). "

import re
# from cobra import *


def parse_gprs(model):
    """ Parse Gene–protein–reaction relationships (GPRs)

        Input: model: cobra.io.core.model.Model

        Output: GPR_file with all gene IDs and GPR_rxns with all related reactions

        """

    # get all gene-protein-reaction rules in the given model
    # a string representation of the GPR rules defined in a readable format
    # For example: - 9962.2 is the gene with the Entrez ID: 9962, that encodes 2 transcripts
    #              - 8706.3 is the gene with the Entrez ID: 8706, that encodes 3 transcripts
    all_GPRS = [react.gene_reaction_rule for react in model.reactions]
    # get number of totally involved genes in each rule
    # possible cases:
    #       (1) 'or' and optionally 'and' --> occurrence of 'or' + 1 (cause above the first gene was not considered while counting the number of 'or' words)
    #       (2)  only 'and' --> put 1
    #       (3)  nor 'or' or 'and' --> put 1
    #       (4) empty cell --> put 0
    num_GPRS = []
    for gprs in all_GPRS:
        if 'and' not in gprs and 'or' in gprs or 'and' in gprs and 'or' in gprs:
            num_GPRS.append(gprs.count("or") + 1)
        elif 'and' in gprs and 'or' not in gprs:
            num_GPRS.append(1)
        elif 'and' not in gprs and 'or' not in gprs and gprs != '':
            num_GPRS.append(1)
        else:
            num_GPRS.append(0)

    num_rows = sum(num_GPRS)  # number of total gene identifiers extracted from GPRs
    num_cols = max(num_GPRS)

    # create an empty array with dimension: 5962x653
    GPR_file = [[0 for x in range(num_cols)] for y in range(num_rows)]
    # create an empty array with dimensions 5962x1
    GPR_rxns = [[0 for x in range(1)] for y in range(num_rows)]

    count = 0
    for i in range(len(all_GPRS)):
        if all_GPRS[i] != '':
            # remove all 'or' words and keep only gene identifies
            rule_GPRS = all_GPRS[i].split('or')
            # remove parenthesis from gene identifies
            rule_GPRS = [re.sub(r'[()]', '', rule) for rule in rule_GPRS]
            for j in range(len(rule_GPRS)):
                # remove 'and', when present in the rule
                GPR = rule_GPRS[j].split('and')
                # fill the GPR file with the gene IDs based on the GPR rules
                # If a rule involves 'AND', the gene IDs will be written in the same row but different column
                GPR_file[count] = GPR
                # GPR_file[count][0:len(GPR)] = GPR
                # fill the GPR_rxns file with the corresponding reaction
                GPR_rxns[count] = model.reactions[i].id
                count += 1

    # remove comma-numbers to match format of Entrez IDs and unnecessary space before and after strings
    GPR_file = [[y.split('_')[0] for y in x] for x in GPR_file]  # --> split('_') for R3con3D, split('.') for humanModel
    GPR_file = [[y.strip() for y in x] for x in GPR_file]

    return GPR_rxns, GPR_file

