__author__ = "Nantia Leonidou"
__description__ = " Pre-processing of data from GEO Samples. Result files will be used to compute ubiquity scores "

import pandas as pd
from cobra.io.sbml import *

#########################
# read generic model
#########################
model = read_sbml_model('dataset/Recon3D.xml')
print('Human Model is imported ...  DONE')
# store genes in generic model
generic_genes = []
for g in model.genes:
    generic_genes.append(g.name)  # no gene name starts with LOC in generic model
# print('Number of genes in generic model (Recon3D): %s' %(len(generic_genes)))
print('GENERIC MODEL GENES:', generic_genes)

os.chdir('dataset/1_GPL570_GSE3397')
#os.chdir('dataset/2_GPL571_GSE6802')

curr_files = os.listdir(os.getcwd())
count = 1
for f in curr_files:
    if '_GSM' in f:
        ####################################
        # read sample data (controls)
        ####################################
        print('File [%s] is analysed ... ' % f)
        t_1 = f
        # read txt file as df
        df_1 = pd.read_csv(t_1, sep='\t', skiprows=3)
        # replace M with A
        df_1['ABS_CALL'] = df_1['ABS_CALL'].replace(['M'], 'A')

        ###################################
        # read file with microarray info
        ###################################
        platform_file_name = 'GPL570-55999.txt'  # --> INFO about this file: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
        #platform_file_name = 'GPL571-17391.txt'

        data = pd.read_csv(platform_file_name, sep='\t', dtype=str, usecols=["ID", "Gene Symbol", "ENTREZ_GENE_ID"])
        # append abs_call column to final df
        abs_calls = df_1['ABS_CALL']
        data = data.join(abs_calls)
        # drop rows with NaN as Entrez IDs
        data = data[data['ENTREZ_GENE_ID'].notna()]

        # print('INITIAL:', data)

        for index, row in data.iterrows():
            curr_symbol = row['Gene Symbol']  # get gene symbols in current df row
            curr_entrez = row['ENTREZ_GENE_ID']  # entrez IDs
            curr_call = row['ABS_CALL']  # presence/absence
            curr_exp_id = row['ID']  # experiment ID

            if '///' in curr_symbol:  # if /// in gene symbol (means that multiple gene symbols exist) --> create separate
                # and individual row for each one

                # remove row from dataframe
                indexNames = data[data['Gene Symbol'] == curr_symbol].index
                # Delete these row indexes from dataFrame
                data = data.drop(indexNames)

                gene_symbols = curr_symbol.split('///')
                entrez_ids = curr_entrez.split('///')
                for sym, ids in zip(gene_symbols, entrez_ids):  # iterate over two lists of same length at the same time
                    # print(sym,ids)
                    # remove space before/after string
                    sym = sym.strip()
                    ids = ids.strip()
                    # create new row to add
                    new_row = {'ID': curr_exp_id, 'Gene Symbol': sym, 'ENTREZ_GENE_ID': ids, 'ABS_CALL': curr_call}
                    # append row to the dataframe
                    data = data.append(new_row, ignore_index=True)

        # print('after duplcating rows:', data)

        #########################################################
        # remove genes from sample that are not in the model
        #########################################################
        for index, row in data.iterrows():
            gene = row['Gene Symbol']
            if gene not in generic_genes:
                # Get row with g as gene symbol
                indexNames = data[data['Gene Symbol'] == gene].index
                # Delete these row indexes from dataFrame
                data = data.drop(indexNames)

        #########################################################################################################
        # add genes to sample, that present in the model but not in the sample and put A (0, not expressed)
        #########################################################################################################
        sample_gene_symbols = list(data['Gene Symbol'])
        # print('sample_gene_symbols', sample_gene_symbols)
        for g in generic_genes:
            if g not in sample_gene_symbols:
                # print('not in sample symbols:', g)
                data = data.append({'Gene Symbol': g, 'ABS_CALL': 'A'}, ignore_index=True)

        # replace: A-->0(not expressed) and P-->1(expressed)
        data['ABS_CALL'] = data['ABS_CALL'].replace({"P": "1", "A": "0"})

        # pd.set_option("display.max_rows", None, "display.max_columns", None)
        # print('FINAL:', data)

        # store whole dataframe to csv file --> use this df to compute ubiquity scores
        data.to_csv('pre_processed_control_' + str(count) + '_table.csv', index=False)
        count += 1

