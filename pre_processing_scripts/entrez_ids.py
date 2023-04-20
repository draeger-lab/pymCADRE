__author__ = "Nantia Leonidou"
__description__ = " Use databse call to find missing entrez IDs "

import os
import pandas as pd
from cobra.io.sbml import *
import requests
import json

os.chdir('dataset/1_GPL570_GSE3397')
#os.chdir('dataset/2_GPL571_GSE6802')

curr_files = os.listdir(os.getcwd())

count = 1
for f in curr_files:
    if f.startswith('pre_processed_control_'):
        print('File [%s] is analysed ... ' % f)
        data = pd.read_csv(f, sep=',')

        is_NaN = data.isnull()
        row_has_NaN = is_NaN.any(axis=1)
        rows_with_NaN = data[row_has_NaN]

        for index, row in rows_with_NaN.iterrows():
            gene = row['Gene Symbol']

            # make database call to find missing Entrez IDs if possible
            try:
                # call gene from BiGG database
                bigg_response = requests.get(
                    "http://bigg.ucsd.edu/api/v2/search?query=" + gene + "&search_type=genes")  # result included Recon3D, Recon1 and more
                bigg_json = bigg_response.content.decode("utf-8")
                info = json.loads(bigg_json)

                res = {}
                for i in range(len(info['results'])):
                    # only for Recon3D
                    if info['results'][i]['model_bigg_id'] == 'Recon3D':
                        if info['results'][i]['name'] == gene:
                            res = info['results'][i]

                entrez_id = res['bigg_id'].split('_')[0].strip()

            except:
                if gene[0].isdigit():
                    entrez_id = gene.split('_')[0].strip()
                else:
                    print('%s can not be found' % gene)

            data.loc[index] = ['NaN', gene, entrez_id, 0]

        # save dataframe
        data.to_csv('final_control_' + str(count) + '.csv', index=False)
        count += 1

# write ONLY entrez IDs to a file
for f in curr_files:
    if f.startswith('final_'):
        print('File [%s] is read ... ' % f)
        data = pd.read_csv(f, sep=',')
        all_entrez = data['ENTREZ_GENE_ID']

        out_fname = os.getcwd().split('dataset')[-1].replace('\\', '')
        print('out_fname', out_fname)
        all_entrez.to_csv(out_fname + '_entrez_ids.csv', index=False)
        print('all_entrez', all_entrez)
        break
