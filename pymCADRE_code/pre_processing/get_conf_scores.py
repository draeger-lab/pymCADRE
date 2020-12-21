__author__ = "Nantia Leonidou"
__description__ = " Obtain confidence scores for reactions from generic model by accessing the API of Virtual Metabolic Human  "

import json
import requests
from cobra.io.sbml import *
import pandas as pd

#############################################################################
#  Confidence Scores for rxns in generic model ##############################
#############################################################################

# read model
model = read_sbml_model('dataset/Recon2.2_formatted.xml')
print('Human Model is imported ... ')

# format reaction names in Recon2
# for idx,g in enumerate(model.reactions):
#     if g.id == 'EX_3aib_D_e':
#         model.reactions[idx].id = 'EX_3aib__D_e'
#     if g.id.startswith('EX_') and g.id.endswith('_b'):
#         new_id = g.id.replace(g.id[-2:],'_e')
#         model.reactions[idx].id = new_id
# # print([i.id for i in model.reactions])

# format reaction names in Recon1
# remove _copy2 reactions as they are the same as _copy1 but reversible
# for idx,g in enumerate(model.reactions):
#     if g.id.endswith('_copy1'):
#         new_id = g.id.replace(g.id[-6:],'')
#         model.reactions[idx].id = new_id
# # print([i.id for i in model.reactions])


# define null, true and false to avoid "NameError: name 'null' is not defined"
null = None
true = True
false = False

# store conf. scores for generic model's reactions, value is None if no conf.score is assigned
confidence_scores = {}

print('Starting database calls ...')
c=0
for idx,r in enumerate(model.reactions):
    print('Reaction ID:', r.id)
    try:
        # make a get request to BiGG models
        # bigg_response = requests.get(
        #     "http://bigg.ucsd.edu/api/v2/models/Recon3D/reactions/" + r.id)
        # bigg_json = bigg_response.content.decode("utf-8")
        # data = json.loads(bigg_json)
        # # store respective old identifier
        # old_id = data["old_identifiers"]

        bigg_response = requests.get("http://bigg.ucsd.edu/api/v2/universal/reactions/" + r.id)
        bigg_json = bigg_response.content.decode("utf-8")
        info = json.loads(bigg_json)
        # print('info:',info)
        old_id = info["old_identifiers"]

        for i in old_id:
            # call reaction from VMH using the respective old identifier
            vmh_res = requests.get(
                "https://www.vmh.life/_api/reactions/?abbreviation=" + i)

            vmh_json = vmh_res.content.decode("utf-8")
            vmh_confidence = json.loads(vmh_json)
            score = 0
            if vmh_res.status_code == 200 and len(vmh_confidence['results']):
                #print(vmh_res.status_code)

                # replace None to 0
                if vmh_confidence['results'][0]['mcs'] is None:
                    vmh_confidence['results'][0]['mcs'] = 0

                if r.id in list(confidence_scores.keys()):
                    print('ALREADY IN', r.id)

                # store confidence scores to the corresponding 'initial' identifier (found in BiGG model)
                score = vmh_confidence['results'][0]['mcs']
                print('Confidence scores:',vmh_confidence['results'][0]['mcs'],i)
                break
            # elif vmh_res.status_code == 200 and len(vmh_confidence['results'])==0:
            #     print('Problem1:',vmh_confidence['results'],old_id,i)
            # else:
            #     print('Problem:',vmh_confidence['results'],old_id,i)

        print(r.id, 'successfull')
        confidence_scores[r.id] = score
        c+=1

    except:
        confidence_scores[r.id] = 500
        print('rxn [%s] not found' % r.id)
        print('-------------------')
        c+=1

print('Database calls are done ...')

print('Dictionary with confidence scores:', confidence_scores)
print('c',c)

# ##################
# ## Save File #####
# ##################
print('Saving output file ... ')
# convert dictionary to dataframe for easier storing
confscores_df = pd.DataFrame(list(confidence_scores.items()), columns=['Reaction ID', 'Confidence Score'])
confscores_df.to_csv('Recon2_confidence_scores.csv', index=False)

print('len conf scores', len(confscores_df['Confidence Score']))
print('Retrieval of confidence scores: DONE ...')
