__author__ = "Nantia Leonidou"
__description__ = " Calculate ubiquity scores of an individual GEO Series"

import os
import pandas as pd

os.chdir('dataset/1_GPL570_GSE3397')
#os.chdir('dataset/2_GPL571_GSE6802')

curr_files = os.listdir(os.getcwd())

present_files = []
abs_calls = []
for f in curr_files:
    if f.startswith('final_control_'):
        present_files.append(f)
        data = pd.read_csv(f)
        abs_calls.append(list(data['ABS_CALL']))
#print(abs_calls)

# ubiquity scores
ubiquity_scores = []
c = len(abs_calls[0]) # all lists in abs_calls have the same length
count = 0
for i in range(c):
    # Change this based on how many files current directory has. Currently it works for four files
    freq = abs_calls[0][i] + abs_calls[1][i] + abs_calls[2][i] + abs_calls[3][i]
    score = freq / len(present_files)
    # keep only two decimals
    formatted = "{:.2f}".format(score)
    ubiquity_scores.append(float(formatted))
print(ubiquity_scores)

ubiq_scores_df = pd.DataFrame(ubiquity_scores)
out_fname = os.getcwd().split('dataset')[-1].replace('\\','')
ubiq_scores_df.to_csv(out_fname+'_ubiquity.csv', index=False, header=False)