import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json

path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data'
path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml'
path_out = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_h1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\cluster_pheno-ml'

os.chdir(path_out)

labels = np.load('labels_nclus_bio-effect.npy')

os.chdir(path_data_file)

metadata_file = 'metadata-pheno-ml.json'

with open(metadata_file) as output_file:
    metadata = json.load(output_file)

drugs = np.array(metadata['drug'])
cells = np.array(metadata['cell'])
concs = np.array(metadata['conc'])

data_labels = {'labels' : labels,'drugs': drugs, 'cells':cells, "concs": concs}

data_labels = pd.DataFrame(data_labels)

data_labels["conf_100"] = np.nan
data_labels["GR"] = np.nan

os.chdir(path_data_file)

segmentation_data = pd.read_csv("growth_curves.csv")

# remove unnecessary columns

segmentation_data = segmentation_data[['Time', 'Conf',  "Drug", "Final_conc_uM", "cell"]]

for drug in set(drugs):
    idx_drug = [True if x == drug else False for x in drugs ]

    drug_conc = concs[idx_drug]
    drug_cell = cells[idx_drug]

    max_conc = max(drug_conc)

    idx_max_conc = [True if x == max_conc else False for x in drug_conc]

    drug_data = segmentation_data[segmentation_data["Drug"] == drug]

    combined = [(cell, conc) for cell in set(drug_cell) for conc in set(drug_conc)]

    for x in combined:
        cell = x[0]
        conc = float(x[1])
        time = 100

        drug_data = segmentation_data[segmentation_data["Drug"] == drug]
        drug_data = drug_data[drug_data["cell"] == x[0]]
        drug_data = drug_data[drug_data["Final_conc_uM"] == conc]
        drug_data = drug_data[drug_data["Time"] == time]

        data_labels.at[(data_labels['cells'] == cell) &
                    (data_labels['drugs'] == drug) &
                    (data_labels['concs'] == str(conc)),'conf_100'] = drug_data['Conf'].mean()

os.chdir(path_fig+'\\confluence_by_cluster')

for drug in set(drugs):

    sub_data = data_labels[data_labels['drugs']==drug]

    noEffect = sub_data[sub_data['labels']==0]['conf_100']
    mixedEffect = sub_data[sub_data['labels']==1]['conf_100']
    cytostatic = sub_data[sub_data['labels']==2]['conf_100']
    cytotoxic = sub_data[sub_data['labels']==3]['conf_100']

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.boxplot([noEffect.dropna(), mixed.dropna(), cytostatic.dropna(), cytotoxic.dropna()],
               labels=['No Effect', 'Mixed Effect', "Cytostatic", "Cytotoxic"])
    plt.ylabel('Confluence')
    plt.title(drug)
    plt.xlabel('Type of effect')
    plt.savefig(drug+'_confluence_by_cluster.png')

    #todo understand why lower concentration is now shown




