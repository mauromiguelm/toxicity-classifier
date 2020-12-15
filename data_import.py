import os, json, re
import numpy as np
import pandas as pd
from sklearn import preprocessing

def get_distance_metrics(path_in, path_out):
    """ Import distance metrics of drugs to control, and convert to array"""

    dist_files = []

    for path, subdirs, files in os.walk(path_in):
        for name in files:
            dist_files.append(os.path.join(path, name))

    dist_files = [x for x in dist_files if re.search("euclidean",x)]

    data = [json.load(open(x)) for x in dist_files]

    drug = []
    cell = []
    conc = []

    store_data = np.empty(shape=(0,100))

    min_length = store_data.shape[1]

    for idx in range(0,len(data),1):

        file = dist_files[idx]

        metadata = re.compile("\\\\")

        metadata = metadata.split(file)

        data_sub = data[idx]

        time = data_sub["time"]

        keep = [x > 0 for x in time]

        del data_sub['time']

        idx_conc = [*data_sub.keys()]

        conc.append(idx_conc)

        len_conc = len(idx_conc)

        drug.append([re.split("_", metadata[3])[0]] * len_conc)

        cell.append([re.split("_", metadata[2])[0]] * len_conc)

        len_data = max(map(len, data_sub))

        data_sub = np.array([x + [None] * (len_data - len(x)) for x in data_sub.values()])

        time = [i for data_sub, i in enumerate(time) if keep[data_sub] == True]

        data_sub = data_sub[:,keep]

        if data_sub.shape[1] < min_length:
            min_length = data_sub.shape[1]
            store_data = store_data[:,1:min_length+1]
        elif data_sub.shape[1] > min_length:
            data_sub = data_sub[:,1:min_length+1]

        store_data = np.append(store_data, data_sub, axis=0)

        del keep

    drug = [item for sublist in drug for item in sublist]
    cell = [item for sublist in cell for item in sublist]
    conc = [item for sublist in conc for item in sublist]

    metadata = {'drug' : drug, 'cell' : cell, 'conc' : conc}

    os.chdir(path_out)

    with open('metadata-pheno-ml.json', 'w') as output_file:
        json.dump(metadata, output_file)

    np.save("dist_combined", store_data)

def scale_results(metadata,
                  output_file,
                  data_file = 'dist_combined.npy',
                  path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data',
                  path_out = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data',
                  apply_denoise = True):
    "scale distance matrix by row, to avoid biass in clustering"

    os.chdir(path_data_file)
    data = np.load(data_file)
    cells = metadata['cell']

    if apply_denoise == True:

        data_pd = pd.DataFrame(data.T)

        data_pd = data_pd.rolling(window=5).median().T

        data = data_pd.iloc[:,4:data_pd.shape[1]-1]

        data = np.asarray(data)

    else:
        pass

    for cell in set(cells):

        idx = [x == cell for x in cells]

        idx = np.array(idx, dtype='bool')

        data_cells = data[idx]

        ncol = data_cells.shape[0]

        nrow = data_cells.shape[1]

        data_cells = data_cells.reshape(ncol*nrow,1)

        scaler = preprocessing.MinMaxScaler().fit(data_cells)

        scaled_data = scaler.transform(data_cells)

        scaled_data = scaled_data.reshape(ncol,nrow)

        data[idx] = scaled_data


    os.chdir(path_out)

    np.save(output_file, data)

if __name__ == '__main__':

    path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data'
    os.chdir(path_data_file)
    metadata_file = 'metadata-pheno-ml.json'

    with open(metadata_file) as output_file:
        metadata = json.load(output_file)

    get_distance_metrics(path_in='//d.ethz.ch/groups/biol/sysbc/sauer_1/users/Mauro/from_Andrei/distances/cropped',
                         path_out='\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data')


    scale_results(data_file = 'dist_combined.npy',
                  output_file = "dist_combined_cell-scaled_denoise",
                  path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data',
                  metadata = metadata,
                  path_out = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data')
