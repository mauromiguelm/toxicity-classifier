
import os, json, re
import numpy as np
from sklearn import preprocessing

def get_distance_metrics(path_in, path_out):
    """ Import distance metrics of drugs to control, and convert to array"""

    dist_files = []

    for path, subdirs, files in os.walk(path_main):
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

def scale_results(data_file = 'dist_combined.npy',
                  path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data',
                  path_out):
    "scale distance matrix by row, to avoid biass in clustering"

    os.chdir(path_data_file)

    data = np.load(data_file)

    scaler = preprocessing.MinMaxScaler().fit(data.T)

    scaled_data = scaler.transform(data.T).T

    os.chdir(path_out)

    np.save("dist_combined_scaled", scaled_data)

if __name__ == '__main__':

    get_distance_metrics('//d.ethz.ch/groups/biol/sysbc/sauer_1/users/Mauro/from_Andrei/distances/cropped',
                         '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data')

    scale_dist_results("dist_combined.npy",
                       '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data',
                       '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data')




