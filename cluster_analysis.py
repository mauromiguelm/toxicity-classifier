import os
import numpy as np
import random
from tslearn.clustering import TimeSeriesKMeans
import matplotlib.pyplot as plt

def run_clustering_methods(data_file = 'dist_combined_scaled.npy',
                           path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data',
                           path_out,
                           path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml',
                           n_clusters = 3):

    os.chdir(path_data_file)

    data = np.load(data_file)

    model = TimeSeriesKMeans(n_clusters= n_clusters, metric="dtw" )
    model.fit(data)

    os.chdir(path_fig)

    plt.hist(x=list(model.labels_))

    plt.xlabel("DTW K-means clusters")

    os.chdir(path_fig)

    plt.savefig("hist-clusters-drug-eff.png")

    for cluster_id in range(0,max(model.labels_)):

        idx = model.labels_ == cluster_id

        data_clustered = data[np.array(idx),]

        for i in random.sample(range(0, data_clustered.shape[0]), 10):
            plt.plot(data_clustered[i])
        plt.savefig('cluster'+str(cluster_id)+'_kmeans.png')
        plt.show()




