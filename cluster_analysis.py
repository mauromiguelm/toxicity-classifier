import os, random, json
import numpy as np
from tslearn.clustering import TimeSeriesKMeans
import matplotlib.pyplot as plt
import seaborn as sns

def run_clustering_methods(data_file = 'dist_combined.npy',
                           path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data',
                           path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml',
                           n_clusters = 4):
    "run clustering method on temporal distance files, and output cluster labels and a few diagnostic plots"
    os.chdir(path_data_file)

    data = np.load(data_file)

    model = TimeSeriesKMeans(n_clusters= n_clusters, metric="dtw" )

    model.fit(data)

    os.chdir(path_fig)

    plt.hist(x=list(model.labels_))

    plt.xlabel("DTW K-means clusters")

    plt.savefig("hist-clusters-drug-eff.png")

    plt.show()

    for cluster_id in range(0,max(model.labels_+1)):
        print(cluster_id)

        idx = model.labels_ == cluster_id

        data_clustered = data[np.array(idx),]

        for i in random.sample(range(0, data_clustered.shape[0]), 10):
            plt.plot(data_clustered[i])
        plt.savefig('cluster'+str(cluster_id)+'_kmeans.png')
        plt.show()

    os.chdir(path_data_file)

    np.save('model_cluster_labels', model.labels_)

    return model.labels_

def drug_centric_analysis(path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data',
                          metadata_file = 'metadata-pheno-ml.json',
                          cluster_labels_file = "model_cluster_labels.npy",
                          path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml'
                          ):
    "run drug-centric analysis, to observe possible differences in drug effect from clustering analysis"

    os.chdir(path_data_file)

    cluster_labels = np.load(cluster_labels_file)

    n_clusters = max(cluster_labels)

    with open(metadata_file) as output_file:
        metadata = json.load(output_file)


    drug_set = set(metadata['drug'])

    drugs = metadata['drug']

    clusters_by_drug = np.empty(shape=(0,n_clusters+1))

    for drug in drug_set:

        idx = [x == drug for x in drugs]

        drug_labels = cluster_labels[idx]

        cluster_freq = []

        for cluster in range(0,n_clusters+1):

            cluster_freq.append(sum(drug_labels == cluster)/len(drug_labels))

        cluster_freq = np.array(cluster_freq).reshape(1,n_clusters+1)

        clusters_by_drug = np.append(clusters_by_drug, cluster_freq, axis = 0)

    heat = sns.heatmap(clusters_by_drug, linewidth= 0.5, yticklabels = drug_set, center=0.3)

    os.chdir(path_fig)

    plt.savefig("cluster_drug_effect.png")

    plt.show()

if __name__ == "__main__":

    run_clustering_methods(data_file = 'dist_combined.npy',
                           path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data',
                           path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml',
                           n_clusters = 4)

    drug_centric_analysis(path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data',
                          metadata_file = 'metadata-pheno-ml.json',
                          cluster_labels_file = "model_cluster_labels.npy",
                          path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml'
                          )