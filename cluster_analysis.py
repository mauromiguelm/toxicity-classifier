import os, random, json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tslearn.clustering import TimeSeriesKMeans
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from sklearn.metrics import davies_bouldin_score

random.seed(10)

def run_clustering_methods(data,
                           n_clusters,
                           path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml',
                           path_out = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\cluster_pheno-ml'):
    "run clustering method on temporal distance files, and output cluster labels and a few diagnostic plots"

    model = TimeSeriesKMeans(n_clusters= n_clusters, metric="dtw" )

    model.fit(data)

    os.chdir(path_fig)

    plt.hist(x=list(model.labels_))

    plt.xlabel("DTW K-means clusters")

    plt.savefig("hist-clusters-drug-eff_scaled_denoise.png")

    plt.show()

    for cluster_id in range(0,max(model.labels_+1)):
        print(cluster_id)

        idx = model.labels_ == cluster_id

        data_clustered = data[np.array(idx),]

        for i in random.sample(range(0, data_clustered.shape[0]), 10):
            plt.plot(data_clustered[i])
        plt.savefig('nclus'+str(n_clusters)+'clusterID-'+str(cluster_id)+'_kmeans_scaled_denoise.png')
        plt.show()

    os.chdir(path_out)

    np.save('nclus'+str(n_clusters)+'model_cluster_labels_scaled_denoise', model.labels_)

    return(model.labels_)

def cluster_eval_metrics(X, labels, metric = 'euclidean'):
    'run evaluation metrics for different number of clusters'

    ss_metric = metrics.silhouette_score(X, labels, metric)

    ch_metric = metrics.calinski_harabasz_score(X, labels)

    db_metric = davies_bouldin_score(X, labels)

    return([ss_metric, ch_metric, db_metric])

def drug_centric_analysis(metadata,
                          cluster_labels,
                          path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml',
                          ):
    "run drug-centric analysis, to observe possible differences in drug effect from clustering analysis"

    n_clusters = max(cluster_labels)

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

    plt.savefig('nclus_'+str(n_clusters)+"_drug_effect_scaled_denoise.png")

    plt.show()

if __name__ == "__main__":

    path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data'
    path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml'
    path_out = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\cluster_pheno-ml'

    data_file = 'dist_combined_scaled_denoise.npy'

    os.chdir(path_data_file)

    data = np.load(data_file)

    metadata_file = 'metadata-pheno-ml.json'

    summary_eval_metrics = np.empty(shape=(0, 3))

    list_nclus = []

    with open(metadata_file) as output_file:
        metadata = json.load(output_file)

    for idx in range(2,10):
        'iterate for different number of clusters'

        list_nclus.append(idx)

        labels = run_clustering_methods(data = data,
                                        path_fig = path_fig,
                                        path_out = path_out,
                                        n_clusters = idx)

        metric = cluster_eval_metrics(X = data,
                                       labels = labels,
                                       metric = 'euclidean')

        metric = np.array(metric).reshape(1,3)


        summary_eval_metrics = np.append(summary_eval_metrics, metric, axis = 0)

        drug_centric_analysis(metadata = metadata,
                              cluster_labels = labels,
                              path_fig = path_fig)

    os.chdir(path_out)

    np.save('eval_metric', summary_eval_metrics)








