import os, random, json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tslearn.clustering import TimeSeriesKMeans
from sklearn import metrics
from sklearn.metrics import davies_bouldin_score

random.seed(10)

def run_clustering_methods(data,
                           n_clusters,
                           path_fig,
                           path_out):
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
                          path_fig,
                          ):
    "run drug-centric analysis, to observe possible differences in drug effect from clustering analysis"

    #TODO Focus on n_cluster = 9 and foccus on the cluster#2, since it shows
    #TODO cytostatic plus cytotoxic.. check if its a concentration dependent
    #TODO or cell line dependent

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

def cell_centric_analysis(metadata,
                          cluster_labels,
                          path_fig):
    drugs = metadata['drug']
    cells = metadata['cell']
    concs = metadata['conc']

    for cell in set(cells):
        "generate cell-centric results for multiple drugs and their concentrations"
        #cell = 'SKMEL2'

        idx = [x == cell for x in cells]
        idx = np.array(idx, dtype='bool')

        drug_sub = np.array(drugs)[idx]
        conc_sub = np.array(concs, dtype='float')[idx]
        label_sub = np.array(cluster_labels)[idx]

        drug_label = np.empty(shape=(0,5))

        drug_name = []

        for drug in set(drug_sub):
            if drug != 'PBS':
                # drug = drug_sub.item(0)
                idx_drug = [x == drug for x in drug_sub]
                idx_drug = np.array(idx_drug, dtype='bool')

                conc_drug = conc_sub[idx_drug]

                if len(conc_drug) == 5:

                    label_drug = label_sub[idx_drug]
                    label_drug = label_drug.reshape(1, 5)

                    drug_label = np.append(drug_label, label_drug, axis=0)

                    drug_name.append(drug)

                else:
                    pass
            else:
                pass

        print('I am working'+cell)

        new_path = path_fig+'\\'+cell

        if os.path.exists(new_path) == False:
            os.makedirs(new_path)
        else:
            pass

        os.chdir(new_path)

        print(new_path)

        heat = sns.heatmap(drug_label,
                           yticklabels = drug_name
                           )


        plt.savefig("heatmap.png")

        plt.show()

        del new_path
















    #TODO plot clusters for each drug, showing cell line specific differences
    #TODO check differences across concentration, to see if some are toxic and if some are citostatic


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

    os.chdir(path_fig)

    # Create a subplot with 1 rows and 3 columns
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.set_size_inches(21, 7)

    ax1.plot(list_nclus, summary_eval_metrics[:,0])

    ax1.set_title("The silhouette score for the various clusters.")
    ax1.set_ylabel("Silhouette score")
    ax1.set_xlabel("Cluster label")

    ax2.plot(list_nclus, summary_eval_metrics[:, 1])

    ax2.set_title("The Calinski-Harabasz score for the various clusters.")
    ax2.set_ylabel("Calinski-Harabasz score")
    ax2.set_xlabel("Cluster label")

    ax3.plot(list_nclus, summary_eval_metrics[:, 2])

    ax3.set_title("The Davies-Bouldin score for the various clusters.")
    ax3.set_ylabel("Davies-Bouldin score")
    ax3.set_xlabel("Cluster label")

    os.chdir(path_fig)

    plt.savefig('cluster_eval_metrics.png')

    plt.show()

    cell_centric_analysis(metadata = metadata,
                          cluster_labels = labels,
                          path_fig = path_fig)



