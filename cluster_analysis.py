import os, random, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tslearn.clustering import TimeSeriesKMeans
from scipy.signal import savgol_filter
from sklearn import metrics
from math import ceil

seed = 10
random.seed(seed)

def run_clustering_methods(data,
                           n_clusters,
                           path_fig,
                           path_out,
                           output_file,
                           hist_plot,
                           cluster_plot,
                           ):
    "run clustering method on temporal distance files, and output cluster labels and a few diagnostic plots"

    model = TimeSeriesKMeans(n_clusters= n_clusters,
                             metric="dtw",
                             random_state=seed)

    model.fit(data)

    os.chdir(path_fig)

    plt.hist(x=list(model.labels_))

    plt.xlabel("DTW K-means clusters")

    plt.savefig("hist-"+hist_plot+".png")

    #plt.show()
    plt.close("all")

    plt.figure()
    sz = data.shape[1]
    for cluster_id in range(0, max(model.labels_ + 1)):

        idx = model.labels_ == cluster_id

        data_clustered = data[np.array(idx),]

        plt.subplot(3, 3, cluster_id + 1)
        for xx in data_clustered:
            plt.plot(xx.ravel(), "k-", alpha=.2)

        plt.plot(savgol_filter(model.cluster_centers_[cluster_id].ravel(), 7, 2), "r-", linewidth = 2.5)
        plt.xlim(0, sz)
        plt.ylim(0, 1.2)
        plt.text(0.55, 0.85, 'Cluster %d' % (cluster_id),
                 transform=plt.gca().transAxes)

    plt.tight_layout()

    plt.savefig('nclus' + str(n_clusters) + 'clusterID-' + str(cluster_id) + cluster_plot + '.png')

    plt.close("all")

    os.chdir(path_out)

    np.save('labels_nclus_' + str(idx), model.labels_)

    return(model.labels_)

def cluster_eval_metrics(X, labels, metric = 'euclidean'):
    'run evaluation metrics for different number of clusters'

    ss_metric = metrics.silhouette_score(X, labels, metric)

    ch_metric = metrics.calinski_harabasz_score(X, labels)

    db_metric = metrics.davies_bouldin_score(X, labels)

    return([ss_metric, ch_metric, db_metric])

def plot_eval_metrics(list_nclus,
                      summary_eval_metrics,
                      path_fig):
    # Create a subplot with 1 rows and 3 columns for visualize cluster evaluation metrics

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.set_size_inches(21, 7)

    ax1.plot(list_nclus, summary_eval_metrics[:, 0])

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

    plt.close("all")

def drug_centric_analysis(metadata,
                          cluster_labels,
                          path_fig,
                          heatmap_label
                          ):
    "run drug-centric analysis, to observe possible differences in drug effect from clustering analysis"

    #TODO add pie chart with summary for clusters
    #TODO pie chart should have max concentration and then first effect
    n_clusters = max(cluster_labels)

    drug_set = list(set(metadata['drug']))

    drugs = metadata['drug']

    clusters_by_drug = np.empty(shape=(0,n_clusters+1))

    pie_size = 1

    ncol = 2
    nrow = ceil(len(drug_set) / ncol)


    fig, ax = plt.subplots(nrow, ncol, figsize=(10, 40))

    count = 0

    for i, ax_row in enumerate(ax):
        for j, axes in enumerate(ax_row):

            drug = drug_set[count]

            idx = [x == drug for x in drugs]
            idx = np.array(idx, dtype='bool')

            drug_labels = np.array(cluster_labels)[idx]

            cluster_freq = []

            for cluster in range(0,n_clusters+1):

                cluster_freq.append(sum(drug_labels == cluster)/len(drug_labels))

            cluster_freq = np.array(cluster_freq).reshape(1,n_clusters+1)

            axes.set_title(str(drug).format(i,j))
            axes.set_yticklabels([])
            axes.set_xticklabels([])

            axes.pie(cluster_freq.sum(axis=0), radius=1,
                      wedgeprops=dict(width=0.3,edgecolor='w'), normalize=True,
                     colors=['xkcd:grey', 'xkcd:apple green', 'xkcd:orange',
                             'xkcd:red']
                     )
            count += 1

    fig.subplots_adjust(wspace=.2)

    #plt.show()
    plt.close("all")

    clusters_by_drug = np.append(clusters_by_drug, cluster_freq, axis = 0)

    heat = sns.heatmap(clusters_by_drug, linewidth= 0.5, yticklabels = drug_set, center=0.3)

    os.chdir(path_fig)

    plt.savefig('nclus_'+str(n_clusters)+heatmap_label+".png")

    #plt.show()
    plt.close("all")

def cell_centric_analysis(metadata,
                          cluster_labels,
                          path_fig,
                          heatmap_label):
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
                #drug = drug_sub.item(0)
                idx_drug = [x == drug for x in drug_sub]
                idx_drug = np.array(idx_drug, dtype='bool')

                conc_drug = conc_sub[idx_drug]

                if len(conc_drug) == 5:

                    label_drug = label_sub[idx_drug]

                    conc_sort = conc_drug.argsort()

                    label_drug = label_drug[conc_sort[::-1]]

                    label_drug = label_drug.reshape(1, 5)

                    drug_label = np.append(drug_label, label_drug, axis=0)

                    drug_name.append(drug)

                else:
                    pass
            else:
                pass

        new_path = path_fig+'\\'+'cell-analysis'

        if os.path.exists(new_path) == False:
            os.makedirs(new_path)
        else:
            pass

        os.chdir(new_path)

        #print(new_path)

        heat = sns.heatmap(drug_label,
                           yticklabels = drug_name,
                           annot=True
                           )

        plt.savefig(cell+"_"+heatmap_label+".png")

        #plt.show()
        plt.close("all")

        del new_path


def drug_conc_centric_analysis(metadata,
                               cluster_labels,
                               path_fig,
                               heatmap_label):

    "run drug-centric analysis, to observe possible differences in drug effect from clustering analysis"

    n_clusters = max(cluster_labels)

    drug_set = list(set(metadata['drug']))

    drugs = metadata['drug']
    concs = metadata['conc']

    clusters_by_drug = np.empty(shape=(0,n_clusters+1))


    os.chdir(path_fig)

    pie_size = 1

    ncol = 2
    nrow = ceil(len(drug_set) / ncol)


    fig, ax = plt.subplots(nrow, ncol, figsize=(10, 40))

    count = 0

    for i, ax_row in enumerate(ax):
        for j, axes in enumerate(ax_row):

            drug = drug_set[count]


            idx = [x == drug for x in drugs]
            idx = np.array(idx, dtype='bool')

            conc_sub = np.array(metadata['conc'])[idx]
            drug_labels = np.array(cluster_labels)[idx]

            cluster_freq = np.empty(shape=(n_clusters+1, 5))

            if drug != 'PBS':
                if len(set(conc_sub)) == 5:

                    for cluster in range(0,n_clusters+1):

                        #cluster = 0
                        idx_cluster = [x == cluster for x in drug_labels]
                        idx_cluster = np.array(idx_cluster, dtype='bool')

                        conc_labels_sub = conc_sub[idx_cluster]
                        count_conc = 0
                        for conc in sorted(set(conc_sub), reverse=False):
                            #conc = conc_sub[1]

                            cluster_freq[cluster,count_conc] = sum(conc_labels_sub == conc) / sum(idx_cluster)

                            count_conc += 1

                    axes.set_title(str(drug).format(i, j))
                    axes.set_yticklabels([])
                    axes.set_xticklabels([])

                    conc_iter = cluster_freq[:, 4].reshape(1, n_clusters + 1)

                    axes.pie(conc_iter.sum(axis=0), radius=1,
                             wedgeprops=dict(width=0.3, edgecolor='w'), normalize=True,
                             colors=['xkcd:grey', 'xkcd:apple green', 'xkcd:orange',
                                     'xkcd:red']
                             )

                    conc_iter = cluster_freq[:, 3].reshape(1, n_clusters + 1)

                    axes.pie(conc_iter.sum(axis=0), radius=0.85,
                             wedgeprops=dict(width=0.3, edgecolor='w'), normalize=True,
                             colors=['xkcd:grey', 'xkcd:apple green', 'xkcd:orange',
                                     'xkcd:red']
                             )

                    conc_iter = cluster_freq[:, 2].reshape(1, n_clusters + 1)

                    axes.pie(conc_iter.sum(axis=0), radius=0.70,
                             wedgeprops=dict(width=0.3, edgecolor='w'), normalize=True,
                             colors=['xkcd:grey', 'xkcd:apple green', 'xkcd:orange',
                                     'xkcd:red']
                             )

                    conc_iter = cluster_freq[:, 1].reshape(1, n_clusters + 1)

                    axes.pie(conc_iter.sum(axis=0), radius=0.55,
                             wedgeprops=dict(width=0.3, edgecolor='w'), normalize=True,
                             colors=['xkcd:grey', 'xkcd:apple green', 'xkcd:orange',
                                     'xkcd:red']
                             )

                    conc_iter = cluster_freq[:, 0].reshape(1, n_clusters + 1)

                    axes.pie(conc_iter.sum(axis=0), radius=0.40,
                             wedgeprops=dict(width=0.3, edgecolor='w'), normalize=True,
                             colors=['xkcd:grey', 'xkcd:apple green', 'xkcd:orange',
                                     'xkcd:red']
                             )

                else:
                    axes.set_title(str(drug).format(i, j))
                    axes.set_yticklabels([])
                    axes.set_xticklabels([])

                    plt.plot()

            else:
                axes.set_title(str(drug).format(i, j))
                axes.set_yticklabels([])
                axes.set_xticklabels([])

                plt.plot()

            count += 1

    plt.tight_layout()
    plt.savefig("drug-conc_pie-plots.png", transparent = True, dpi = 1200)
    #plt.show()
    plt.close("all")

    clusters_by_drug = np.append(clusters_by_drug, cluster_freq, axis = 0)

    heat = sns.heatmap(clusters_by_drug, linewidth= 0.5, yticklabels = drug_set, center=0.3)

    os.chdir(path_fig)

    plt.savefig('nclus_'+str(n_clusters)+heatmap_label+".png")

    #plt.show()
    plt.close("all")

def conc_centric_analysis(metadata,
                          cluster_labels,
                          path_fig):
    drugs = metadata['drug']
    cells = metadata['cell']
    concs = metadata['conc']

    conc_summary_by_label = np.zeros(shape=(5, len(set(cluster_labels))))

    for drug in set(drugs):
        "generate concentration-centric results"
        #drug = drugs[100]

        idx = [x == drug for x in drugs]
        idx = np.array(idx, dtype='bool')

        conc_sub = np.array(concs, dtype='float')[idx]
        label_sub = np.array(cluster_labels)[idx]

        drug_label = np.empty(shape=(0,5))

        drug_name = []

        if drug != 'PBS':

            if len(set(conc_sub)) == 5:

                iter_conc = 0

                for conc in sorted(set(conc_sub), reverse=True):
                    #conc = sorted(set(conc_sub), reverse=True)[3]

                    idx_conc = [x == conc for x in conc_sub]
                    idx_conc = np.array(idx_conc, dtype='bool')

                    label_for_conc = label_sub[idx_conc]

                    label_per_concentration = [0]*len(set(cluster_labels))

                    for label in sorted(set(label_sub), reverse=False):
                        label_per_concentration[label] = label_per_concentration[label] + sum(label_for_conc == label)

                    conc_summary_by_label[iter_conc,:] = conc_summary_by_label[iter_conc,:] +label_per_concentration

                    iter_conc = iter_conc + 1

                    label_drug = label_sub[idx_conc]

                    label_drug = label_drug.reshape(1, 5)

                    drug_label = np.append(drug_label, label_drug, axis=0)

                    drug_name.append(drug)
            else:
                pass
        else:
            pass

        new_path = path_fig+'\\'+ drug

        if os.path.exists(new_path) == False:
            os.makedirs(new_path)
        else:
            pass

        os.chdir(new_path)

        heat = sns.heatmap(drug_label,
                           yticklabels = drug_name
                           )

        plt.savefig("heatmap.png")

        #plt.show()
        plt.close("all")

        del new_path


def biological_inference_of_clusters(chosen_cluster ,
                                     path_data_file,
                                     path_fig
                                     ):
        os.chdir(path_data_file)

        labels_eff = np.load(chosen_cluster)

        noEffect            =   [0, 2, 5, 6, 7] #as 0
        cytostatic          =   [3]             #as 1
        weak_cytotoxic      =   [4, 8]          #as 2
        strong_cytotoxic    =   [1]             #as 3

        labels_eff = ["NoEff" if x in set(noEffect) else x for x in labels_eff]
        labels_eff = ["Cytostatic" if x in set(cytostatic) else x for x in labels_eff]
        labels_eff = ["Weak Cytotoxic" if x in set(weak_cytotoxic) else x for x in labels_eff]
        labels_eff = ["Strong Cytotoxic" if x in set(strong_cytotoxic) else x for x in labels_eff]

        labels_eff = [0 if x == "NoEff" else x for x in labels_eff]
        labels_eff = [1 if x == "Cytostatic" else x for x in labels_eff]
        labels_eff = [2 if x == "Weak Cytotoxic" else x for x in labels_eff]
        labels_eff = [3 if x == "Strong Cytotoxic" else x for x in labels_eff]


        drug_centric_analysis(metadata = metadata,
                              cluster_labels = labels_eff,
                              path_fig = path_fig,
                              heatmap_label="_drug_effect2_cell-drug_scaled_denoise_eff")


if __name__ == "__main__":

    path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data'
    path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml'
    path_out = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\cluster_pheno-ml'

    data_file = 'dist_combined_celldrug-scaled_denoise.npy'

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
                                        n_clusters = idx,
                                        output_file = "model_cluster_labels_celldrug-scaled_denoise",
                                        hist_plot = "clusters-celldrug-scaled_denoise",
                                        cluster_plot = '_kmeans_cell-drugscaled_denoise_avg'
                                        )

        metric = cluster_eval_metrics(X = data,
                                       labels = labels,
                                       metric = 'euclidean')

        metric = np.array(metric).reshape(1,3)

        summary_eval_metrics = np.append(summary_eval_metrics, metric, axis = 0)

        drug_centric_analysis(metadata = metadata,
                              cluster_labels = labels,
                              path_fig = path_fig,
                              heatmap_label="_drug_effect2_cell-drug_scaled_denoise")

        cell_centric_analysis(metadata=metadata,
                              cluster_labels=labels,
                              path_fig=path_fig,
                              heatmap_label="-heatmap_cell-drug_scaled2")

        #TODO fix conc_centric_analysis

        # conc_centric_analysis(metadata=metadata,
        #                       cluster_labels=labels,
        #                       path_fig=path_fig)

    np.save('eval_metric', summary_eval_metrics)

    os.chdir(path_fig)

    plot_eval_metrics(list_nclus,
                      summary_eval_metrics,
                      path_fig)

    biological_inference_of_clusters(chosen_cluster = "labels_nclus_9.npy",
                                     path_data_file = path_out,
                                     path_fig
                                     )


