import os, random, json
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from PIL import Image
from tslearn.clustering import TimeSeriesKMeans
from scipy.signal import savgol_filter
from sklearn import metrics
from random import sample
from math import ceil

seed = 10
random.seed(seed)

def run_clustering_methods(data,
                           n_clusters,
                           path_fig,
                           path_out,
                           hist_plot,
                           cluster_plot,
                           ):
    "run clustering method on temporal distance files, and output cluster labels and a few diagnostic plots"

    model = TimeSeriesKMeans(n_clusters= n_clusters,
                             metric="dtw",
                             random_state=seed)

    model.fit(data)

    os.chdir(path_fig)

    ax = sns.histplot(data= model.labels_,
                      kde=True,
                      discrete = True
                      )

    ax.set(xlabel='DTW K-means clusters={}'.format(str(n_clusters)))

    plt.savefig("hist-" + hist_plot + 'cluster_n-' + str(n_clusters) + ".svg",
                transparent = True, dpi = 1200)

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

    plt.savefig('nclus' + str(n_clusters) + cluster_plot + '.svg')

    plt.close("all")

    os.chdir(path_out)

    np.save('labels_nclus_' + str(n_clusters), model.labels_)



    return(model.labels_)

def cluster_eval_metrics(X,
                         labels,
                         metric = 'euclidean'):
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

    plt.savefig('cluster_eval_metrics.svg')

    plt.close("all")

def drug_centric_analysis(metadata,
                          cluster_labels,
                          path_fig,
                          heatmap_label
                          ):
    "run drug-centric analysis, to observe possible differences in drug effect from clustering analysis"

    n_clusters = max(cluster_labels)

    drug_set = set(metadata['drug'])

    drugs = metadata['drug']

    clusters_by_drug = np.empty(shape=(0, n_clusters + 1))

    for drug in drug_set:
        if drug != "PBS":
            # drug = drug_set[1]
            idx = [x == drug for x in drugs]

            idx = np.array(idx, dtype='bool')

            drug_labels = np.array(cluster_labels)[idx]

            cluster_freq = []

            for cluster in range(0, n_clusters + 1):

                cluster_freq.append(sum(drug_labels == cluster) / len(drug_labels))

            cluster_freq = np.array(cluster_freq).reshape(1, n_clusters + 1)

            clusters_by_drug = np.append(clusters_by_drug, cluster_freq, axis=0)
        else:
            pass

    drug_names = [x for x in drug_set if x != "PBS"]

    heat = sns.heatmap(data = clusters_by_drug,
                       linewidth=0.5,
                       yticklabels=drug_names,
                       cmap="YlOrBr")

    plt.tight_layout()

    os.chdir(path_fig)

    plt.savefig('nclus_'+str(n_clusters+1)+heatmap_label+".svg", dpi=1200)

    plt.close("all")

def cell_freq_by_cluster(metadata,
                         cluster_labels,
                         path_fig,
                         heatmap_label
                         ):
    "run drug-centric analysis, to observe possible differences in drug effect from clustering analysis"

    n_clusters = max(cluster_labels)

    cell_set = set(metadata['cell'])

    cells = metadata['cell']

    clusters_by_cell = np.empty(shape=(0, n_clusters + 1))

    for cell in cell_set:

        idx = [x == cell for x in cells]

        idx = np.array(idx, dtype='bool')

        cell_labels = np.array(cluster_labels)[idx]

        cluster_freq = []

        for cluster in range(0, n_clusters + 1):

            cluster_freq.append(sum(cell_labels == cluster) / len(cell_labels))

        cluster_freq = np.array(cluster_freq).reshape(1, n_clusters + 1)

        clusters_by_cell = np.append(clusters_by_cell, cluster_freq, axis=0)

    cell_names = [x for x in cell_set]

    heat = sns.heatmap(data = clusters_by_cell,
                       linewidth=0.5,
                       yticklabels=cell_names,
                       cmap="YlOrBr")

    plt.tight_layout()

    os.chdir(path_fig)

    plt.savefig('nclus_'+str(n_clusters+1)+heatmap_label+".svg", dpi=1200)

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

        plt.savefig(cell+"_"+heatmap_label+".svg")

        #plt.show()
        plt.close("all")

        del new_path

def max_conc_analysis(metadata,
                       cluster_labels,
                       path_fig,
                       heatmap_label):
    drugs = metadata['drug']
    cells = metadata['cell']
    concs = metadata['conc']

    iter_cells = set(cells)
    iter_drugs = set(drugs)

    comb_array = np.empty(shape = (len(iter_cells), len(iter_drugs)))
    comb_array = pd.DataFrame(comb_array, columns=iter_drugs, index=iter_cells)


    for cell in iter_cells:
        "generate cell-centric results for multiple drugs and their concentrations"
        #cell = 'SKMEL2'

        idx = [x == cell for x in cells]
        idx = np.array(idx, dtype='bool')

        drug_sub = np.array(drugs)[idx]
        conc_sub = np.array(concs, dtype='float')[idx]
        label_sub = np.array(cluster_labels)[idx]

        for drug in iter_drugs:
            if drug != 'PBS':
                #drug = "Methotrexate"
                idx_drug = [x == drug for x in drug_sub]
                idx_drug = np.array(idx_drug, dtype='bool')

                conc_of_drug = conc_sub[idx_drug]
                label_of_drug = label_sub[idx_drug]

                max_conc = max(conc_of_drug)

                label_at_max_conc = label_of_drug[conc_of_drug == max_conc]

                comb_array[drug][cell] = label_at_max_conc

    comb_array = comb_array.drop("PBS", axis = 1)

    os.chdir(path_fig)

    my_colors = ['xkcd:grey', 'xkcd:orange', 'xkcd:apple green','xkcd:red']

    drug_names = [x for x in iter_drugs if x != "PBS"]

    heat = sns.clustermap(comb_array,
                          yticklabels=iter_cells,
                          xticklabels=drug_names,
                          cmap = my_colors,
                          metric = 'correlation'
                          )

    heat.cax.set_visible(False)

    plt.tight_layout()

    plt.savefig(heatmap_label+".svg", transparent = True, dpi = 1200)

def drug_conc_centric_analysis(metadata,
                               cluster_labels,
                               path_fig):

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
                             colors=['xkcd:grey', 'xkcd:orange', 'xkcd:apple green',
                                     'xkcd:red']
                             )

                    conc_iter = cluster_freq[:, 3].reshape(1, n_clusters + 1)

                    axes.pie(conc_iter.sum(axis=0), radius=0.85,
                             wedgeprops=dict(width=0.3, edgecolor='w'), normalize=True,
                             colors=['xkcd:grey', 'xkcd:orange', 'xkcd:apple green',
                                     'xkcd:red']
                             )

                    conc_iter = cluster_freq[:, 2].reshape(1, n_clusters + 1)

                    axes.pie(conc_iter.sum(axis=0), radius=0.70,
                             wedgeprops=dict(width=0.3, edgecolor='w'), normalize=True,
                             colors=['xkcd:grey', 'xkcd:orange', 'xkcd:apple green',
                                     'xkcd:red']
                             )

                    conc_iter = cluster_freq[:, 1].reshape(1, n_clusters + 1)

                    axes.pie(conc_iter.sum(axis=0), radius=0.55,
                             wedgeprops=dict(width=0.3, edgecolor='w'), normalize=True,
                             colors=['xkcd:grey', 'xkcd:orange', 'xkcd:apple green',
                                     'xkcd:red']
                             )

                    conc_iter = cluster_freq[:, 0].reshape(1, n_clusters + 1)

                    axes.pie(conc_iter.sum(axis=0), radius=0.40,
                             wedgeprops=dict(width=0.3, edgecolor='w'), normalize=True,
                             colors=['xkcd:grey', 'xkcd:orange', 'xkcd:apple green',
                                     'xkcd:red']
                             )

                else:
                    pass
                    # axes.set_title(str(drug).format(i, j))
                    # axes.set_yticklabels([])
                    # axes.set_xticklabels([])
                    #
                    # plt.plot()

            else:
                pass
                # axes.set_title(str(drug).format(i, j))
                # axes.set_yticklabels([])
                # axes.set_xticklabels([])

                # plt.plot()


            count += 1

    plt.tight_layout()
    plt.savefig("drug-conc_pie-plots.svg", transparent = True, dpi = 1200)
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

        plt.savefig("heatmap.svg")

        #plt.show()
        plt.close("all")

        del new_path


def biological_inference_of_clusters(chosen_cluster ,
                                     path_fig
                                     ):

        labels_eff = chosen_cluster

        noEffect            =   [0,4, 5, 7] #as 0
        mixedEffect         =   [2,6]       #as 1
        cytostatic          =   [3]         #as 2
        cytotoxic           =   [1]         #as 3

        labels_eff = ["No Effect" if x in set(noEffect) else x for x in labels_eff]
        labels_eff = ["Mixed Effect" if x in set(mixedEffect) else x for x in labels_eff]
        labels_eff = ["Cytostatic" if x in set(cytostatic) else x for x in labels_eff]
        labels_eff = ["Cytotoxic" if x in set(cytotoxic) else x for x in labels_eff]

        labels_eff = [0 if x == "No Effect" else x for x in labels_eff]
        labels_eff = [1 if x == "Mixed Effect" else x for x in labels_eff]
        labels_eff = [2 if x == "Cytostatic" else x for x in labels_eff]
        labels_eff = [3 if x == "Cytotoxic" else x for x in labels_eff]

        drug_centric_analysis(metadata = metadata,
                              cluster_labels = labels_eff,
                              path_fig = path_fig,
                              heatmap_label="_drug_effect2_cell-drug_scaled_denoise_eff")
        return(labels_eff)

def get_drug_position(path_data_file):
    os.chdir(path_data_file)

    MSP1 = pd.read_excel('randomized_layout_1MSP_batch2.xls')
    MSP1['MSP_plate'] = "P1"
    MSP2 = pd.read_excel('randomized_layout_2MSP_batch2.xls')
    MSP2['MSP_plate'] = "P2"

    MSP_full = [MSP1, MSP2]

    MSP_full = pd.concat(MSP_full)

    return(MSP_full)


def get_raw_images(metadata,
                   cluster_labels,
                   output_path,
                   n_images_per_cluster
                   ):

    drug_map = get_drug_position('\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data')
    drugs = metadata['drug']
    cells = metadata['cell']
    concs = metadata['conc']

    cells_image_dir = []

    for n in range(1,8):
        batch = "batch_{}/".format(str(n))
        batch_in = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users/Mauro/Cell_culture_data/190310_LargeScreen/imageData/" + batch

        cell_in_batch_dir = [x.path for x in os.scandir(batch_in) if x.is_dir()]

        [cells_image_dir.append(x) for x in cell_in_batch_dir]

    for cluster in range(0,max(cluster_labels)+1):
        #cluster = 2
        #subset list of clusters from metadata

        output_path_cluster = output_path+"/cluster_{}/".format(str(cluster))

        if not os.path.exists(output_path_cluster):
            os.makedirs(output_path_cluster)
        else:
            pass

        idx = [x == cluster for x in labels]
        idx = np.array(idx, dtype='bool')

        drug_sub = np.array(drugs)[idx]
        cell_sub = np.array(cells)[idx]
        conc_sub = np.array(concs)[idx]

        #take 15 random samples from the clusters

        idx_random = sample(range(0,len(drug_sub)), n_images_per_cluster)

        drug_rand = np.array(drug_sub)[idx_random]
        cell_rand = np.array(cell_sub)[idx_random]
        conc_rand = np.array(conc_sub)[idx_random]

        for idx in range(0,len(drug_rand)):
            #idx = 1

            drug_to_map = drug_rand[idx]
            cell_to_map = cell_rand[idx]
            conc_to_map =  conc_rand[idx]

            #convert concentration into well (available in drug map)

            d1 = np.array(drug_map["Drug"] == drug_to_map, dtype='bool')
            d2 = np.array(drug_map["Final_conc_uM"].round(5) == float(conc_to_map), dtype = 'bool')

            well_to_search = drug_map[d1 & d2 ]

            plate_to_search = well_to_search["MSP_plate"].unique()

            batch_matched = re.compile(".*({}).*".format(cell_to_map))

            if len(plate_to_search) > 1:
                exit()
            else:
                pass


            cell_folder = [x for x in cells_image_dir if batch_matched.search(x)]

            if plate_to_search== "P1":
                match_pattern = ".*({}).*".format("P1")
                plate_to_match = re.compile(match_pattern)
                plate_to_match = [x for x in cell_folder if plate_to_match.search(x)]
            else:
                pass

            if plate_to_search == "P2":
                match_pattern = ".*({}).*".format("P2")
                plate_to_match = re.compile(match_pattern)
                plate_to_match = [x for x in cell_folder if plate_to_match.search(x)]
            else:
                pass

            os.chdir(plate_to_match[0])

            well_position_list = []

            for x in range(0,len(well_to_search)):
                #x = 0

                col = well_to_search['Column'].to_list()[x]
                row = well_to_search.iloc[x]['Row']
                well_position = row + str(col)

                well_position_list.append(well_position)

            well_position_list = ["_{}_".format(x) for x in well_position_list]

            images_to_move = []

            for x in os.listdir():

                if bool(re.search("|".join(well_position_list), x)) == True:
                    images_to_move.append(x)
                else:
                    pass

            for image in images_to_move:

                os.chdir(plate_to_match[0])

                image_file = Image.open(image)

                os.chdir(output_path_cluster)

                image_file.save(image)

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
        #idx = 8

        list_nclus.append(idx)

        labels = run_clustering_methods(data = data,
                                        path_fig = path_fig,
                                        path_out = path_out,
                                        n_clusters = idx,
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

        cell_freq_by_cluster(metadata = metadata,
                             cluster_labels = labels,
                             path_fig = path_fig,
                             heatmap_label = "freq-by-cell"
                             )

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
    os.chdir(path_out)

    labels = np.load("labels_nclus_8.npy")

    get_raw_images(metadata=metadata,
                   cluster_labels=labels,
                   output_path='\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml\\images_by_cluster',
                   n_images_per_cluster=3
                   )

    labels_eff = biological_inference_of_clusters(chosen_cluster = labels,
                                                  path_fig = path_fig)

    os.chdir(path_out)

    np.save('labels_nclus_bio-effect', labels_eff)

    get_raw_images(metadata = metadata,
                   cluster_labels = labels_eff,
                   output_path='\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\pheno-ml\\images_by_cluster_eff',
                   n_images_per_cluster=3
                   )

    drug_centric_analysis(metadata=metadata,
                          cluster_labels=labels_eff,
                          path_fig=path_fig,
                          heatmap_label="_drug_biological_labels")

    drug_conc_centric_analysis(metadata = metadata,
                               cluster_labels=labels_eff,
                               path_fig = path_fig)

    max_conc_analysis(metadata = metadata,
                      cluster_labels = labels_eff,
                      path_fig = path_fig,
                      heatmap_label= "max_conc_gram")




