import itertools
import multiprocessing as mp
import time

import matplotlib.patches as mpatches
import numpy as np
import seaborn as sn
from Bio.SVDSuperimposer import SVDSuperimposer
from matplotlib import pyplot as plt
from pymol import cmd
from scipy import cluster
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score

from utilities import generate_colormap

n_best = 20
structure_coords = dict()
counter: mp.Value
window = None


def cm_to_inch(x):
    return x / 2.54


def hierarchy_optimization(X, logger, height, desired_clusters, method='complete'):
    result_label = dict()
    cut_heights = dict()

    range_n_clusters = list(range(2, len(X)))
    Z = cluster.hierarchy.linkage(squareform(X), optimal_ordering=True, method=method)

    centroid_clusters = []
    names = range(len(X))
    clusters = cluster.hierarchy.cut_tree(Z, n_clusters=range_n_clusters)
    
    if desired_clusters is not None:
        cluster_labels = clusters[:, desired_clusters - 2].flatten()

        silhouette_avg = silhouette_score(X, cluster_labels, metric='precomputed')
        result_label.setdefault(desired_clusters, (silhouette_avg, cluster_labels))
        tmp = cluster.hierarchy.dendrogram(Z, p=desired_clusters, truncate_mode='lastp', no_plot=True)
        last_h = min(list(filter(lambda x: x != 0, [item for sublist in tmp['dcoord'] for item in sublist])))
        cut_heights.setdefault(desired_clusters, last_h)

        for cluster_id in range(desired_clusters):
            nameList = list(zip(names, cluster_labels))
            mask = np.array([i == cluster_id for i in cluster_labels])
            idx = np.argmin(sum(X[:, mask][mask, :]))
            sublist = [name for (name, label) in nameList if label == cluster_id]
            centroid_clusters.append(sublist[idx])
    elif height is not None:
        for i, n_clusters in enumerate(range_n_clusters):
            logger.progress((i / len(range_n_clusters)) * 100)
            cluster_labels = clusters[:, i].flatten()

            silhouette_avg = silhouette_score(X, cluster_labels, metric='precomputed')
            result_label.setdefault(n_clusters, (silhouette_avg, cluster_labels))

            tmp = cluster.hierarchy.dendrogram(Z, p=n_clusters, truncate_mode='lastp', no_plot=True)
            last_h = min(list(filter(lambda x: x != 0, [item for sublist in tmp['dcoord'] for item in sublist])))
            cut_heights.setdefault(n_clusters, last_h)

            if last_h <= height or n_clusters == len(range_n_clusters) - 1:
                for cluster_id in range(n_clusters):
                    nameList = list(zip(names, cluster_labels))
                    mask = np.array([i == cluster_id for i in cluster_labels])
                    idx = np.argmin(sum(X[:, mask][mask, :]))
                    sublist = [name for (name, label) in nameList if label == cluster_id]
                    centroid_clusters.append(sublist[idx])
                break

    logger.close_progress()
    return result_label, Z, cut_heights, centroid_clusters


def cluster_distribution_heatmap(logger, pdb_id, method, rmsd_val=None, desired_clusters=None, x_len=50):
    logger.disable_window()

    X = load_rmsd_dis_matrix(logger, pdb_id)

    labels, Z, cut_heights, centroid_cluster = hierarchy_optimization(X, logger, height=rmsd_val,
                                                                      desired_clusters=desired_clusters, method=method)

    if desired_clusters is not None and desired_clusters not in labels:
        logger.log("The number of cluster has to be in the range 2 - {} (inclusive)".format(max(labels)))
        return
    else:
        logger.log('Number of clusters for selected RMSD cut: {}'.format(len(centroid_cluster)), warning=True)

    n_clusters = desired_clusters if desired_clusters is not None else len(centroid_cluster)

    labels = labels[n_clusters][1]

    def remove_spines(ax):
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['bottom'].set_color('none')

    plt.close()
    plt.style.use('default')
    fig, ax1 = plt.subplots(dpi=70)

    pad_size = x_len - len(labels) % x_len if len(labels) % x_len != 0 else 0
    labels = labels.astype('float32')
    best = np.pad(labels, (0, pad_size), mode='constant', constant_values=np.nan)
    best = np.reshape(best, (int(len(best) / x_len), x_len))

    cmap = generate_colormap(n_clusters)

    sn.set(font_scale=0.8)
    ax1 = sn.heatmap(best, cmap=cmap, linewidths=.5, linecolor='white', ax=ax1, square=True, cbar=True, annot=False)
    ax1.set_yticks(np.arange(0.5, best.shape[0] + 0.5))
    ax1.set_yticklabels(np.arange(1, best.size + 1, x_len), rotation=0, fontsize=10)
    ax1.set_xticks(np.arange(4.5, best.shape[1] + 0.5, 5))
    ax1.set_xticklabels(range(5, best.shape[1] + 1, 5))

    handles = []
    labels = []
    for i in range(n_clusters):
        handles.append(mpatches.Patch(color=cmap.colors[i]))
        labels.append(centroid_cluster[i] + 1)

    def flip(items, ncol):
        return itertools.chain(*[items[i::ncol] for i in range(ncol)])

    plt.legend(handles=flip(handles, 20), labels=flip(labels, 20), bbox_to_anchor=(0., 0., 1., -.2), loc='upper left',
               ncol=20, mode="expand", borderaxespad=0., facecolor="white", handlelength=1.2, handleheight=1.2,
               fontsize='large')

    remove_spines(ax1)
    logger.enable_window()

    plt.xlabel('States', fontdict={'size': 13, })
    plt.title("Structure {} - {} clusters".format(pdb_id, len(set(labels))),
              fontsize=14)
    plt.tight_layout()
    plt.grid(False)
    plt.show(block=False)


def hierarchy_cut_plot(logger, pdb_id, method, rmsd_val=None, desired_clusters=None):
    X = load_rmsd_dis_matrix(logger, pdb_id)

    result_labels, Z, cut_heights, centroids = hierarchy_optimization(X, logger, rmsd_val, desired_clusters,
                                                                      method=method)

    if desired_clusters is not None:
        if desired_clusters in result_labels:
            silh_val = result_labels[desired_clusters][0]
            y_val = cut_heights[desired_clusters]
        else:
            logger.log("The number of cluster has to be in the range 2 - {} (inclusive)".format(max(result_labels)))
            return
    else:
        logger.log('Number of clusters for selected RMSD cut: {}'.format(len(centroids)), warning=True)
        silh_val = result_labels[len(centroids)][0]
        y_val = cut_heights[len(centroids)]

    n_clusters = desired_clusters if desired_clusters is not None else len(centroids)
    plt.close()
    plt.style.use('default')
    plt.figure(1, figsize=(12, 9))
    plt.axhline(y=y_val, linestyle="--", zorder=0, linewidth=1.3,
                label="{} clusters, silh: {:.3f}".format(n_clusters, silh_val))

    R = cluster.hierarchy.dendrogram(Z, no_plot=True, p=n_clusters, truncate_mode='lastp')

    temp = {R["leaves"][ii]: (centroids[ii] + 1, R["ivl"][ii] if '(' in R["ivl"][ii] else '(1)') for ii in range(len(R["leaves"]))}

    def llf(xx):
        return "{} - {}".format(*temp[xx])

    cluster.hierarchy.dendrogram(Z, p=n_clusters, truncate_mode='lastp', leaf_label_func=llf)

    plt.ylim(bottom=y_val - 0.5)
    plt.ylabel('RMSD (Å)')
    plt.suptitle("RMSD clustering", fontsize=14, fontweight='bold')
    plt.legend(loc='upper right', framealpha=1, prop={'size': 9})
    plt.tight_layout()
    plt.show(block=False)


def get_rmsd(args):
    sup = SVDSuperimposer()
    i, j = args
    sup.set(structure_coords[i], structure_coords[j])
    global counter
    with counter.get_lock():
        counter.value += 1
    sup.run()
    return i, j, sup.get_rms()


def load_structure_coords(filename):
    coords = dict()
    current_model = 0

    with open(filename, 'r') as file:
        file.readline()
        file.readline()
        for line in file:
            line = line.strip().split(' ')
            if len(line) == 4:
                x, y, z = line[1:]
                coords.setdefault(current_model, [])
                coords[current_model].append(np.asarray([float(x), float(y), float(z)], dtype=np.float32))
            else:
                current_model += 1
                file.readline()  # Skip second line of header

    for model in coords.keys():
        coords[model] = np.asarray(coords[model], dtype=np.float32)

    return coords


def init(args):
    global counter
    counter = args


def compute_rmsd_dist_matrix(logger, pdb_id):
    mtrx_file = "/tmp/ring/{}.npy".format(pdb_id)

    logger.log("Loading structure")

    filename = "/tmp/ring/{}.xyz".format(pdb_id)

    global structure_coords

    structure_coords = load_structure_coords(filename)

    n_models = len(structure_coords.keys())
    X = np.zeros((n_models, n_models))
    indexes = np.tril_indices(n_models)

    args = [(i, j,) for i, j in zip(indexes[0], indexes[1])]

    logger.log("Computing distance matrix")

    counter = mp.Value('i', 0)
    with mp.Pool(mp.cpu_count() - 1, initializer=init, initargs=(counter,)) as p:
        results = p.map_async(get_rmsd, iterable=args, chunksize=100)
        while not results.ready():
            val = counter.value / len(indexes[0]) * 100
            logger.progress(val)
            time.sleep(1)

    logger.close_progress()
    del structure_coords

    for i, j, v in results.get():
        X[i, j] = v
    X += X.transpose()
    np.fill_diagonal(X, 0)
    np.save(mtrx_file, X)
    logger.log('Done')


def load_rmsd_dis_matrix(logger, pdb_id):
    mtrx_file = "/tmp/ring/{}.npy".format(pdb_id)
    logger.log("Loading distance matrix")
    X = np.load(mtrx_file)
    return X


def cluster_states_obj(logger, pdb_id, method, rmsd_val=None, desired_clusters=None):
    logger.disable_window()
    X = load_rmsd_dis_matrix(logger, pdb_id)

    logger.log("Operation started, please wait")

    result_labels, _, cut_heights, repr_labels = hierarchy_optimization(X, logger, height=rmsd_val,
                                                                        desired_clusters=desired_clusters,
                                                                        method=method)

    if desired_clusters is not None and desired_clusters not in result_labels:
        logger.log("The number of cluster has to be in the range 2 - {} (inclusive)".format(max(result_labels)))
        return

    obj_name = "{}_cl".format(pdb_id.strip('_ca'))
    cmd.delete(obj_name)

    repr_labels = np.asarray(repr_labels) + 1

    for state in repr_labels:
        cmd.create(obj_name, pdb_id.strip('_ca'), source_state=state, target_state=-1, copy_properties=True)

    logger.log("Created new object with representative states from original object")
    logger.log("States used from original object: {}".format(sorted(repr_labels)))
    logger.enable_window()


if __name__ == '__main__':
    class Logger:
        def __init__(self):
            pass

        @staticmethod
        def log(s, warning=False, error=False):
            print(s)

        @staticmethod
        def progress(n):
            print(n)

        def close_progress(self):
            pass

        def disable_window(self):
            pass

        def enable_window(self):
            pass


    temporary = Logger()
    cluster_distribution_heatmap(temporary, "trj_ca", 'complete', desired_clusters=25)
