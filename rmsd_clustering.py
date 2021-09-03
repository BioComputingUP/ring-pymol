import multiprocessing as mp
import os
import time

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


def hierarchy_optimization(X, get_heights=False, get_center_label=False):
    result_label = dict()
    cut_heights = dict()

    range_n_clusters = list(range(2, len(X)))
    Z = cluster.hierarchy.linkage(squareform(X), optimal_ordering=True, method='complete')

    centroid_clusters = dict()
    names = range(len(X))
    clusters = cluster.hierarchy.cut_tree(Z, n_clusters=range_n_clusters)
    for i, n_clusters in enumerate(range_n_clusters):
        cluster_labels = clusters[:, i].flatten()

        silhouette_avg = silhouette_score(X, cluster_labels, metric='precomputed')
        result_label.setdefault(n_clusters, (silhouette_avg, cluster_labels))

        if get_heights:
            tmp = cluster.hierarchy.dendrogram(Z, p=n_clusters, truncate_mode='lastp', no_plot=True)
            last_h = min(list(filter(lambda x: x != 0, [item for sublist in tmp['dcoord'] for item in sublist])))
            cut_heights.setdefault(n_clusters, last_h)

        if get_center_label:
            for counter in range(n_clusters):
                nameList = list(zip(names, cluster_labels))
                mask = np.array([i == counter for i in cluster_labels])
                idx = np.argmin(sum(X[:, mask][mask, :]))
                sublist = [name for (name, label) in nameList if label == counter]
                centroid_clusters.setdefault(n_clusters, []).append(sublist[idx])

    return result_label, Z, cut_heights, centroid_clusters


def cluster_distribution_heatmap(logger, pdb_id, rmsd_val=None, desired_clusters=None, x_len=50):
    logger.disable_window()

    X = get_rmsd_dist_matrix(logger, pdb_id)

    labels, Z, cut_heights, _ = hierarchy_optimization(X, get_heights=True)

    if desired_clusters is not None:
        if desired_clusters not in labels:
            logger.log("The number of cluster has to be in the range 2 - {} (inclusive)".format(max(labels)))
            return
    else:
        point = max(labels)
        for n, h in cut_heights.items():
            if h <= rmsd_val:
                point = n
                break
        logger.log('Number of clusters for selected RMSD cut: {}'.format(point), warning=True)

    n_clusters = desired_clusters if desired_clusters is not None else point

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

    ax1 = sn.heatmap(best, cmap=cmap, linewidths=.5, linecolor='white', ax=ax1, square=True, cbar=False, annot=False)
    ax1.set_yticks(np.arange(0.5, best.shape[0] + 0.5))
    ax1.set_yticklabels(np.arange(1, best.size + 1, x_len), rotation=0, fontsize=10)
    ax1.set_xticks(np.arange(0.5, best.shape[1] + 0.5))
    ax1.set_xticklabels(range(1, best.shape[1] + 1))

    remove_spines(ax1)
    logger.enable_window()

    plt.xlabel('States')
    plt.title("Structure {} - {} clusters".format(pdb_id, len(set(labels))),
              fontsize=14)
    plt.tight_layout()
    plt.grid(False)
    plt.show()


def hierarchy_cut_plot(logger, pdb_id, rmsd_val=None, desired_clusters=None):
    logger.disable_window()
    X = get_rmsd_dist_matrix(logger, pdb_id)

    result_labels, Z, cut_heights, _ = hierarchy_optimization(X, get_heights=True)

    if desired_clusters is not None:
        if desired_clusters in result_labels:
            silh_val = result_labels[desired_clusters][0]
            y_val = cut_heights[desired_clusters]
        else:
            logger.log("The number of cluster has to be in the range 2 - {} (inclusive)".format(max(result_labels)))
            return
    else:
        point = max(result_labels)
        for n, h in cut_heights.items():
            if h <= rmsd_val:
                point = n
                break
        logger.log('Number of clusters for selected RMSD cut: {}'.format(point), warning=True)
        silh_val = result_labels[point][0]
        y_val = cut_heights[point]

    logger.enable_window()

    n_clusters = desired_clusters if desired_clusters is not None else point
    plt.close()
    plt.style.use('default')
    plt.axhline(y=y_val, linestyle="--", zorder=0,
                label="{} clusters, silh: {:.3f}".format(n_clusters,
                                                         silh_val),
                linewidth=1.3)
    cluster.hierarchy.dendrogram(Z, distance_sort=True, p=n_clusters, truncate_mode='lastp')
    plt.ylabel('RMSD (Å)')
    plt.suptitle("RMSD clustering", fontsize=14, fontweight='bold')
    plt.legend(loc='upper right', framealpha=1, prop={'size': 9})
    plt.tight_layout()
    plt.show()


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


def get_rmsd_dist_matrix(logger, pdb_id):
    mtrx_file = "/tmp/ring/{}.npy".format(pdb_id)

    if not os.path.exists(mtrx_file):
        logger.disable_window()
        logger.log("Loading structure")

        filename = "/tmp/ring/{}.xyz".format(pdb_id)

        global structure_coords

        structure_coords = load_structure_coords(filename)

        n_models = len(structure_coords.keys())
        X = np.zeros((n_models, n_models))
        indexes = np.tril_indices(n_models)

        args = [(i, j,) for i, j in zip(indexes[0], indexes[1])]

        logger.log("Computing dist matrix")

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
        logger.enable_window()
    else:
        logger.log("Loading distance matrix")
        X = np.load(mtrx_file)
    return X


def cluster_states_obj(logger, pdb_id, rmsd_val=None, desired_clusters=None):
    logger.disable_window()
    X = get_rmsd_dist_matrix(logger, pdb_id)

    logger.log("Operation started, please wait")

    result_labels, _, cut_heights, repr_labels = hierarchy_optimization(X, get_heights=True, get_center_label=True)

    if desired_clusters is not None:
        if desired_clusters not in result_labels:
            logger.log("The number of cluster has to be in the range 2 - {} (inclusive)".format(max(result_labels)))
            return
    else:
        point = max(result_labels)
        for n, h in cut_heights.items():
            if h <= rmsd_val:
                point = n
                break
        logger.log('Number of clusters for selected RMSD cut: {}'.format(point), warning=True)

    n_clusters = desired_clusters if desired_clusters is not None else point

    obj_name = "{}_cl".format(pdb_id.strip('_ca'))
    cmd.delete(obj_name)

    for state in repr_labels[n_clusters]:
        cmd.create(obj_name, pdb_id.strip('_ca'), source_state=state + 1, target_state=-1, copy_properties=True)

    logger.log("Created new object with representative states from original object")
    logger.log("States used from original object: {}".format(sorted(repr_labels[n_clusters])))
    logger.enable_window()


if __name__ == '__main__':
    class Logger:
        def __init__(self):
            pass

        def log(self, s, warning=False, error=False):
            print(s)


    temporary = Logger()
    hierarchy_cut_plot(temporary, "2h9r_ca", rmsd_val=3.5)
