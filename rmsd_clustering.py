import multiprocessing as mp
import os
import time

import matplotlib.cm as cm
import numpy as np
import seaborn as sn
from Bio.SVDSuperimposer import SVDSuperimposer
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import cluster
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score

n_best = 20
structure_coords = dict()
counter: mp.Value
window = None


def cm_to_inch(x):
    return x / 2.54


def hierarchy_optimization(X):
    result_label = dict()
    cut_heights = dict()

    range_n_clusters = range(2, len(X))
    Z = cluster.hierarchy.linkage(squareform(X), optimal_ordering=True, method='complete')

    clusters = cluster.hierarchy.cut_tree(Z, n_clusters=range_n_clusters)
    for i, n_clusters in enumerate(range_n_clusters):
        cluster_labels = clusters[:, i].flatten()

        silhouette_avg = silhouette_score(X, cluster_labels, metric='precomputed')
        result_label.setdefault(n_clusters, (silhouette_avg, cluster_labels))

        tmp = cluster.hierarchy.dendrogram(Z, p=n_clusters, truncate_mode='lastp', no_plot=True)
        last_h = min(list(filter(lambda x: x != 0, [item for sublist in tmp['dcoord'] for item in sublist])))
        cut_heights.setdefault(n_clusters, last_h)

    return result_label, Z, cut_heights


def get_cluster_labels_for_rmsd(logger, pdb_id, rmsd_val):
    mtrx_file = "/tmp/ring/{}.npy".format(pdb_id)
    if not os.path.exists(mtrx_file):
        logger.log('Run the clustering calculation first', error=True)
        return

    X = get_rmsd_dist_matrix(logger, pdb_id)
    result_labels, _, n_cluster_height = hierarchy_optimization(X)
    point = 0
    for n, h in n_cluster_height.items():
        if h <= rmsd_val:
            point = n
            logger.log('Number of clusters for selected RMSD cut: {}'.format(n), warning=True)
    return result_labels[point][2]


def get_cluster_labels_for_n_cluster(logger, pdb_id, n_cluster):
    mtrx_file = "/tmp/ring/{}.npy".format(pdb_id)
    if not os.path.exists(mtrx_file):
        logger.log('Run the clustering calculation first', error=True)
        return

    X = get_rmsd_dist_matrix(logger, pdb_id)
    result_labels, *_ = hierarchy_optimization(X)

    if n_cluster not in result_labels.keys():
        logger.log('Number of clusters is not in the range 2 - {}'.format(max(result_labels.keys())), error=True)
        raise ValueError

    return result_labels[n_cluster][2]


def cluster_distribution_heatmap(logger, pdb_id, x_len=50):
    mtrx_file = "/tmp/ring/{}.npy".format(pdb_id)
    if not os.path.exists(mtrx_file):
        logger.log('Run the clustering calculation first', error=True)
        return

    X = get_rmsd_dist_matrix(logger, pdb_id)

    labels, Z, _ = hierarchy_optimization(X)
    labels = sorted([(v1, v2) for (v1, v2) in labels.values()], key=lambda x: x[0], reverse=True)
    labels = labels[2][1]

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

    colors = []
    n_clusters = len(set(labels.tolist()))
    for i in range(n_clusters):
        colors.append(cm.nipy_spectral((float(i) + 1) / (n_clusters + 1)))
    cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

    sn.set(font_scale=0.8)

    ax1 = sn.heatmap(best, cmap=cmap, linewidths=.5, linecolor='white', ax=ax1, square=True, cbar=False, annot=False)
    ax1.set_yticks(np.arange(0.5, best.shape[0] + 0.5))
    ax1.set_yticklabels(np.arange(1, best.size + 1, x_len), rotation=0, fontsize=10)
    ax1.set_xticks(np.arange(0.5, best.shape[1] + 0.5))
    ax1.set_xticklabels(range(1, best.shape[1] + 1))

    remove_spines(ax1)

    plt.xlabel('States')
    plt.title("Structure {} - {} clusters".format(pdb_id, len(set(labels))),
              fontsize=14)
    plt.tight_layout()
    plt.grid(False)
    plt.show()


def hierarchy_cut_plot(logger, pdb_id):
    mtrx_file = "/tmp/ring/{}.npy".format(pdb_id)
    if not os.path.exists(mtrx_file):
        logger.log('Run the clustering calculation first', error=True)
        return

    X = get_rmsd_dist_matrix(logger, pdb_id)

    result_labels, Z, cut_heights = hierarchy_optimization(X)
    result_labels = sorted([(n_cluster, silh_val) for n_cluster, (silh_val, _) in result_labels.items()],
                           key=lambda x: x[1], reverse=True)
    colors = []

    plt.close()
    plt.style.use('default')
    for i in range(2, n_best):
        color = cm.nipy_spectral(float(i + 1) / (n_best + 1))
        colors.append(color)
        plt.axhline(y=cut_heights[i], color=color, linestyle="--", zorder=0,
                    label="{} clusters, silh: {:.3f}".format(result_labels[i][0], result_labels[i][1]),
                    linewidth=1.3)
    cluster.hierarchy.dendrogram(Z, distance_sort=True, p=50, truncate_mode='lastp')
    plt.ylabel('RMSD (Å)')
    plt.suptitle("RMSD clustering", fontsize=14, fontweight='bold')
    plt.legend(loc='upper right', framealpha=1, prop={'size': 6})
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
    else:
        logger.log("Loading distance matrix")
        X = np.load(mtrx_file)
    return X


if __name__ == '__main__':
    class Logger:
        def __init__(self):
            pass

        def log(self, s, warning=False, error=False):
            print(s)


    logger = Logger()
    cluster_distribution_heatmap(logger, "trj_ca", 50)
