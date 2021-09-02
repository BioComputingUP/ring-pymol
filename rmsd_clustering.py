import multiprocessing as mp
import os

import matplotlib.cm as cm
import numpy as np
import seaborn as sn
from Bio.SVDSuperimposer import SVDSuperimposer
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap
from scipy import cluster
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score
from tqdm import tqdm

n_best = 3
structure_coords = dict()


def hierarchy_optimization(X, range_n_clusters):
    result_label = dict()

    Z = cluster.hierarchy.linkage(squareform(X), optimal_ordering=True, method='complete')

    clusters = cluster.hierarchy.cut_tree(Z, n_clusters=range_n_clusters)
    for i, n_clusters in enumerate(range_n_clusters):
        cluster_labels = clusters[:, i].flatten()

        silhouette_avg = silhouette_score(X, cluster_labels, metric='precomputed')
        # print("Hierarchical avg silhouette for {} clusters: {}".format(n_clusters, silhouette_avg))
        result_label.setdefault(n_clusters, (silhouette_avg, cluster_labels))

    return result_label, Z


def cluster_distribution_heatmap(logger, pdb_id, x_len=30):
    X = get_rmsd_dist_matrix(logger, pdb_id)

    range_n_clusters = range(2, len(X))
    labels, Z = hierarchy_optimization(X, range_n_clusters)
    labels = sorted([(v1, v2) for (v1, v2) in labels.values()], key=lambda x: x[0], reverse=True)
    labels = labels[2][1]

    def remove_spines(ax):
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['bottom'].set_color('none')

    plt.close()
    plt.style.use('default')
    ax1: Axes = plt.subplot()
    ax2: Axes = ax1.twinx()

    pad_size = x_len - len(labels) % x_len if len(labels) % x_len != 0 else 0
    labels = labels.astype('float32')
    best = np.pad(labels, (0, pad_size), mode='constant', constant_values=np.nan)
    best = np.reshape(best, (int(len(best) / x_len), x_len))

    colors = []
    n_clusters = len(set(labels.tolist()))
    for i in range(n_clusters):
        colors.append(cm.nipy_spectral((float(i) + 1) / (n_clusters + 1)))
    cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

    sn.set(font_scale=1)

    ax1 = sn.heatmap(best, cmap=cmap, linewidths=.5, linecolor='white', ax=ax1, square=True, cbar=False, annot=True,
                     fmt='.4g')
    ax1.set_yticks(np.arange(0.5, best.size / x_len + 0.5))
    ax1.set_yticklabels(np.arange(1, best.size + 1, x_len), rotation=0, fontsize=14)

    remove_spines(ax1)
    remove_spines(ax2)

    ax2.yaxis.tick_right()
    ax2.set_ylim(ax1.get_ylim())
    ax2.set_xticks([])
    ax2.set_yticks(np.arange(0.5, best.size / x_len + 0.5))
    ax2.set_yticklabels(np.arange(x_len, best.size + 1, x_len), rotation=0, fontsize=14)
    plt.tight_layout()
    plt.suptitle("Structure {} - {} clusters".format(pdb_id, len(set(labels))),
                 fontsize=18, fontweight='bold')
    plt.grid(False)
    plt.show()


def hierarchy_cut_plot(logger, pdb_id):
    X = get_rmsd_dist_matrix(logger, pdb_id)

    range_n_clusters = range(2, len(X))
    result_labels, Z = hierarchy_optimization(X, range_n_clusters)
    result_labels = sorted([(v1, v2) for (v1, v2) in result_labels.values()], key=lambda x: x[0], reverse=True)

    cut_heights = []
    for i in range(n_best):
        n_clusters = len(set(result_labels[i][1]))
        tmp = cluster.hierarchy.dendrogram(Z, p=n_clusters, truncate_mode='lastp', no_plot=True)
        cut_heights.append(tmp['dcoord'][0][1])
    colors = []

    plt.close()
    plt.style.use('default')
    for i, height in enumerate(cut_heights):
        color = cm.nipy_spectral(float(i + 1) / (len(cut_heights) + 1))
        colors.append(color)
        plt.axhline(y=height, color=color, linestyle="--", zorder=0,
                    label="Cl {} silh: {:.3f}".format(len(set(result_labels[i][1])), result_labels[i][0]),
                    linewidth=2.)
    cluster.hierarchy.dendrogram(Z, distance_sort=True)
    plt.suptitle("RMSD clustering", fontsize=14, fontweight='bold')
    plt.legend(loc='upper right', framealpha=1)
    plt.tight_layout()
    plt.show()


def f(args):
    sup = SVDSuperimposer()
    i, j = args
    sup.set(structure_coords[i], structure_coords[j])
    sup.run()
    return i, j, sup.get_rms()


def load_structure_coords(filename):
    coords = dict()
    current_model = 0

    with open(filename, 'r') as file:
        file.readline()
        file.readline()
        for line in tqdm(file):
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


def get_rmsd_dist_matrix(logger, pdb_id):
    mtrx_file = "/tmp/ring/{}.npy".format(pdb_id)

    if not os.path.exists(mtrx_file):
        logger.log("Loading structure")

        filename = "/tmp/ring/{}.xyz".format(pdb_id)

        logger.log("Getting coords")
        global structure_coords

        structure_coords = load_structure_coords(filename)

        n_models = len(structure_coords.keys())
        X = np.zeros((n_models, n_models))
        indexes = np.tril_indices(n_models)

        args = [(i, j,) for i, j in zip(indexes[0], indexes[1])]

        logger.log("Computing dist matrix")

        with mp.Pool(mp.cpu_count() - 1) as p:
            results = p.map(f, iterable=args, chunksize=100)

        del structure_coords

        for i, j, v in results:
            X[i, j] = v
        X += X.transpose()
        np.fill_diagonal(X, 0)
        np.save(mtrx_file, X)
    else:
        logger.log("Loading dist matrix")
        X = np.load(mtrx_file)
    return X
