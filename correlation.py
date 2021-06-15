import pandas as pd
import numpy as np
import seaborn as sn
from matplotlib import pyplot as plt
from tqdm import tqdm
from scipy.stats import pearsonr


def calculate_correlation(obj, frames, min_presence=0.05, max_presence=0.95, coeff_thresh=0.5, p_thresh=0.05):
    contacts_sparse = dict()
    for j in tqdm(range(1, frames + 1)):
        try:
            df = pd.read_csv('/tmp/ring/{}.cif_cm_m{}_HBOND'.format(obj, j), sep=' ', header=None)
        except FileNotFoundError:
            print("Run ring on all the states first!")
            return
        names = df[0]
        names = [x.split('_')[0][:-1] for x in names]
        df = df.iloc[:, 1:]
        matrix = df.values
        matrix[np.triu_indices(matrix.shape[0])] = 0
        for i in np.argwhere(matrix > 0):
            contacts_sparse.setdefault((names[i[0]], names[i[1]]), [])
            contacts_sparse[(names[i[0]], names[i[1]])].append((j - 1, matrix[i[0], i[1]]))

    to_pop = []
    for k, v in contacts_sparse.items():
        presence = len(v) / frames
        if max_presence < presence or presence < min_presence:
            to_pop.append(k)

    for k in to_pop:
        contacts_sparse.pop(k)

    z = np.zeros((len(contacts_sparse), frames))
    for i, v in enumerate(contacts_sparse.values()):
        for j, contact in v:
            z[i, j] = contact

    coeffs_matr = np.empty((z.shape[0], z.shape[0]))
    coeffs_matr.fill(np.nan)
    for i, j in zip(np.tril_indices(z.shape[0])[0], np.tril_indices(z.shape[0])[1]):
        corr_coeff, p_val = pearsonr(z[i], z[j])
        if p_val < p_thresh and (corr_coeff > coeff_thresh or corr_coeff < -coeff_thresh):
            coeffs_matr[i, j] = corr_coeff

    ticks = ["{}_{}".format(x, y) for (x, y) in contacts_sparse.keys()]
    hm = sn.heatmap(coeffs_matr, xticklabels=ticks, yticklabels=ticks)
    hm.set_xticklabels(hm.get_xmajorticklabels(), fontsize=8)
    hm.set_yticklabels(hm.get_ymajorticklabels(), fontsize=8)
    plt.show()
    return ticks, coeffs_matr


if __name__ == '__main__':
    obj = '2l3y'
    frames = 52

    a, b = calculate_correlation(obj, frames)

    indexes = np.argwhere(np.isnan(b))
