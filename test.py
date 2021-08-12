intTypeMap = {
        "IONIC"    : "blue",
        "SSBOND"   : "yellow",
        "PIPISTACK": "orange",
        "PICATION" : "red",
        "HBOND"    : "cyan",
        "VDW"      : "gray50",
        "IAC"      : "white"
}

if __name__ == '__main__':
    obj = "2h9r"
    frames = 10
    min_presence = 0.05
    max_presence = 0.95
    coeff_thresh = 0.5
    p_thresh = 0.3,
    import seaborn as sn
    from matplotlib import pyplot as plt

    try:
        import pandas as pd
        import numpy as np
        from scipy.stats import pearsonr
    except ImportError:
        print("To run this you have to install pandas, numpy and scipy in python")
        exit(0)

    all_cm = dict()
    try:
        for int_type in intTypeMap.keys():
            all_cm[int_type] = pd.read_csv('/tmp/ring/md/{}.cm_{}'.format(obj, int_type), sep=' ', header=None)
    except FileNotFoundError:
        exit(0)

    contacts_sparse = dict()
    for j in range(1, frames + 1):
        for int_type in intTypeMap.keys():
            df = all_cm[int_type][all_cm[int_type][0] == j]
            names = df[1]
            names = [x.replace(':_:', ':') for x in names]
            df = df.iloc[:, 2:]
            matrix = df.values
            matrix[np.triu_indices(matrix.shape[0])] = 0
            for i in np.argwhere(matrix > 0):
                int_id1 = names[i[0]].split(':')
                int_id1 = (int_id1[0], int(int_id1[1]), int_id1[2])
                int_id2 = names[i[1]].split(':')
                int_id2 = (int_id2[0], int(int_id2[1]), int_id2[2])
                tmp = tuple(sorted([int_id1, int_id2]))
                contacts_sparse.setdefault(tmp, dict())
                contacts_sparse[tmp].setdefault(j - 1, 0)
                contacts_sparse[tmp][j - 1] += 1

    to_pop = []
    for k, v in contacts_sparse.items():
        presence = len(v) / frames
        if max_presence < presence or presence < min_presence:
            to_pop.append(k)

    for k in to_pop:
        contacts_sparse.pop(k)

    z = np.zeros((len(contacts_sparse), frames))
    for i, contacts_for_frame in enumerate(contacts_sparse.values()):
        for j, contacts in contacts_for_frame.items():
            z[i, j] = contacts

    coeffs_matr = np.ones((z.shape[0], z.shape[0])) * np.nan
    p_matr = np.ones((z.shape[0], z.shape[0])) * np.nan

    for i in range(z.shape[0]):
        for j in range(z.shape[0]):
            if i != j:
                corr_coeff, p_val = pearsonr(z[i], z[j])
                if p_val < p_thresh and (corr_coeff > coeff_thresh or corr_coeff < -coeff_thresh):
                    coeffs_matr[i, j] = corr_coeff
                    p_matr[i, j] = p_val

    ticks = ["{}/{}/{} - {}/{}/{}".format(x[0], x[1], x[2], y[0], y[1], y[2]) for (x, y) in contacts_sparse.keys()]
    ticks = np.array(ticks)
    hm = sn.heatmap(coeffs_matr, xticklabels=ticks, yticklabels=ticks)
    hm.set_xticklabels(hm.get_xmajorticklabels(), fontsize=8)
    hm.set_yticklabels(hm.get_ymajorticklabels(), fontsize=8)
