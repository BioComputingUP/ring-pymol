from PyQt5.QtGui import QColor
from pymol.cgo import *

intTypeMap = {
        "IONIC"    : "blue",
        "SSBOND"   : "yellow",
        "PIPISTACK": "orange",
        "PICATION" : "red",
        "HBOND"    : "cyan",
        "VDW"      : "gray50",
        "IAC"      : "white"
}

intColorMap = {
        "IONIC"    : "blue",
        "SSBOND"   : "yellow",
        "PIPISTACK": "orange",
        "PICATION" : "red",
        "HBOND"    : "cyan",
        "VDW"      : "gray",
        "IAC"      : "black"
}


def get_freq(obj):
    conn_freq = dict()
    for inter in intTypeMap.keys():
        try:
            conn_freq.setdefault(inter, dict())
            with open("/tmp/ring/md/{}.gfreq_{}".format(obj, inter), 'r') as f:
                for line in f:
                    edge1, _, edge2, perc = line.split('\t')
                    edge1 = edge1.rsplit(':', 2)[0]
                    edge2 = edge2.rsplit(':', 2)[0]
                    conn_freq[inter].setdefault((edge1, edge2), float(perc))
        except FileNotFoundError:
            raise FileNotFoundError
    return conn_freq


def get_freq_combined(obj, bond, interchain=False, intrachain=False):
    import math
    conn_freq = dict()
    try:
        with open("/tmp/ring/md/{}.gfreq_{}".format(obj, bond), 'r') as f:
            for line in f:
                edge1, _, edge2, perc = line.split('\t')
                edge1 = edge1.replace(":_:", ":")
                edge2 = edge2.replace(":_:", ":")
                chain1 = edge1.split(':')[0]
                chain2 = edge2.split(':')[0]
                if interchain and chain1 == chain2:
                    continue
                if intrachain and chain1 != chain2:
                    continue
                conn_freq.setdefault(edge1, [])
                conn_freq[edge1].append(float(perc))
    except FileNotFoundError:
        raise FileNotFoundError
    for k, v in conn_freq.items():
        conn_freq[k] = 1 - math.prod([(1 - x) for x in v])

    return conn_freq


def draw_links(interactions, color, object_name, coords, state):
    from pymol import cmd

    tup_color = []
    if type(color) is str:
        try:
            tup_color = list(map(float, color.replace('(', '').replace(')', '').split(',')))
        except ValueError:
            tup_color = list(cmd.get_color_tuple(color))
    elif type(color) is list or type(color) is tuple:
        tup_color = list(color)

    not_present = 0
    obj = [BEGIN, LINES, COLOR] + tup_color
    for interaction in interactions:
        valid = True
        # obj.extend([BEGIN, LINES, COLOR] + tup_color)
        if "," in interaction[0]:
            coord1 = ([float(x) for x in interaction[0].split(',')],)
        else:
            try:
                coord1 = (coords[interaction[0]],)
            except KeyError:
                valid = False

        if "," in interaction[1]:
            coord2 = ([float(x) for x in interaction[1].split(',')],)
        else:
            try:
                coord2 = (coords[interaction[1]],)
            except KeyError:
                valid = False

        if valid:
            for x, y in zip(coord1, coord2):
                obj.extend([VERTEX] + x + [VERTEX] + y)
        else:
            not_present += 1
    obj.append(END)
    cmd.load_cgo(obj, object_name, state=state, zoom=False)
    return not_present


def calculate_correlation(obj, frames, min_presence=0.05, max_presence=0.95, coeff_thresh=0.5, p_thresh=0.3,
                          int_type="HBOND"):
    # import seaborn as sn
    # from matplotlib import pyplot as plt
    try:
        import pandas as pd
        import numpy as np
        from scipy.stats import pearsonr
    except ImportError:
        print("To run this you have to install pandas, numpy and scipy in python")
        return

    all_cm = dict()
    try:
        if int_type == "ALL":
            for interaction in intTypeMap.keys():
                all_cm[interaction] = pd.read_csv('/tmp/ring/md/{}.cm_{}'.format(obj, interaction), sep=' ',
                                                  header=None)
        else:
            all_cm[int_type] = pd.read_csv('/tmp/ring/md/{}.cm_{}'.format(obj, int_type), sep=' ', header=None)
    except FileNotFoundError:
        return

    contacts_sparse = dict()
    for j in range(1, frames + 1):
        if int_type == "ALL":
            interactions = intTypeMap.keys()
        else:
            interactions = [int_type]
        for interaction in interactions:
            df = all_cm[interaction][all_cm[interaction][0] == j]
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
    # hm = sn.heatmap(coeffs_matr, xticklabels=ticks, yticklabels=ticks)
    # hm.set_xticklabels(hm.get_xmajorticklabels(), fontsize=8)
    # hm.set_yticklabels(hm.get_ymajorticklabels(), fontsize=8)
    # plt.show()
    return ticks, coeffs_matr, p_matr


def get_bg_fg_colors(color):
    if color == 1:
        bk_color = QColor(232, 231, 252)
        fg_color = QColor(0, 0, 0)
        color = 2
    else:
        bk_color = QColor(255, 255, 255)
        fg_color = QColor(0, 0, 0)
        color = 1
    return bk_color, fg_color, color


def is_selection(string):
    return string[0] == "(" and string[-1] == ")"
