import math
from typing import Dict, List, Union

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


class Node:
    def __init__(self, *args):
        if len(args) == 1:
            self.init_string(*args)
        else:
            self.init_args(*args)

    def init_string(self, string_id: str):
        if ':' in string_id:
            ids = string_id.strip().split(':')
        else:
            ids = string_id.strip().split('/')

        self.chain: str = ids[0]
        self.resi: int = int(ids[1])
        self.ins = None
        self.resn = None

        if len(ids) > 2:
            if len(ids[2]) == 3:
                self.resn: str = ids[2]
            else:
                self.ins: str = ids[2]
            if len(ids) > 3:
                self.resn: str = ids[3]

        if self.ins == '_':
            self.ins = None

    def init_args(self, chain: str, resi: Union[int, str], resn: str = None, ins: str = None):
        self.chain: str = chain
        self.resi: int = int(resi)
        self.ins: str = ins
        self.resn: str = resn

    def __lt__(self, other):
        if self.ins:
            return self.chain < other.chain or \
                   self.chain == other.chain and self.resi < other.resi or \
                   self.chain == other.chain and self.resi == other.resi and self.ins < other.ins
        return self.chain < other.chain or self.chain == other.chain and self.resi < other.resi

    def __le__(self, other):
        return self == other or self < other

    def __gt__(self, other):
        return other < self

    def __ge__(self, other):
        return self == other or other < self

    def __eq__(self, other):
        if not other:
            return False
        base = self.chain == other.chain and self.resi == other.resi and self.ins == other.ins
        if self.resn and other.resn:
            base = base and self.resn == other.resn
        return base

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        if self.resn:
            if self.ins:
                return "{}/{}/{}/{}".format(self.chain, self.resi, self.ins, self.resn)
            else:
                return "{}/{}/{}".format(self.chain, self.resi, self.resn)
        elif self.ins:
            return "{}/{}/{}".format(self.chain, self.resi, self.ins)
        return "{}/{}".format(self.chain, self.resi)

    def __hash__(self):
        if self.ins is not None:
            return hash((self.chain, self.resi, self.ins))
        return hash((self.chain, self.resi))

    def id_repr(self):
        return "{}/{}".format(self.chain, self.resi)

    def id_tuple(self):
        return self.chain, self.resi


class Edge:
    def __init__(self, *args):
        if len(args) == 2:
            self.init_nodes(*args)
        else:
            self.init_list(*args)

    def init_nodes(self, node1: Node, node2: Node):
        self.node1: Node = node1
        self.node2: Node = node2

    def init_list(self, sorted_node_list: List[Node]):
        if len(sorted_node_list) != 2:
            raise ValueError("Cannot create an Edge with more than two nodes")
        self.node1: Node = sorted_node_list[0]
        self.node2: Node = sorted_node_list[1]

    def __lt__(self, other):
        return self.node1 < other.node1 or (self.node1 == other.node1 and self.node2 < other.node2)

    def __le__(self, other):
        return self == other or self < other

    def __gt__(self, other):
        return other < self

    def __ge__(self, other):
        return self == other or other < self

    def __eq__(self, other):
        if not other:
            return False
        return (self.node1 == other.node1 and self.node2 == other.node2) or (
                self.node1 == other.node2 and self.node2 == other.node1)

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return "{} - {}".format(self.node1, self.node2)

    def __hash__(self):
        return hash((self.node1, self.node2))


def get_freq(obj, interchain=False, intrachain=False) -> Dict[str, Dict[Edge, float]]:
    conn_freq = dict()
    for inter in intTypeMap.keys():
        conn_freq.setdefault(inter, dict())
        with open("/tmp/ring/md/{}.gfreq_{}".format(obj, inter), 'r') as f:
            for line in f:
                node1, _, node2, perc = line.split('\t')
                node1 = Node(node1)
                node2 = Node(node2)
                edge = Edge(node1, node2)

                if intrachain and node1.chain != node2.chain:
                    continue
                if interchain and node1.chain == node2.chain:
                    continue

                conn_freq[inter].setdefault(edge, float(perc))
    return conn_freq


def get_freq_combined(obj, bond, interchain=False, intrachain=False, key_string=False):
    conn_freq = dict()
    try:
        with open("/tmp/ring/md/{}.gfreq_{}".format(obj, bond), 'r') as f:
            for line in f:
                node1, _, node2, perc = line.split('\t')
                node1 = Node(node1)
                node2 = Node(node2)
                if intrachain and node1.chain != node2.chain:
                    continue
                if interchain and node1.chain == node2.chain:
                    continue
                if not key_string:
                    conn_freq.setdefault(node1, [])
                    conn_freq[node1].append(float(perc))
                else:
                    conn_freq.setdefault(str(node1), [])
                    conn_freq[str(node1)].append(float(perc))

    except FileNotFoundError:
        raise FileNotFoundError
    for k, v in conn_freq.items():
        conn_freq[k] = 1 - math.prod([(1 - x) for x in v])

    return conn_freq


def get_freq_combined_all_interactions(obj):
    conn_freq = dict()
    for inter in intTypeMap.keys():
        with open("/tmp/ring/md/{}.gfreq_{}".format(obj, inter), 'r') as f:
            for line in f:
                node1, _, node2, perc = line.split('\t')
                node1 = Node(node1)
                node2 = Node(node2)
                edge = Edge(node1, node2)
                if node1.chain != node2.chain:
                    conn_freq.setdefault(edge, [])
                    conn_freq[edge].append(float(perc))

    all_freq = dict()
    for k, v in conn_freq.items():
        all_freq[k] = 1 - math.prod([(1 - x) for x in v])

    return all_freq


def get_node_names_ordered(obj):
    node_list = []
    with open("/tmp/ring/{}.cif_ringNodes".format(obj), 'r') as f:
        f.readline()
        for line in f:
            node_id, *_, model = line.strip().split("\t")
            if model == "1":
                node_list.append(Node(node_id))
            else:
                return node_list
    return node_list


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
            nodes = df[1]
            nodes = [Node(x) for x in nodes]
            df = df.iloc[:, 2:]
            matrix = df.values
            matrix[np.triu_indices(matrix.shape[0])] = 0
            for i in np.argwhere(matrix > 0):
                node1 = nodes[i[0]]
                node2 = nodes[i[1]]
                edge = Edge(sorted([node1, node2]))
                contacts_sparse.setdefault(edge, dict())
                contacts_sparse[edge].setdefault(j - 1, 0)
                contacts_sparse[edge][j - 1] += 1

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

    ticks = [str(edge) for edge in contacts_sparse.keys()]
    ticks = np.array(ticks)
    # hm = sn.heatmap(coeffs_matr, xticklabels=ticks, yticklabels=ticks)
    # hm.set_xticklabels(hm.get_xmajorticklabels(), fontsize=8)
    # hm.set_yticklabels(hm.get_ymajorticklabels(), fontsize=8)
    # plt.show()
    return list(contacts_sparse.keys()), coeffs_matr, p_matr


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
