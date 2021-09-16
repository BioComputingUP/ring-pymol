import json
import math
from json import JSONEncoder
from typing import Dict, List, Union

import matplotlib.cm as cm
import networkx as nx
import numpy as np
import pandas as pd
from PyQt5.QtGui import QColor
from matplotlib.colors import ListedColormap
from pymol.cgo import *
from scipy.stats import pearsonr

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


def _default(self, obj):
    return getattr(obj.__class__, "to_json", _default.default)(obj)


_default.default = JSONEncoder().default
JSONEncoder.default = _default


class Node:
    def __init__(self, *args):
        if len(args) == 1:
            self.init_string(*args)
        else:
            self.init_args(*args)
        self.atom = None

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

    def to_json(self):
        return self.__repr__()


class Edge:
    def __init__(self, *args):
        self.node1 = None
        self.node2 = None
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

    obj = [BEGIN, LINES, COLOR] + tup_color
    for interaction in interactions:
        valid = True
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
    obj.append(END)
    cmd.load_cgo(obj, object_name, state=state, zoom=False)


def calculate_correlation(obj, frames, min_presence=0.05, max_presence=0.95, coeff_thresh=0.5, p_thresh=0.3,
                          int_type="HBOND"):
    all_cm = dict()
    nodes = []
    if int_type == "ALL":
        to_read = intTypeMap.keys()
    else:
        to_read = [int_type]

    for interaction in to_read:
        all_cm[interaction] = pd.read_csv('/tmp/ring/md/{}.cm_{}'.format(obj, interaction), sep=' ',
                                          header=None)
        if len(nodes) == 0:
            nodes = all_cm[interaction][all_cm[interaction][0] == 1][1]
            nodes = [Node(x) for x in nodes]

    conn_freq = get_freq(obj)
    contacts_sparse = dict()
    for frame in range(0, frames):
        for interaction in to_read:
            df = all_cm[interaction][all_cm[interaction][0] == frame + 1]
            df = df.iloc[:, 2:]
            matrix = df.values
            matrix[np.triu_indices(matrix.shape[0])] = 0
            for i, j in np.argwhere(matrix > 0):
                node1 = nodes[i]
                node2 = nodes[j]
                edge = Edge(sorted([node1, node2]))
                if min_presence < conn_freq[interaction][edge] < max_presence:
                    contacts_sparse.setdefault(edge, dict())
                    contacts_sparse[edge].setdefault(frame, 0)
                    contacts_sparse[edge][frame] += 1

    z = np.zeros((len(contacts_sparse), frames))
    for i, contacts_for_frame in enumerate(contacts_sparse.values()):
        for j, contacts in contacts_for_frame.items():
            z[i, j] = contacts

    coeffs_matr = np.ones((z.shape[0], z.shape[0])) * np.nan
    p_matr = np.ones((z.shape[0], z.shape[0])) * np.nan

    indexes = np.triu_indices(z.shape[0], k=1)
    for i, j in zip(indexes[0], indexes[1]):

        corr_coeff, p_val = pearsonr(z[i], z[j])
        if p_val < p_thresh and (corr_coeff > coeff_thresh or corr_coeff < -coeff_thresh):
            coeffs_matr[i, j] = corr_coeff
            p_matr[i, j] = p_val
    p_matr[np.tril_indices_from(p_matr)] = p_matr[np.triu_indices_from(p_matr)]
    coeffs_matr[np.tril_indices_from(coeffs_matr)] = coeffs_matr[np.triu_indices_from(coeffs_matr)]
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


def generate_colormap(number_of_distinct_colors: int = 80):
    if number_of_distinct_colors == 0:
        number_of_distinct_colors = 80

    number_of_distinct_colors_min = max(8, number_of_distinct_colors)

    number_of_shades = 7
    number_of_distinct_colors_with_multiply_of_shades = int(
            math.ceil(number_of_distinct_colors_min / number_of_shades) * number_of_shades)

    # Create an array with uniformly drawn floats taken from <0, 1) partition
    linearly_distributed_nums = np.arange(
            number_of_distinct_colors_with_multiply_of_shades) / number_of_distinct_colors_with_multiply_of_shades

    # We are going to reorganise monotonically growing numbers in such way that there will be single array with saw-like pattern
    #     but each saw tooth is slightly higher than the one before
    # First divide linearly_distributed_nums into number_of_shades sub-arrays containing linearly distributed numbers
    arr_by_shade_rows = linearly_distributed_nums.reshape(number_of_shades,
                                                          number_of_distinct_colors_with_multiply_of_shades // number_of_shades)

    # Transpose the above matrix (columns become rows) - as a result each row contains saw tooth with values slightly higher than row above
    arr_by_shade_columns = arr_by_shade_rows.T

    # Keep number of saw teeth for later
    number_of_partitions = arr_by_shade_columns.shape[0]

    # Flatten the above matrix - join each row into single array
    nums_distributed_like_rising_saw = arr_by_shade_columns.reshape(-1)

    # HSV colour map is cyclic (https://matplotlib.org/tutorials/colors/colormaps.html#cyclic), we'll use this property
    initial_cm = cm.hsv(nums_distributed_like_rising_saw)

    lower_partitions_half = number_of_partitions // 2
    upper_partitions_half = number_of_partitions - lower_partitions_half

    # Modify lower half in such way that colours towards beginning of partition are darker
    # First colours are affected more, colours closer to the middle are affected less
    lower_half = lower_partitions_half * number_of_shades
    for i in range(3):
        initial_cm[0:lower_half, i] *= np.arange(0.2, 1, 0.8 / lower_half)

    # Modify second half in such way that colours towards end of partition are less intense and brighter
    # Colours closer to the middle are affected less, colours closer to the end are affected more
    for i in range(3):
        for j in range(upper_partitions_half):
            modifier = np.ones(number_of_shades) - initial_cm[lower_half + j * number_of_shades: lower_half + (
                    j + 1) * number_of_shades, i]
            modifier = j * modifier / upper_partitions_half
            initial_cm[lower_half + j * number_of_shades: lower_half + (j + 1) * number_of_shades, i] += modifier

    return ListedColormap(initial_cm, N=number_of_distinct_colors)


def remap(value, low1, high1, low2, high2):
    return low2 + (value - low1) * (high2 - low2) / (high1 - low1)


def export_network_graph(model):
    G = nx.MultiGraph()

    # Add the nodes to the graph
    file_pth = "/tmp/ring/" + model + ".cif_ringNodes"
    df = pd.read_csv(file_pth, sep='\t')
    df = df.groupby('NodeId').mean()

    for (nodeId, _, degree, *_) in df.itertuples(index=True):
        node = Node(nodeId)
        G.add_node(node, degree=round(degree, 3), chain=node.chain, resi=node.resi, resn=node.resn)

    # Add the edges to the graph
    file_pth = "/tmp/ring/" + model + ".cif_ringEdges"
    df = pd.read_csv(file_pth, sep='\t')

    distance_dict = dict()
    mean_distance = df.groupby(['NodeId1', 'NodeId2', 'Interaction']).mean()
    for (nodeId, distance, *_) in mean_distance.itertuples(index=True, name='Distance'):
        nodeId1, nodeId2, interaction = nodeId
        intType = interaction.split(":")[0]
        node1 = Node(nodeId1)
        node2 = Node(nodeId2)
        edge = Edge(node1, node2)
        distance_dict.setdefault(intType, dict()).setdefault(edge, distance)

    conn_freq = get_freq(model)

    sawn = set()
    df = df.groupby(["NodeId1", "Interaction", "NodeId2"]).sum()
    for (ids, *_) in df.itertuples(index=True):
        nodeId1, interaction, nodeId2 = ids
        intType = interaction.split(":")[0]
        node1 = Node(nodeId1)
        node2 = Node(nodeId2)
        edge = Edge(node1, node2)
        key = (edge, intType)
        if key not in sawn:
            G.add_edge(node1, node2, interaction=intType, frequency=round(conn_freq[intType][edge], 3),
                       distance=round(distance_dict[intType][edge], 3))
            sawn.add(key)

    with open("/tmp/ring/{}.json".format(model), 'w+') as f:
        json.dump(nx.cytoscape_data(G), f)


if __name__ == '__main__':
    # export_network_graph('2h9r')
    calculate_correlation("trj_cl", 20)
