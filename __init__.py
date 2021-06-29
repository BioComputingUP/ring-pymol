# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.
import pandas as pd
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


def __init_plugin__(app=None):
    """
    Add an entry to the PyMOL "Plugin" menu
    """
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Ring plugin', run_plugin_gui)


def draw_links(interactions, color, object_name, coords, state):
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
    import pandas as pd
    import numpy as np
    # import seaborn as sn
    # from matplotlib import pyplot as plt
    from scipy.stats import pearsonr

    contacts_sparse = dict()

    try:
        all_cm = pd.read_csv('/tmp/ring/md/{}.cm_{}'.format(obj, int_type), sep=' ', header=None)
    except FileNotFoundError:
        return

    for j in range(1, frames + 1):
        df = all_cm[all_cm[0] == j]
        names = df[1]
        names = [x.replace(':_:', ':') for x in names]
        df = df.iloc[:, 2:]
        matrix = df.values
        matrix[np.triu_indices(matrix.shape[0])] = 0
        for i in np.argwhere(matrix > 0):
            if names[i[0]].split(':')[0] < names[i[1]].split(':')[0]:
                tmp = (names[i[0]], names[i[1]])
            elif names[i[0]].split(':')[0] > names[i[1]].split(':')[0]:
                tmp = (names[i[1]], names[i[0]])
            else:
                if int(names[i[0]].split(':')[1]) <= int(names[i[1]].split(':')[1]):
                    tmp = (names[i[0]], names[i[1]])
                else:
                    tmp = (names[i[1]], names[i[0]])

            contacts_sparse.setdefault(tmp, [])
            contacts_sparse[tmp].append((j - 1, matrix[i[0], i[1]]))

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
    p_matr = np.empty((z.shape[0], z.shape[0]))
    p_matr.fill(np.nan)

    for i in range(z.shape[0]):
        for j in range(z.shape[0]):
            if i != j:
                corr_coeff, p_val = pearsonr(z[i], z[j])
                if p_val < p_thresh and (corr_coeff > coeff_thresh or corr_coeff < -coeff_thresh):
                    coeffs_matr[i, j] = corr_coeff
                    p_matr[i, j] = p_val

    ticks = ["{} - {}".format(x, y) for (x, y) in contacts_sparse.keys()]
    ticks = np.array(ticks)
    # hm = sn.heatmap(coeffs_matr, xticklabels=ticks, yticklabels=ticks)
    # hm.set_xticklabels(hm.get_xmajorticklabels(), fontsize=8)
    # hm.set_yticklabels(hm.get_ymajorticklabels(), fontsize=8)
    # plt.show()
    return ticks, coeffs_matr, p_matr


# global reference to avoid garbage collection of our dialog
dialog = None
dialog_corr = None
dialog_freq = None
correlations = dict()


def run_plugin_gui():
    """
    Open our custom dialog
    """
    global dialog
    global dialog_corr
    global dialog_freq

    if dialog is None and dialog_corr is None:
        dialog_corr = open_correlated_window()
        dialog_freq = open_frequency_window()
        dialog = make_dialog()
    dialog.adjustSize()
    dialog.show()


def open_correlated_window():
    from pymol.Qt import QtWidgets
    from pymol.Qt import QtCore
    from pymol.Qt.utils import loadUi

    dialog_corr = QtWidgets.QDialog()

    dialog_corr.setWindowFlags(dialog_corr.windowFlags() & QtCore.Qt.WindowMinimizeButtonHint)

    uifile = os.path.join(os.path.dirname(__file__), 'correlated.ui')
    loadUi(uifile, dialog_corr)

    return dialog_corr


def open_frequency_window():
    from pymol.Qt import QtWidgets
    from pymol.Qt import QtCore
    from pymol.Qt.utils import loadUi

    dialog_freq = QtWidgets.QDialog()
    dialog_freq.setWindowFlags(dialog_freq.windowFlags() & QtCore.Qt.WindowMinimizeButtonHint)

    uifile = os.path.join(os.path.dirname(__file__), 'frequency.ui')
    loadUi(uifile, dialog_freq)

    return dialog_freq


def make_dialog():
    from pymol import cmd
    from os import environ

    from pymol.Qt import QtWidgets
    from pymol.Qt import QtCore
    from PyQt5.QtGui import QColor
    from pymol.Qt.utils import loadUi

    import datetime

    app = QtWidgets.QApplication([])

    # create a new Window
    dialog = QtWidgets.QDialog()

    dialog.setWindowFlags(dialog.windowFlags() & QtCore.Qt.WindowMinimizeButtonHint)

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'plugin.ui')
    form = loadUi(uifile, dialog)

    def log(s: str, timed=True, error=False, warning=False, process=True):
        now = datetime.datetime.now()
        now = "{}:{}:{:02}".format(now.hour, now.minute, now.second)
        out_s = ""
        if timed:
            out_s += "{}  - ".format(now)
        if error:
            out_s += "ERROR: "
        if warning:
            out_s += "WARNING: "
        out_s += s

        form.console_log.insertItem(0, out_s)
        if error:
            form.console_log.item(0).setForeground(QColor(237, 67, 55))
        if warning:
            form.console_log.item(0).setForeground(QColor(255, 204, 0))

        if process:
            app.processEvents()

    def disable_window():
        dialog.main.setEnabled(False)
        app.processEvents()

    def enable_window():
        dialog.main.setEnabled(True)
        app.processEvents()

    def run():
        from pymol import stored

        import subprocess

        obj_name = form.pymol_obj.currentText()
        if not obj_name:
            log("Please select a Pymol object first!", error=True)
            return
        if not os.path.exists('/tmp/ring'):
            os.mkdir('/tmp/ring')

        file_pth = "/tmp/ring/" + obj_name + ".cif"

        form.main.setEnabled(False)
        form.button_start.setText("Running...")
        app.processEvents()

        log("Exporting pymol object {} in cif format ({})".format(obj_name, file_pth))
        cmd.save(filename=file_pth, selection=obj_name, state=0)
        log("Exporting done")

        stored.state = ''
        cmd.iterate_state(state=-1, selection=obj_name, expression='stored.state=state')

        nStates = cmd.count_states(obj_name)

        edge_policy = ""
        if form.best_edge.isChecked():
            edge_policy = "--best_edge"
        if form.multi_edge.isChecked():
            edge_policy = "--multi_edge"
        if form.all_edge.isChecked():
            edge_policy = "--all_edges"

        seq_sep = str(form.seq_separation.value())

        len_hbond = str(form.len_hbond.value())
        len_pica = str(form.len_pica.value())
        len_pipi = str(form.len_pipi.value())
        len_salt = str(form.len_salt.value())
        len_ss = str(form.len_ss.value())
        len_vdw = str(form.len_vdw.value())

        try:
            p = subprocess.Popen(
                    [form.ring_path.text(), "-i", file_pth, "--out_dir", "/tmp/ring/", "-g", seq_sep,
                     "-o", len_salt, "-s", len_ss, "-k", len_pipi, "-a", len_pica, "-b", len_hbond, "-w", len_vdw,
                     "--all_chains", edge_policy, "--all_models"], stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE, universal_newlines=True)
        except FileNotFoundError:
            log("Ring path is not correct!", error=True)
            form.main.setEnabled(True)
            return

        log("Ring generation started")

        while p.poll() is None:
            line = p.stderr.readline()
            if line != "":
                if "model" in line:
                    nModel = int(line.split("model ")[1].strip())
                    form.button_start.setText("Running on model {} | {:.2%}".format(nModel, nModel / nStates))
            app.processEvents()

        form.button_start.setText("Start Ring")
        log("Ring generation finished")
        visualize(first=True, log_iter=True)

    def browse_ring_exe():
        filename = QtWidgets.QFileDialog.getOpenFileNames(dialog, "Select Ring executable")[0][0]

        if filename:
            form.ring_path.setText(filename)

    def refresh_sele():
        form.sele_names.blockSignals(True)
        selections = cmd.get_names('public_selections')
        selections.extend(list(filter(lambda x: x.split('_')[-1][-3:] != 'cgo',
                                      cmd.get_names('public_nongroup_objects'))))
        selections = sorted(selections)
        if selections != sorted([form.sele_names.itemText(i) for i in range(form.sele_names.count())]):
            form.sele_names.clear()
            form.sele_names.addItems(selections)
        form.sele_names.blockSignals(False)

    def refresh_obj():
        form.pymol_obj.blockSignals(True)
        selections = sorted(list(filter(lambda x: x.split('_')[-1][-3:] != 'cgo',
                                        cmd.get_names('public_nongroup_objects'))))
        if selections != sorted([form.pymol_obj.itemText(i) for i in range(form.pymol_obj.count())]):
            form.pymol_obj.clear()
            form.pymol_obj.addItems(selections)
            form.sele_obj.clear()
            form.sele_obj.addItems(selections)
            form.sele_obj_2.clear()
            form.sele_obj_2.addItems(selections)
        form.pymol_obj.blockSignals(False)

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
                log("Please run Ring on the selected object first!", error=True)
                return None
        return conn_freq

    def get_freq_combined(obj, bond):
        import math
        conn_freq = dict()
        try:
            with open("/tmp/ring/md/{}.gfreq_{}".format(obj, bond), 'r') as f:
                for line in f:
                    edge1, _, edge2, perc = line.split('\t')
                    edge1 = edge1.replace(":_:", ":")
                    edge2 = edge2.replace(":_:", ":")
                    conn_freq.setdefault(edge1, [])
                    chain1 = edge1.split(':')[0]
                    chain2 = edge2.split(':')[0]
                    if form.interchain_2.isChecked() and chain1 == chain2:
                        continue
                    if form.intrachain_2.isChecked() and chain1 != chain2:
                        continue
                    conn_freq[edge1].append(float(perc))
        except FileNotFoundError:
            log("Please run Ring on the selected object first!", error=True)
            return None
        for k, v in conn_freq.items():
            conn_freq[k] = 1 - math.prod([(1 - x) for x in v])

        return conn_freq

    def visualize(first=False, selection=None, color=None, int_type=None, pair_set=None, block=True, log_iter=False):
        from pymol import stored

        if block:
            form.main.setEnabled(False)
            app.processEvents()
        if selection:
            obj = selection
        else:
            if first:
                obj = form.pymol_obj.currentText()
            else:
                obj = form.sele_names.currentText()

        if obj == '':
            log("Please provide a selection", error=True)
            form.main.setEnabled(True)
            return

        cmd.delete(obj + "_edges")

        states = int(cmd.count_states(selection=obj))
        stored.model = ""
        cmd.iterate(obj, 'stored.model = model')
        file_pth = "/tmp/ring/" + stored.model + ".cif_ringEdges"

        if not os.path.exists(file_pth):
            log("Please run Ring on the selected object first!", error=True)
            return

        interactions_per_state = pd.read_csv(file_pth, sep='\t')

        for state in range(1, states + 1):
            df = interactions_per_state[interactions_per_state.Model == state]

            stored.chain_resi = set()

            conn_freq = get_freq(stored.model)

            if conn_freq is None:
                form.main.setEnabled(True)
                app.processEvents()
                return

            cmd.iterate(obj, 'stored.chain_resi.add((chain, resi))')
            stored.coords = dict()
            cmd.iterate_state(state=state, selection=obj,
                              expression='stored.coords["chain {} and resi {} and name {}".format(chain, resi, name)] = [x,y,z]')

            interactions_per_type = dict()

            for (nodeId1, interaction, nodeId2, _, _, _, atom1, atom2, *_) in df.itertuples(index=False):

                intType, intSubType = interaction.split(":")
                chain1, pos1, i1, res1 = nodeId1.split(":")
                chain2, pos2, i2, res2 = nodeId2.split(":")
                if chain1 == '.':
                    chain1 = ''
                if chain2 == '.':
                    chain2 = ''

                try:
                    freq = conn_freq[intType][("{}:{}".format(chain1, pos1), "{}:{}".format(chain2, pos2))] * 100
                except KeyError:
                    freq = 0.5 * 100

                if not selection:
                    if form.interchain_1.isChecked() and chain1 == chain2:
                        continue
                    if form.intrachain_1.isChecked() and chain1 != chain2:
                        continue
                tmp1 = ("{}:{}".format(chain1, pos1), "{}:{}".format(chain2, pos2))
                tmp2 = ("{}:{}".format(chain2, pos2), "{}:{}".format(chain1, pos1))
                if (chain1, pos1) in stored.chain_resi and (chain2, pos2) in stored.chain_resi \
                        and (form.min_freq.value() <= freq <= form.max_freq.value() or selection) \
                        and ((selection and int_type == intType and
                              (tmp1 in pair_set or tmp2 in pair_set)) or not selection):
                    interactions_per_type.setdefault(intType, [])

                    if intType == "PIPISTACK" or intType == "IONIC":
                        t = tuple()
                        if "," in atom1:
                            t += (atom1,)
                        else:
                            t += ("chain {} and resi {} and name {}".format(chain1, pos1, atom1),)
                        if "," in atom2:
                            t += (atom2,)
                        else:
                            t += ("chain {} and resi {} and name {}".format(chain2, pos2, atom2),)
                        interactions_per_type[intType].append(t)
                    else:
                        interactions_per_type[intType].append(
                                ("chain {} and resi {} and name {}".format(chain1, pos1, atom1),
                                 "chain {} and resi {} and name {}".format(chain2, pos2, atom2)))

            not_present = 0
            for intType, interactions in interactions_per_type.items():
                not_present += draw_links(interactions,
                                          object_name=obj + "_" + intType + "_cgo" if not selection else selection + "_cgo",
                                          color=intTypeMap[intType] if not color else color,
                                          coords=stored.coords,
                                          state=state)

            if log_iter:
                log_s = "Interactions state {}: ".format(state)
                for intType in sorted(interactions_per_type.keys()):
                    log_s += "{} {}, ".format(intType, len(interactions_per_type[intType]))
                log(log_s.rstrip(', '), timed=False, process=False)

            if not_present > 0:
                log("{} connections not displayed because atoms not present".format(not_present), warning=True,
                    process=False)

        if not selection:
            if form.check_hide_others.isChecked():
                cmd.hide(selection="*_edges")

            members = ""
            for k in intTypeMap.keys():
                members += " {}_{}_cgo".format(obj, k)
            cmd.group(obj + "_edges", members=members)
        else:
            cmd.hide(selection="*_edges")

        if block:
            form.main.setEnabled(True)
            app.processEvents()

    def inter_freq_analysis():
        obj = form.sele_obj.currentText()
        if obj == '':
            log("Please provide a selection", error=True)
            return

        disable_window()

        inter = ""
        if form.hbond.isChecked():
            inter = "HBOND"
        if form.ionic.isChecked():
            inter = "IONIC"
        if form.pipistack.isChecked():
            inter = "PIPISTACK"
        if form.pication.isChecked():
            inter = "PICATION"
        if form.vdw.isChecked():
            inter = "VDW"
        if form.ssbond.isChecked():
            inter = "SSBOND"
        if form.iac.isChecked():
            inter = "IAC"
        conn_freq = get_freq_combined(obj, inter)

        cmd.delete("filtered_resi")
        for k, v in conn_freq.items():
            if form.min_freq_2.value() / 100 <= v <= form.max_freq_2.value() / 100:
                cmd.select(name="filtered_resi", selection="/{}//{}/{}/".format(obj, k.split(':')[0], k.split(':')[1]),
                           merge=1)
            else:
                conn_freq[k] = 0.001

        if conn_freq is not None:
            myspace = {'dict_freq': conn_freq}
            express = "b=dict_freq['{}:{}:{}'.format(chain,resi,resn)] if '{}:{}:{}'.format(chain,resi,resn) " \
                      "in dict_freq.keys() else 0.001"
            cmd.alter_state(-1, obj, expression=express, space=myspace)
            cmd.spectrum("b", "white yellow orange red", obj, minimum=0.001, maximum=1.0)

        log("Selection of filtered residues based on bond type and frequency created")
        enable_window()

    def inter_freq_table():
        import numpy as np

        obj = form.sele_obj.currentText()
        if obj == '':
            log("Please provide a selection", error=True)
            return
        disable_window()
        log("Creation of residue interaction frequency table")
        freq_bond = dict()
        for bondType in intTypeMap.keys():
            freq_bond.setdefault(bondType, get_freq_combined(obj, bondType))

        all_resi = set()
        for bondType, freqs in freq_bond.items():
            all_resi.update(set(freqs.keys()))

        mtr = np.zeros((len(all_resi), 7))
        for i, resi in enumerate(all_resi):
            for j, bondType in enumerate(["PIPISTACK", "PICATION", "IONIC", "HBOND", "SSBOND", "VDW", "IAC"]):
                try:
                    mtr[i, j] = freq_bond[bondType][resi] * 100
                except KeyError:
                    pass

        right = pd.DataFrame(mtr, columns=["PIPISTACK", "PICATION", "IONIC", "HBOND", "SSBOND", "VDW", "IAC"])
        resi_chain = [x.split(':')[0] for x in all_resi]
        resi_number = [int(x.split(':')[1]) for x in all_resi]
        center = pd.DataFrame([all_resi]).transpose()
        left = pd.DataFrame([resi_chain, resi_number]).transpose().rename(columns={0: 'chain', 1: 'number'})
        df = pd.concat([left, center, right], axis=1)
        tableWidget = dialog_freq.freqTable
        tableWidget.setRowCount(0)

        tableWidget.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)

        df = df.sort_values(by=['chain', 'number'], ascending=[True, True]).reset_index()
        prevEdge = None
        color = 2
        for (rowPosition, (_, row)) in enumerate(df.iterrows()):
            tableWidget.insertRow(rowPosition)
            for i, item in enumerate(row.to_list()[3:]):
                if i == 0 and item != prevEdge:
                    if color == 1:
                        bk_color = QColor(173, 173, 173)
                        fg_color = QColor(0, 0, 0)
                        color = 2
                    else:
                        bk_color = QColor(121, 121, 121)
                        fg_color = QColor(255, 255, 255)
                        color = 1

                if i == 0:
                    tableWidget.setItem(rowPosition, i, QtWidgets.QTableWidgetItem(item))
                else:
                    wItem = QtWidgets.QTableWidgetItem()
                    wItem.setData(QtCore.Qt.DisplayRole, float(item))
                    tableWidget.setItem(rowPosition, i, wItem)
                tableWidget.item(rowPosition, i).setBackground(bk_color)
                tableWidget.item(rowPosition, i).setForeground(fg_color)
                tableWidget.item(rowPosition, i).setTextAlignment(QtCore.Qt.AlignCenter)
                if i == 0:
                    prevEdge = item
        tableWidget.viewport().update()
        dialog_freq.show()
        enable_window()

    def slider_radius_change():
        value = form.radius_value.value()
        cmd.set("cgo_line_width", value)

    def slider_transp_change():
        value = form.transp_value.value() / 100
        cmd.set("cgo_transparency", value)

    def correlation_obj():
        obj = form.sele_obj_2.currentText()
        if obj == '':
            log("Please select an object first!", error=True)
            return
        disable_window()
        log("Calculation of correlated interactions started")

        inter = ""
        if form.hbond_2.isChecked():
            inter = "HBOND"
        if form.ionic_2.isChecked():
            inter = "IONIC"
        if form.pipistack_2.isChecked():
            inter = "PIPISTACK"
        if form.pication_2.isChecked():
            inter = "PICATION"
        if form.vdw_2.isChecked():
            inter = "VDW"
        if form.ssbond_2.isChecked():
            inter = "SSBOND"
        if form.iac_2.isChecked():
            inter = "IAC"

        coeff_thr = float(form.coeff_thr.value())
        p_thr = float(form.p_thr.value())

        max_presence = float(form.max_presence.value()) / 100
        min_presence = float(form.min_presence.value()) / 100

        states = cmd.count_states(obj)
        if states < 2:
            log("Correlation cannot be calculated if the number of states is {}".format(states), error=True)
            enable_window()
            return
        try:
            selections, corr_matr, p_matr = calculate_correlation(obj, states, int_type=inter,
                                                                  coeff_thresh=coeff_thr,
                                                                  p_thresh=p_thr, max_presence=max_presence,
                                                                  min_presence=min_presence)
        except TypeError:
            log("Run ring on all the states first!", error=True)
            enable_window()
            return

        correlations.setdefault(obj, dict())
        correlations[obj][inter] = (selections, corr_matr, p_matr)

        create_table(obj, inter)

        enable_window()

    def create_table(obj, inter):
        import numpy as np
        import pandas as pd

        selections, corr_matr, p_matr = correlations[obj][inter]
        with np.errstate(divide='ignore', invalid='ignore'):
            indexes = np.argwhere(~np.isnan(corr_matr))

        edge1s = [selections[x] for x in [y[0] for y in indexes]]
        edge2s = [selections[x] for x in [y[1] for y in indexes]]
        corr_vals = [corr_matr[i, j] for (i, j) in indexes]
        p_vals = [p_matr[i, j] for (i, j) in indexes]

        tableWidget = dialog_corr.corrTable
        tableWidget.setRowCount(0)
        rowPosition = tableWidget.rowCount()  # necessary even when there are no rows in the table
        prevEdge = None
        color = 2

        edge1_chains1 = [x.split(':')[0] for x in edge1s]
        edge1_resi1 = [int(x.split(':')[1]) for x in edge1s]
        edge1_chains2 = [x.split('- ')[1].split(':')[0] for x in edge1s]
        edge1_resi2 = [int(x.split(':')[3]) for x in edge1s]
        df = pd.DataFrame(
                [edge1s, edge2s, corr_vals, p_vals, edge1_chains1, edge1_resi1, edge1_chains2, edge1_resi2]).transpose()
        df = df.sort_values([4, 5, 6, 7, 2], ascending=(True, True, True, True, False))

        log("{} {} interactions correlates/anti-correlates".format(len(df), inter))

        tableWidget.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        for _, row in df.iterrows():
            x, y, z, w, *_ = row.to_list()
            if x != prevEdge:
                if color == 1:
                    bk_color = QColor(173, 173, 173)
                    fg_color = QColor(0, 0, 0)
                    color = 2
                else:
                    bk_color = QColor(121, 121, 121)
                    fg_color = QColor(255, 255, 255)
                    color = 1

            tableWidget.insertRow(rowPosition)

            tableWidget.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(x))
            tableWidget.setItem(rowPosition, 1, QtWidgets.QTableWidgetItem(y))

            zItem = QtWidgets.QTableWidgetItem()
            zItem.setData(QtCore.Qt.EditRole, float(z))
            tableWidget.setItem(rowPosition, 2, zItem)

            wItem = QtWidgets.QTableWidgetItem()
            wItem.setData(QtCore.Qt.EditRole, float(w))
            tableWidget.setItem(rowPosition, 3, wItem)
            for i in range(4):
                tableWidget.item(rowPosition, i).setBackground(bk_color)
                tableWidget.item(rowPosition, i).setForeground(fg_color)
                tableWidget.item(rowPosition, i).setTextAlignment(QtCore.Qt.AlignCenter)

            if z > 0:
                tableWidget.item(rowPosition, 2).setForeground(QColor(0, 0, 255))
            else:
                tableWidget.item(rowPosition, 2).setForeground(QColor(255, 0, 0))

            prevEdge = x
            rowPosition += 1
        tableWidget.viewport().update()
        dialog_corr.show()

    def show_corr():
        obj = form.sele_obj_2.currentText()
        if obj == '':
            log("Please select an object first!", error=True)
            return
        disable_window()

        cmd.delete("edge1 edge2 edge1_cgo edge2_cgo")
        inter = ""
        if form.hbond_2.isChecked():
            inter = "HBOND"
        if form.ionic_2.isChecked():
            inter = "IONIC"
        if form.pipistack_2.isChecked():
            inter = "PIPISTACK"
        if form.pication_2.isChecked():
            inter = "PICATION"
        if form.vdw_2.isChecked():
            inter = "VDW"
        if form.ssbond_2.isChecked():
            inter = "SSBOND"
        if form.iac_2.isChecked():
            inter = "IAC"

        table = dialog_corr.corrTable
        selection = table.selectionModel().selectedRows()

        sele_set = set()
        corr_set = set()
        anti_set = set()
        for i in range(len(selection)):
            row = selection[i].row()
            edge1 = table.item(row, 0).text()
            edge2 = table.item(row, 1).text()
            corr_val = float(table.item(row, 2).text())

            resi1, resi2 = edge1.split(' - ')
            resi1, resi2 = resi1[:-4], resi2[:-4]
            resi1_c, resi1_n = resi1.split(':')
            resi2_c, resi2_n = resi2.split(':')

            resi3, resi4 = edge2.split(' - ')
            resi3, resi4 = resi3[:-4], resi4[:-4]
            resi3_c, resi3_n = resi3.split(':')
            resi4_c, resi4_n = resi4.split(':')

            cmd.select("edge1", "/{}//{}/{} or /{}//{}/{}".format(obj, resi1_c, resi1_n,
                                                                  obj, resi2_c, resi2_n), merge=1)
            sele_set.add((resi1, resi2))

            cmd.select("edge2", "/{}//{}/{} or /{}//{}/{}".format(obj, resi3_c, resi3_n,
                                                                  obj, resi4_c, resi4_n), merge=1)
            if corr_val > 0:
                corr_set.add((resi3, resi4))
            else:
                anti_set.add((resi3, resi4))

        visualize(selection="edge1", color="white", int_type=inter, pair_set=sele_set, block=False)
        visualize(selection="edge2", color="blue", int_type=inter, pair_set=corr_set, block=False)
        visualize(selection="edge2", color="red", int_type=inter, pair_set=anti_set, block=False)
        log("Selection edge1 contains all the residues from the edges selected in the first column", timed=False)
        log("Selection edge2 contains all the residues from the edges selected in the second column", timed=False)
        log("CGO objects edge1_cgo and edge2_cgo are the selected edges from the first and second column respectively",
            timed=False)
        log("Interactions in blue are the one correlating, and in red the ones that anti-correlates with the respective white interactions",
            timed=False)
        enable_window()

    def tab_click_update():
        refresh_obj()
        refresh_sele()

    form.ring_path.setText("{}/.ring/bin/Ring-md".format(environ['HOME']))

    # Execute Ring
    form.button_start.clicked.connect(run)
    form.ring_exec_button.clicked.connect(browse_ring_exe)

    # Update view
    form.visualize_btn.clicked.connect(lambda: visualize(log_iter=True))
    form.radius_value.valueChanged.connect(slider_radius_change)
    form.transp_value.valueChanged.connect(slider_transp_change)

    # Putty repr
    form.show_resi_analysis.clicked.connect(inter_freq_analysis)
    form.show_resi_analysis_table.clicked.connect(inter_freq_table)

    form.timer = QtCore.QTimer()
    form.timer.timeout.connect(tab_click_update)
    form.timer.start(1500)

    # Correlation
    form.calc_corr.clicked.connect(correlation_obj)
    dialog_corr.show_corr.clicked.connect(show_corr)

    return dialog
