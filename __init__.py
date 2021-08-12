# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.
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
    # import seaborn as sn
    # from matplotlib import pyplot as plt
    try:
        import pandas as pd
        import numpy as np
        from scipy.stats import pearsonr
    except ImportError:
        print("To run this you have to install pandas, numpy and scipy in python")
        return

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
            int_id1 = names[i[0]].split(':')
            int_id1 = (int_id1[0], int(int_id1[1]), int_id1[2])
            int_id2 = names[i[1]].split(':')
            int_id2 = (int_id2[0], int(int_id2[1]), int_id2[2])
            tmp = tuple(sorted([int_id1, int_id2]))
            contacts_sparse.setdefault(tmp, [])
            contacts_sparse[tmp].append(j - 1)

    to_pop = []
    for k, v in contacts_sparse.items():
        presence = len(v) / frames
        if max_presence < presence or presence < min_presence:
            to_pop.append(k)

    for k in to_pop:
        contacts_sparse.pop(k)

    z = np.zeros((len(contacts_sparse), frames))
    for i, v in enumerate(contacts_sparse.values()):
        for j in v:
            z[i, j] = 1

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


# global reference to avoid garbage collection of our dialog
dialog = None
dialog_corr = None
dialog_freq = None
correlations = dict()
prev_launch_config = dict()
prev_sele = ""


def run_plugin_gui():
    """
    Open our custom dialog
    """
    global dialog
    global dialog_corr
    global dialog_freq

    if dialog is None and dialog_corr is None and dialog_freq is None:
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
        now = "{}:{:02}:{:02}".format(now.hour, now.minute, now.second)
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

    try:
        import numpy as np
        import pandas as pd
    except ImportError:
        log("Please install numpy and pandas to use this plugin")

    def run():
        from pymol import stored

        import subprocess

        obj_name = form.selections_list.currentText()

        # If it is a selection then try to just visualize
        if obj_name[0] == "(" and obj_name[-1] == ")":
            visualize(log_iter=True)
            return

        obj_name = obj_name.lstrip('(').rstrip(')')

        if not obj_name:
            log("Please select a Pymol object first!", error=True)
            return
        if not os.path.exists('/tmp/ring'):
            os.mkdir('/tmp/ring')

        stored.chains = ""
        cmd.iterate(obj_name, 'stored.chains += chain')
        if stored.chains == "":
            log("Pymol object does not contain chain name, please set it before running Ring. (Use alter)", error=True)
            return

        file_pth = "/tmp/ring/" + obj_name + ".cif"

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

        current_run_config = {"-g"   : seq_sep,
                              "-o"   : len_salt,
                              "-s"   : len_ss,
                              "-k"   : len_pipi,
                              "-a"   : len_pica,
                              "-b"   : len_hbond,
                              "-w"   : len_vdw,
                              "edges": edge_policy}

        if obj_name in prev_launch_config.keys() and prev_launch_config[obj_name] == current_run_config:
            visualize(log_iter=True)
            return

        form.main.setEnabled(False)
        form.visualize_btn.setText("Running ring...")
        app.processEvents()

        log("Exporting pymol object {} in cif format ({})".format(obj_name, file_pth))
        cmd.save(filename=file_pth, selection=obj_name, state=0)
        log("Exporting done")

        try:
            p = subprocess.Popen(
                    [form.ring_path.text(), "-i", file_pth, "--out_dir", "/tmp/ring/", "-g", seq_sep,
                     "-o", len_salt, "-s", len_ss, "-k", len_pipi, "-a", len_pica, "-b", len_hbond, "-w", len_vdw,
                     "--all_chains", edge_policy, "--all_models"], stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE, universal_newlines=True)

            prev_launch_config[obj_name] = current_run_config
        except FileNotFoundError:
            log("Ring path is not correct!", error=True)
            form.visualize_btn.setText("Show")
            form.main.setEnabled(True)
            return

        log("Ring generation started")

        while p.poll() is None:
            line = p.stderr.readline()
            if line != "":
                if "model" in line:
                    nModel = int(line.split("model ")[1].strip())
                    log("Running on model {} | {:.2%}".format(nModel, nModel / nStates))
            app.processEvents()

        form.visualize_btn.setText("Show")
        log("Ring generation finished")
        visualize(log_iter=True)

    def browse_ring_exe():
        filename = QtWidgets.QFileDialog.getOpenFileNames(dialog, "Select Ring executable")[0][0]

        if filename:
            form.ring_path.setText(filename)

    def refresh_sele():
        form.selections_list.blockSignals(True)
        present = set(form.selections_list.itemText(i) for i in range(form.selections_list.count()))
        selections = set(map(lambda x: "(" + x + ")", cmd.get_names('public_selections')))
        selections.update(list(filter(lambda x: x.split('_')[-1][-3:] != 'cgo',
                                      cmd.get_names('public_nongroup_objects'))))
        not_present_anymore = present - selections
        new_selections = selections - present

        for sele in not_present_anymore:
            for sele_in, idx in [(form.selections_list.itemText(i), i) for i in range(form.selections_list.count())]:
                if sele == sele_in:
                    form.selections_list.removeItem(idx)
                    break

        for sele in new_selections:
            form.selections_list.addItem(sele, form.selections_list.count())

        form.selections_list.blockSignals(False)
        global prev_sele

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
                    chain1 = edge1.split(':')[0]
                    chain2 = edge2.split(':')[0]
                    if form.interchain.isChecked() and chain1 == chain2:
                        continue
                    if form.intrachain.isChecked() and chain1 != chain2:
                        continue
                    conn_freq.setdefault(edge1, [])
                    conn_freq[edge1].append(float(perc))
        except FileNotFoundError:
            log("Please run Ring on the selected object first!", error=True)
            return None
        for k, v in conn_freq.items():
            conn_freq[k] = 1 - math.prod([(1 - x) for x in v])

        return conn_freq

    def visualize(selection=None, color=None, int_type=None, pair_set=None, block=True, log_iter=False):
        from pymol import stored

        if block:
            form.main.setEnabled(False)
            app.processEvents()
        if selection:
            obj = selection
        else:
            obj = form.selections_list.currentText()

        if obj == '':
            log("Please provide a selection", error=True)
            form.main.setEnabled(True)
            return

        is_selection = obj[0] == "(" and obj[-1] == ")"

        obj = obj.lstrip('(').rstrip(')')

        cmd.delete(obj + "_edges")
        cmd.delete(obj + "_nodes")

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

                try:
                    freq = conn_freq[intType][("{}:{}".format(chain1, pos1), "{}:{}".format(chain2, pos2))] * 100
                except KeyError:
                    freq = 0.5 * 100

                if not selection:
                    if form.interchain.isChecked() and chain1 == chain2:
                        continue
                    if form.intrachain.isChecked() and chain1 != chain2:
                        continue
                tmp1 = ("{}/{}".format(chain1, pos1), "{}/{}".format(chain2, pos2))
                tmp2 = ("{}/{}".format(chain2, pos2), "{}/{}".format(chain1, pos1))
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
                other = ""
                for intType in sorted(interactions_per_type.keys()):
                    other += "{} {}, ".format(intType, len(interactions_per_type[intType]))
                log_s += other if len(other) > 1 else "No interaction for this state, maybe check the filters?"
                log(log_s.rstrip(', '), timed=False, process=False, warning=len(other) == 0)

            if not_present > 0:
                log("{} connections not displayed because atoms not present".format(not_present), warning=True,
                    process=False)

        cmd.hide(selection="*_edges")
        if not selection:
            members = ""
            for k in intTypeMap.keys():
                members += " {}_{}_cgo".format(obj, k)
            cmd.group(obj + "_edges", members=members)
            log("Created group {} for interaction edges".format(obj + "_edges"), timed=False)

            if not is_selection:
                members = ""
                for bond in intTypeMap.keys():
                    sele = "{}_{}_resi".format(obj, bond)
                    members += " {}".format(sele)
                    freqs = get_freq_combined(stored.model, bond)
                    for edge, freq in freqs.items():
                        if form.min_freq.value() <= freq * 100 <= form.max_freq.value():
                            cmd.select(sele,
                                       selection="chain {} and resi {}".format(edge.split(':')[0], edge.split(':')[1]),
                                       merge=1)
                cmd.group(obj + "_nodes", members=members)
                log("Created group {} for interaction nodes".format(obj + "_edges"), timed=False)

        # Set transp and radius after updating the CGOs
        slider_radius_change()
        slider_transp_change()
        if block:
            form.main.setEnabled(True)
            app.processEvents()

    def inter_freq_analysis():
        obj = form.selections_list.currentText().lstrip('(').rstrip(')')
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

        if conn_freq is not None:
            myspace = {'dict_freq': conn_freq}
            express = "b=dict_freq['{}:{}:{}'.format(chain,resi,resn)] if '{}:{}:{}'.format(chain,resi,resn) " \
                      "in dict_freq.keys() else 0.001"
            cmd.alter_state(-1, obj, expression=express, space=myspace)
            cmd.spectrum("b", "white yellow orange red", obj, minimum=0.001, maximum=1.0)

        log("Residues with interaction of type {} colored based on the frequency of contact".format(inter))
        enable_window()

    def inter_freq_table():
        import numpy as np

        obj = form.selections_list.currentText().lstrip('(').rstrip(')')
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
        tableWidget.setSortingEnabled(False)

        tableWidget.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)

        df = df.sort_values(by=['chain', 'number'], ascending=[True, True]).reset_index()
        prevEdge = None
        color = 2
        for (rowPosition, (_, row)) in enumerate(df.iterrows()):
            tableWidget.insertRow(rowPosition)
            for i, item in enumerate(row.to_list()[3:]):
                if i == 0 and item != prevEdge:
                    if color == 1:
                        bk_color = QColor(232, 231, 252)
                        fg_color = QColor(0, 0, 0)
                        color = 2
                    else:
                        bk_color = QColor(255, 255, 255)
                        fg_color = QColor(0, 0, 0)
                        color = 1

                if i == 0:
                    tableWidget.setItem(rowPosition, i, QtWidgets.QTableWidgetItem(item))
                else:
                    wItem = QtWidgets.QTableWidgetItem()
                    wItem.setData(QtCore.Qt.DisplayRole, round(float(item), 2))
                    tableWidget.setItem(rowPosition, i, wItem)
                tableWidget.item(rowPosition, i).setBackground(bk_color)
                tableWidget.item(rowPosition, i).setForeground(fg_color)
                tableWidget.item(rowPosition, i).setTextAlignment(QtCore.Qt.AlignCenter)
                if i == 0:
                    prevEdge = item
        tableWidget.setSortingEnabled(True)
        tableWidget.viewport().update()
        dialog_freq.show()
        enable_window()

    def sele_selected_resi_freq_table():
        tableWidget = dialog_freq.freqTable
        indexes = tableWidget.selectionModel().selectedRows()
        for index in sorted(indexes):
            chain, resi = tableWidget.item(index.row(), 0).text().split(":")[0:2]
            cmd.select("sele_row", selection="chain {} and resi {}".format(chain, resi), merge=1)
            log("Updated selection sele_row with the residue selected in the frequency table")

    def slider_radius_change():
        value = form.radius_value.value()
        cmd.set("cgo_line_width", value)

    def slider_transp_change():
        value = form.transp_value.value() / 100
        cmd.set("cgo_transparency", value)

    def correlation_obj():
        obj = form.selections_list.currentText().lstrip('(').rstrip(')')
        if obj == '':
            log("Please select an object first!", error=True)
            return
        disable_window()
        log("Calculation of correlated interactions started")

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
            for inter in intTypeMap.keys():
                selections, corr_matr, p_matr = calculate_correlation(obj, states, int_type=inter,
                                                                      coeff_thresh=coeff_thr,
                                                                      p_thresh=p_thr, max_presence=max_presence,
                                                                      min_presence=min_presence)
                correlations.setdefault(obj, dict())
                correlations[obj][inter] = (selections, corr_matr, p_matr)

        except TypeError:
            log("Run ring on all the states first!", error=True)
            enable_window()
            return

        create_table(obj)

        enable_window()

    def create_table(obj):
        df = pd.DataFrame()

        for inter in intTypeMap.keys():
            selections, corr_matr, p_matr = correlations[obj][inter]
            with np.errstate(divide='ignore', invalid='ignore'):
                indexes = np.argwhere(~np.isnan(corr_matr))

            edge1s = [selections[x] for x in [y[0] for y in indexes]]
            edge2s = [selections[x] for x in [y[1] for y in indexes]]
            inter_labels = [inter for x in edge1s]
            corr_vals = [corr_matr[i, j] for (i, j) in indexes]
            p_vals = [p_matr[i, j] for (i, j) in indexes]

            tableWidget = dialog_corr.corrTable
            tableWidget.setRowCount(0)
            tableWidget.setSortingEnabled(False)

            rowPosition = tableWidget.rowCount()  # necessary even when there are no rows in the table
            prevEdge = None
            color = 2

            edge1_chains1 = [x.split('/')[0] for x in edge1s]
            edge1_resi1 = [int(x.split('/')[1]) for x in edge1s]
            edge1_chains2 = [x.split('- ')[1].split('/')[0] for x in edge1s]
            edge1_resi2 = [int(x.split('/')[3]) for x in edge1s]
            df = df.append(pd.DataFrame([edge1s, inter_labels, edge2s, corr_vals, p_vals, edge1_chains1,
                                         edge1_resi1, edge1_chains2, edge1_resi2]).transpose())

        df = df.sort_values([5, 6, 7, 8, 3], ascending=(True, True, True, True, False))

        log("{} {} interactions correlates/anti-correlates".format(len(df), inter))

        tableWidget.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        for _, row in df.iterrows():
            x, p, y, z, w, *_ = row.to_list()
            if x != prevEdge:
                if color == 1:
                    bk_color = QColor(232, 231, 252)
                    fg_color = QColor(0, 0, 0)
                    color = 2
                else:
                    bk_color = QColor(255, 255, 255)
                    fg_color = QColor(0, 0, 0)
                    color = 1

            tableWidget.insertRow(rowPosition)

            tableWidget.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(x))
            tableWidget.setItem(rowPosition, 1, QtWidgets.QTableWidgetItem(p))
            tableWidget.setItem(rowPosition, 2, QtWidgets.QTableWidgetItem(y))

            zItem = QtWidgets.QTableWidgetItem()
            zItem.setData(QtCore.Qt.EditRole, round(float(z), 2))
            tableWidget.setItem(rowPosition, 3, zItem)

            wItem = QtWidgets.QTableWidgetItem()
            wItem.setData(QtCore.Qt.EditRole, float(w))
            tableWidget.setItem(rowPosition, 4, wItem)
            for i in range(5):
                tableWidget.item(rowPosition, i).setBackground(bk_color)
                tableWidget.item(rowPosition, i).setForeground(fg_color)
                tableWidget.item(rowPosition, i).setTextAlignment(QtCore.Qt.AlignCenter)

            if z > 0:
                tableWidget.item(rowPosition, 3).setForeground(QColor(0, 0, 255))
            else:
                tableWidget.item(rowPosition, 3).setForeground(QColor(255, 0, 0))

            prevEdge = x
            rowPosition += 1
        tableWidget.setSortingEnabled(True)
        tableWidget.viewport().update()
        dialog_corr.show()

    def show_corr():
        obj = form.selections_list.currentText().lstrip('(').rstrip(')')
        if obj == '':
            log("Please select at least one row first!", error=True)
            return
        disable_window()

        cmd.delete("edge1 edge2 edge1_cgo edge2_cgo")

        table = dialog_corr.corrTable
        selection = table.selectionModel().selectedRows()

        sele_set = set()
        corr_set = set()
        anti_set = set()
        for i in range(len(selection)):
            row = selection[i].row()
            edge1 = table.item(row, 0).text()
            inter = table.item(row, 1).text()
            edge2 = table.item(row, 2).text()
            corr_val = float(table.item(row, 3).text())

            resi1, resi2 = edge1.split(' - ')
            resi1, resi2 = resi1[:-4], resi2[:-4]
            resi1_c, resi1_n = resi1.split('/')
            resi2_c, resi2_n = resi2.split('/')

            resi3, resi4 = edge2.split(' - ')
            resi3, resi4 = resi3[:-4], resi4[:-4]
            resi3_c, resi3_n = resi3.split('/')
            resi4_c, resi4_n = resi4.split('/')

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

    def resi_plot():
        from pymol import stored

        try:
            from matplotlib import pyplot as plt
        except ImportError:
            log("To run this feature you have to install matplotlib for the python version that PyMol is using",
                error=True)
            return

        obj = form.selections_list.currentText()
        if len(obj) == 0 or obj[0] != "(" or obj[-1] != ")":
            log("Please select a selection on the box above to use this feature", error=True)
            return

        states = int(cmd.count_states(selection=obj))
        if states == 1:
            log("You can use this feature only on multi-state objects", error=True)
            return

        stored.model = ""
        cmd.iterate(obj, 'stored.model = model')
        file_pth = "/tmp/ring/" + stored.model + ".cif_ringEdges"

        if not os.path.exists(file_pth):
            log("Before this you need to run Ring-md on the whole object first. Select it above and press the Show button",
                error=True)
            return

        stored.chain_resi = set()
        cmd.iterate(obj, 'stored.chain_resi.add((chain, int(resi), resn))')
        if len(stored.chain_resi) != 2:
            log("You need to create a selection with exactly two residues to use this feature", error=True)
            return

        resi1, resi2 = list(stored.chain_resi)
        resi1_name = resi1[0] + "/" + str(resi1[1]) + "/" + resi1[2]
        resi1 = (resi1[0], resi1[1])
        resi2_name = resi2[0] + "/" + str(resi2[1]) + "/" + resi2[2]
        resi2 = (resi2[0], resi2[1])
        interaction_distance = dict()

        interactions_per_state = pd.read_csv(file_pth, sep='\t')
        for state in range(1, states + 1):
            df = interactions_per_state[interactions_per_state.Model == state]
            for (nodeId1, interaction, nodeId2, distance, _, _, atom1, atom2, *_) in df.itertuples(index=False):

                intType, intSubType = interaction.split(":")
                chain1, pos1, *_ = nodeId1.split(":")
                chain2, pos2, *_ = nodeId2.split(":")
                resi11 = (chain1, int(pos1))
                resi22 = (chain2, int(pos2))

                if (resi1 == resi11 and resi2 == resi22) or (resi2 == resi11 and resi1 == resi22):
                    interaction_distance.setdefault(intType, np.ones(states) * 999)
                    interaction_distance[intType][state - 1] = min(interaction_distance[intType][state - 1],
                                                                   float(distance))

        for inter in interaction_distance.keys():
            interaction_distance[inter] = list(map(lambda x: x if x != 999 else np.nan, interaction_distance[inter]))

        plt.close()
        something = False
        for inter in interaction_distance.keys():
            something = True
            plt.scatter(np.arange(1, states + 1), interaction_distance[inter], label=inter, c=intColorMap[inter],
                        marker='.')

        if something:
            plt.title("{} - {}".format(resi1_name, resi2_name))
            plt.grid()
            plt.ylim(bottom=0, top=max(filter(lambda x: x != np.nan,
                                              [item for sublist in list(interaction_distance.values()) for item in
                                               sublist])) + 1)
            plt.xlim(left=0, right=states + 1)
            plt.xlabel("State")
            plt.ylabel("Distance (Å)")
            plt.legend()
            plt.tight_layout()
            plt.show()
        else:
            log("No interactions found between the two selected residues", warning=True)

    def chain_graph():
        from pymol import stored
        try:
            from matplotlib import pyplot as plt
            import matplotlib.patches as mpatches
            import networkx as nx
        except ImportError:
            log("To run this feature you have to install matplotlib and networkx and for the python version "
                "that PyMol is using", error=True)
            return

        obj = form.selections_list.currentText()
        if len(obj) == 0:
            log("Please select an object on the box above to use this feature", error=True)
            return

        if obj[0] == "(" and obj[-1] == ")":
            log("Please select an object on the box above to use this feature", error=True)
            return

        plt.close()
        stored.model = ""
        cmd.iterate(obj, 'stored.model = model')
        file_pth = "/tmp/ring/" + stored.model + ".cif_ringEdges"
        if not os.path.exists(file_pth):
            log("Before this you need to run Ring-md on the object first. Select it above and press the Show button",
                error=True)
            return

        if cmd.get_chains(stored.model) == 1:
            log("Only one chain is present in the selected object", warning=True)

        G = nx.MultiGraph()
        edges = dict()
        interactions_first_state = pd.read_csv(file_pth, sep='\t')
        interactions_first_state = interactions_first_state[interactions_first_state.Model == 1]
        for (nodeId1, interaction, nodeId2, *_) in interactions_first_state.itertuples(index=False):
            chain1, *_ = nodeId1.split(":")
            chain2, *_ = nodeId2.split(":")
            intType, *_ = interaction.split(":")
            keyyyy = tuple(sorted([chain1, chain2]))
            edges.setdefault(keyyyy, dict())
            edges[keyyyy].setdefault(intType, 0)
            edges[keyyyy][intType] += 1
            if edges[keyyyy][intType] == 1:
                G.add_edge(chain1, chain2, type=intType)

        seen = dict()
        present_interaction = set()
        pos = nx.kamada_kawai_layout(G)
        nx.draw_networkx_nodes(G, pos)
        nx.draw_networkx_labels(G, pos)
        ax = plt.gca()
        for e in G.edges(data=True):
            seen.setdefault((e[0], e[1]), 0)
            present_interaction.add(e[2]["type"])
            ax.annotate("",
                        xy=pos[e[0]], xycoords='data',
                        xytext=pos[e[1]], textcoords='data',
                        arrowprops=dict(arrowstyle="<->", color=intColorMap[e[2]["type"]],
                                        shrinkA=5, shrinkB=5,
                                        patchA=None, patchB=None,
                                        connectionstyle="arc3,rad={}".format(0.2 * seen[(e[0], e[1])]),
                                        ),
                        )
            seen[(e[0], e[1])] += 1

        handler_list = []
        for k, v in intColorMap.items():
            if k in present_interaction:
                handler_list.append(mpatches.Patch(color=v, label=k))

        plt.legend(handles=handler_list, loc='best')

        plt.title("Chain interaction graph for {}".format(stored.model))
        plt.axis('off')
        plt.show()

    def inter_heatmap():
        from pymol import stored
        try:
            from matplotlib import pyplot as plt
            import seaborn as sn
        except ImportError:
            log("To run this feature you have to install matplotlib and seaborn and for the python version "
                "that PyMol is using", error=True)
            return

        obj = form.selections_list.currentText()
        if len(obj) == 0:
            log("Please select an object on the box of the view/filter tab to use this feature", error=True)
            return

        if obj[0] == "(" and obj[-1] == ")":
            log("Please select an object on the box of the view/filter tab to use this feature", error=True)
            return

        sele_inter = form.interaction_sele.currentText()

        stored.model = ""
        cmd.iterate(obj, 'stored.model = model')
        if cmd.get_chains(stored.model) == 1:
            log("Only one chain is present in the selected object, cannot use this tool", error=True)
            return

        file_pth = "/tmp/ring/md/" + stored.model + ".gfreq_{}".format(sele_inter)
        if not os.path.exists(file_pth):
            log("Before this you need to run Ring-md on the object first. Select it in the View/Filter tab and press "
                "the Show button", error=True)
            return

        contact_freq = dict()
        order = []
        present = set()
        with open(file_pth, 'r') as f:
            for line in f:
                edge1, _, edge2, perc = line.split('\t')
                edge1 = edge1.replace(":_:", ":").replace(":", "/")
                edge2 = edge2.replace(":_:", ":").replace(":", "/")
                chain1 = edge1.split('/')[0]
                chain2 = edge2.split('/')[0]
                if chain1 != chain2:
                    contact_freq.setdefault((edge1, edge2), perc)
                    if edge1 not in present:
                        present.add(edge1)
                        order.append(edge1)

        if len(order) == 0:
            log("No interaction of this type found", error=True)
            return

        matr = np.zeros((len(order), len(order))) * np.nan
        for i, edge1 in enumerate(order):
            for j, edge2 in enumerate(order):
                try:
                    matr[i, j] = contact_freq[(edge1, edge2)]
                    matr[j, i] = contact_freq[(edge1, edge2)]
                except KeyError:
                    pass
        plt.close()
        ax = plt.subplot()
        sn.heatmap(matr, vmin=0, vmax=1, xticklabels=order, yticklabels=order, cmap='viridis', ax=ax)
        change_chain = dict()
        for i, x in enumerate(order):
            change_chain.setdefault(x.split('/')[0], i)
        ax.hlines(list(change_chain.values())[1:], *ax.get_xlim(), colors=["k"])
        ax.vlines(list(change_chain.values())[1:], *ax.get_ylim(), colors=["k"])
        plt.title("{} interchain interactions".format(sele_inter))
        plt.tight_layout()
        plt.show()

    form.ring_path.setText("{}/.ring/bin/Ring-md".format(environ['HOME']))

    # Execute Ring
    form.visualize_btn.clicked.connect(run)
    form.ring_exec_button.clicked.connect(browse_ring_exe)

    # Update view
    form.radius_value.valueChanged.connect(slider_radius_change)
    form.transp_value.valueChanged.connect(slider_transp_change)

    # Residue based analysis repr
    form.bttn_color_resi_freq.clicked.connect(inter_freq_analysis)
    form.bttn_table_freq.clicked.connect(inter_freq_table)
    dialog_freq.freqTable.itemClicked.connect(sele_selected_resi_freq_table)

    # Interaction based analysis
    form.calc_corr.clicked.connect(correlation_obj)
    dialog_corr.show_corr.clicked.connect(show_corr)

    # Analysis
    form.resi_plot.clicked.connect(resi_plot)
    form.chain_graph.clicked.connect(chain_graph)
    form.show_inter_heatmap.clicked.connect(inter_heatmap)

    # Misc
    form.timer = QtCore.QTimer()
    form.timer.timeout.connect(refresh_sele)
    form.timer.start(1500)

    return dialog
