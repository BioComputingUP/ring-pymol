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
    for j in range(1, frames + 1):
        try:
            df = pd.read_csv('/tmp/ring/{}.cif_cm_m{}_{}'.format(obj, j, int_type), sep=' ', header=None)
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
    indexes = np.tril_indices(z.shape[0], -1)
    for i, j in zip(indexes[0], indexes[1]):
        corr_coeff, p_val = pearsonr(z[i], z[j])
        if p_val < p_thresh and (corr_coeff > coeff_thresh or corr_coeff < -coeff_thresh):
            coeffs_matr[i, j] = corr_coeff

    ticks = ["{}_{}".format(x, y) for (x, y) in contacts_sparse.keys()]
    ticks = np.array(ticks)
    # hm = sn.heatmap(coeffs_matr, xticklabels=ticks, yticklabels=ticks)
    # hm.set_xticklabels(hm.get_xmajorticklabels(), fontsize=8)
    # hm.set_yticklabels(hm.get_ymajorticklabels(), fontsize=8)
    # plt.show()
    return ticks, coeffs_matr


# global reference to avoid garbage collection of our dialog
dialog = None
correlations = dict()


def run_plugin_gui():
    """
    Open our custom dialog
    """
    global dialog

    if dialog is None:
        dialog = make_dialog()
    dialog.show()


def make_dialog():
    from pymol import cmd
    from os import environ

    from pymol.Qt import QtWidgets
    from pymol.Qt import QtCore
    from pymol.Qt.utils import loadUi

    app = QtWidgets.QApplication([])

    # create a new Window
    dialog = QtWidgets.QDialog()

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'plugin.ui')
    form = loadUi(uifile, dialog)

    def run():
        from pymol import stored

        import subprocess

        obj_name = form.pymol_obj.currentText()
        if not obj_name:
            print("Please select a Pymol object first!")
            return
        if not os.path.exists('/tmp/ring'):
            os.mkdir('/tmp/ring')

        file_pth = "/tmp/ring/" + obj_name + ".cif"

        if form.all_states.isChecked():
            state = 0
        else:
            state = -1

        form.setEnabled(False)
        form.button_start.setText("Running...")
        app.processEvents()

        cmd.save(filename=file_pth, selection=obj_name, state=state)

        stored.state = ''
        cmd.iterate_state(state=-1, selection=obj_name, expression='stored.state=state')

        nStates = cmd.count_states(obj_name) if form.all_states.isChecked() else 1

        edge_policy = ""
        if form.best_edge.isChecked():
            edge_policy = "--best_edge"
        if form.multi_edge.isChecked():
            edge_policy = "--multi_edge"
        if form.all_edge.isChecked():
            edge_policy = "--all_edges"

        all_states = ''
        model = []
        if form.all_states.isChecked():
            all_states = "--all_models"
        else:
            model = ["-m", str(int(stored.state))]

        seq_sep = str(form.seq_separation.value())
        try:
            p = subprocess.Popen(
                    [form.ring_path.text(), "-i", file_pth, "--out_dir", "/tmp/ring/", "-g", seq_sep,
                     "--write_all", "--md", "--all_chains", edge_policy, all_states] + model, stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE, universal_newlines=True)
        except FileNotFoundError:
            print("Ring path is not correct!")
            form.setEnabled(True)
            return

        while p.poll() is None:
            line = p.stderr.readline()
            if line != "":
                if "model" in line:
                    nModel = int(line.split("model ")[1].strip())
                    form.button_start.setText("Running on model {} | {:.2%}".format(nModel, nModel / nStates))
            app.processEvents()

        form.button_start.setText("Start Ring")
        form.setEnabled(True)
        app.processEvents()
        visualize(first=True)

    def browse_ring_exe():
        filename = QtWidgets.QFileDialog.getOpenFileNames(dialog, "Select Ring executable")[0][0]

        if filename:
            form.ring_path.setText(filename)

    def refresh_sele():
        form.sele_names.blockSignals(True)
        selections = cmd.get_names('public_selections')
        selections.extend(
                list(filter(lambda x: ("_" in x and x.split('_')[1] not in intTypeMap.keys()) or "_" not in x,
                            cmd.get_names('public_nongroup_objects'))))
        selections = sorted(selections)
        if selections != sorted([form.sele_names.itemText(i) for i in range(form.sele_names.count())]):
            form.sele_names.clear()
            form.sele_names.addItems(selections)
        form.sele_names.blockSignals(False)

    def refresh_obj():
        form.pymol_obj.blockSignals(True)
        selections = sorted(
                list(filter(lambda x: ("_" in x and x.split('_')[1] not in intTypeMap.keys()) or "_" not in x,
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
                with open("/tmp/ring/{}.gfreq_{}".format(obj, inter), 'r') as f:
                    for line in f:
                        edge1, _, edge2, perc = line.split('\t')
                        edge1 = edge1.rsplit(':', 2)[0]
                        edge2 = edge2.rsplit(':', 2)[0]
                        conn_freq[inter].setdefault((edge1, edge2), float(perc))
            except FileNotFoundError:
                print("Please run Ring on the selected object first!")
                return None
        return conn_freq

    def get_freq_combined(obj, bond):
        import math
        conn_freq = dict()
        try:
            with open("/tmp/ring/{}.gfreq_{}".format(obj, bond), 'r') as f:
                for line in f:
                    edge1, _, edge2, perc = line.split('\t')
                    edge1 = edge1.rsplit(':', 2)[0]
                    edge2 = edge2.rsplit(':', 2)[0]
                    conn_freq.setdefault(edge1, [])
                    chain1 = edge1.split(':')[0]
                    chain2 = edge2.split(':')[0]
                    if form.interchain_2.isChecked() and chain1 == chain2:
                        continue
                    if form.intrachain_2.isChecked() and chain1 != chain2:
                        continue
                    conn_freq[edge1].append(float(perc))
        except FileNotFoundError:
            print("Please run Ring on the selected object first!")
            return None
        for k, v in conn_freq.items():
            conn_freq[k] = 1 - math.prod([(1 - x) for x in v])

        return conn_freq

    def visualize(first=False, selection=None, color=None, int_type=None, pair_set=None):
        from pymol import stored
        form.setEnabled(False)
        app.processEvents()
        if selection:
            obj = selection
        else:
            if first:
                obj = form.pymol_obj.currentText()
            else:
                obj = form.sele_names.currentText()

        if obj == '':
            print("Please provide a selection")
            return

        cmd.delete(obj + "_interact")

        states = int(cmd.count_states(selection=obj))

        for state in range(1, states + 1):
            current_state = cmd.get_state()
            if not form.all_states.isChecked() and state != current_state and not selection:
                continue
            stored.model = ""
            cmd.iterate(obj, 'stored.model = model')
            stored.chain_resi = set()

            conn_freq = get_freq(stored.model)

            if conn_freq is None:
                form.setEnabled(True)
                app.processEvents()
                return

            cmd.iterate(obj, 'stored.chain_resi.add((chain, resi))')
            stored.coords = dict()
            cmd.iterate_state(state=state, selection=obj,
                              expression='stored.coords["chain {} and resi {} and name {}".format(chain, resi, name)] = [x,y,z]')

            file_pth = "/tmp/ring/" + stored.model + ".cif_ringEdges" + "_m" + str(state)

            interactions_per_type = dict()

            if not os.path.exists(file_pth):
                print("Please run Ring on the selected object first!")
                return

            with open(file_pth, 'r') as f:
                next(f)
                for line in f:
                    nodeId1, interaction, nodeId2, distance, angle, energy, atom1, \
                    atom2, donor, positive, cation, orientation = line[0:-1].split("\t")[0:12]
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
                                          object_name=obj + "_" + intType if not selection else selection + "_inter",
                                          color=intTypeMap[intType] if not color else color, coords=stored.coords,
                                          state=state)

            if not_present > 0:
                print("{} connections not displayed because atoms not present".format(not_present))

        if not selection:
            if form.check_hide_others.isChecked():
                cmd.hide(selection="*_interact")

            members = ""
            for k in intTypeMap.keys():
                members += " {}_{}".format(obj, k)
            cmd.group(obj + "_interact", members=members)
        else:
            cmd.hide(selection="*_interact")

        form.setEnabled(True)
        app.processEvents()

    def inter_proba_colors():
        obj = form.sele_obj.currentText()
        if obj == '':
            print("Please provide a selection")
            return

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
            express = "b=dict_freq['{}:{}'.format(chain,resi)] if '{}:{}'.format(chain,resi) in dict_freq.keys() else 0.001"
            cmd.alter_state(-1, obj, expression=express, space=myspace)
            cmd.spectrum("b", "white yellow orange red", obj, minimum=0.001, maximum=1.0)

    def slider_radius_change():
        value = form.radius_slider.value() * 0.3
        cmd.set("cgo_line_width", value)

    def slider_transp_change():
        value = form.transp_slider.value() * 0.08
        cmd.set("cgo_transparency", value)

    def correlation_obj():
        import numpy as np
        obj = form.sele_obj_2.currentText()
        if obj == '':
            print("Please select an object first!")
            return

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

        form.contact_corr_list.clear()
        form.contact_anti_list.clear()
        form.contact_corr_list.setEnabled(False)
        form.contact_anti_list.setEnabled(False)
        app.processEvents()
        coeff_thr = float(form.coeff_thr.value())
        p_thr = float(form.p_thr.value())

        max_presence = float(form.max_presence.value()) / 100
        min_presence = float(form.min_presence.value()) / 100

        states = cmd.count_states(obj)
        if states < 2:
            print("Correlation cannot be calculated if the number of states is {}".format(states))
            return

        try:
            selections, matrix = calculate_correlation(obj, states, int_type=inter,
                                                       coeff_thresh=coeff_thr,
                                                       p_thresh=p_thr, max_presence=max_presence,
                                                       min_presence=min_presence)
        except TypeError:
            return

        correlations.setdefault(obj, dict())
        correlations[obj][inter] = (selections, matrix)

        with np.errstate(divide='ignore', invalid='ignore'):
            indexes_corr = np.argwhere(~np.isnan(matrix) * matrix > 0)
            indexes_anti = np.argwhere(~np.isnan(matrix) * matrix < 0)
        indexes_corr = set([x[1] for x in indexes_corr])
        indexes_anti = set([x[1] for x in indexes_anti])

        valid_sele_corr = [selections[x] for x in indexes_corr]
        valid_sele_anti = [selections[x] for x in indexes_anti]

        corr_items = ["ALL"] + ["{} - {}".format(x.split('_')[0], x.split('_')[1]) for x in
                                sorted(valid_sele_corr, key=lambda x: int(x.split('_')[0].split(':')[1]))]
        anti_items = ["ALL"] + ["{} - {}".format(x.split('_')[0], x.split('_')[1]) for x in
                                sorted(valid_sele_anti, key=lambda x: int(x.split('_')[0].split(':')[1]))]
        form.contact_corr_list.addItems(corr_items if len(corr_items) > 1 else ["No correlated connections"])
        form.contact_anti_list.addItems(anti_items if len(anti_items) > 1 else ["No anti-correlated connections"])
        form.contact_corr_list.setEnabled(len(corr_items) > 1)
        form.contact_anti_list.setEnabled(len(anti_items) > 1)
        form.Correlate.setEnabled(len(corr_items) > 1)
        form.anticorrelate.setEnabled(len(anti_items) > 1)

    def correlates(anti=False):
        import numpy as np

        obj = form.sele_obj_2.currentText()
        if obj == '':
            print("Please select an object first!")
            return

        if not anti:
            contact = form.contact_corr_list.currentText()
        else:
            contact = form.contact_anti_list.currentText()

        if contact == '':
            print("Please select a contact")
            return

        cmd.color('green', selection=obj)
        cmd.delete("selected or selected_inter or corr_inter or anti-corr_inter or corr or anti-corr")

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

        if contact != "ALL":
            contact = "{}_{}".format(contact.split(' - ')[0], contact.split(' - ')[1])

        try:
            selections, matrix = correlations[obj][inter]
        except KeyError:
            print("Run ring on all states of the object first!")
            return

        with np.errstate(divide='ignore', invalid='ignore'):
            indexes = np.argwhere((matrix != 0) * ~np.isnan(matrix))

        if contact != "ALL":
            request_index = np.argwhere(selections == contact)

            resi1, resi2 = contact.split('_')
            resi1_c, resi1_n = resi1.split(':')
            resi2_c, resi2_n = resi2.split(':')
            cmd.select("selected", "/{}//{}/{} or /{}//{}/{}".format(obj, resi1_c, resi1_n,
                                                                     obj, resi2_c, resi2_n))
            sele_set = {(resi1, resi2)}

        corr_set = set()
        anticorr_set = set()
        for index in indexes:
            index_corr = None
            if contact == "ALL":
                index_corr = index[0]
            else:
                if index[0] == request_index:
                    index_corr = index[1]
                if index[1] == request_index:
                    index_corr = index[0]
            if index_corr:
                resi1, resi2 = selections[index_corr].split('_')
                resi1_c, resi1_n = resi1.split(':')
                resi2_c, resi2_n = resi2.split(':')
                if matrix[index[0], index[1]] > 0 and not anti:
                    cmd.select("corr", "/{}//{}/{} or /{}//{}/{}".format(obj, resi1_c, resi1_n,
                                                                         obj, resi2_c, resi2_n),
                               merge=1)
                    corr_set.add((resi1, resi2))
                elif matrix[index[0], index[1]] < 0 and anti:
                    cmd.select("anti-corr", "/{}//{}/{} or /{}//{}/{}".format(obj, resi1_c, resi1_n,
                                                                              obj, resi2_c, resi2_n),
                               merge=1)
                    anticorr_set.add((resi1, resi2))

        if contact != "ALL":
            visualize(selection="selected", color="blue", int_type=inter, pair_set=sele_set)
            cmd.color(selection="selected", color="blue")

        if not anti:
            visualize(selection="corr", color="yellow", int_type=inter, pair_set=corr_set)
            cmd.color(selection="corr", color="yellow")
        if anti:
            visualize(selection="anti-corr", color="red", int_type=inter, pair_set=anticorr_set)
            cmd.color(selection="anti-corr", color="red")

    def tab_click_update():
        refresh_obj()
        refresh_sele()

    form.ring_path.setText("{}/.ring/bin/Ring".format(environ['HOME']))

    # Execute Ring
    form.button_start.clicked.connect(run)
    form.ring_exec_button.clicked.connect(browse_ring_exe)

    # Update view
    form.visualize_btn.clicked.connect(visualize)
    form.radius_slider.valueChanged.connect(slider_radius_change)
    form.transp_slider.valueChanged.connect(slider_transp_change)

    # Putty repr
    form.show_putty.clicked.connect(inter_proba_colors)

    form.timer = QtCore.QTimer()
    form.timer.timeout.connect(tab_click_update)
    form.timer.start(1500)

    # Correlation
    form.Correlate.clicked.connect(correlates)
    form.anticorrelate.clicked.connect(lambda: correlates(anti=True))
    form.calc_corr.clicked.connect(correlation_obj)
    form.Correlate.setEnabled(False)
    form.anticorrelate.setEnabled(False)
    form.contact_corr_list.setEnabled(False)
    form.contact_anti_list.setEnabled(False)

    return dialog
