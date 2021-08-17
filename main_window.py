import datetime
import os
from os import environ

import numpy as np
import pandas as pd
from PyQt5.QtGui import QColor
from pymol import cmd
from pymol.Qt import QtCore, QtWidgets
from pymol.Qt.utils import loadUi

from correlation_window import CorrelationDialog
from frequency_window import FreqDialog
from utilities import calculate_correlation, draw_links, get_freq, get_freq_combined, intTypeMap, intColorMap, \
    is_selection


class MainDialog(QtWidgets.QDialog):
    def __init__(self, app=None, parent=None):
        super(MainDialog, self).__init__(parent)

        self.setWindowFlags(self.windowFlags() & QtCore.Qt.WindowMinimizeButtonHint)

        # populate the Window from our *.ui file which was created with the Qt Designer
        uifile = os.path.join(os.path.dirname(__file__), 'GUIs/plugin.ui')
        self.widg = loadUi(uifile, self)
        self.app = app

        self.prev_launch_config = dict()
        self.correlations = dict()

        self.corr_dialog = CorrelationDialog(self)
        self.freq_dialog = FreqDialog(self)

        if os.path.exists("{}/.ring/bin/Ring-md".format(environ['HOME'])):
            self.widg.ring_path.setText("{}/.ring/bin/Ring-md".format(environ['HOME']))
        else:
            self.log("Ring-md path not found in current directory, please set it manually", warning=True)

        # Execute Ring
        self.widg.visualize_btn.clicked.connect(self.run)
        self.widg.ring_exec_button.clicked.connect(self.browse_ring_exe)

        # Update view
        self.widg.radius_value.valueChanged.connect(self.slider_radius_change)
        self.widg.transp_value.valueChanged.connect(self.slider_transp_change)

        # Residue based analysis repr
        self.widg.bttn_color_resi_freq.clicked.connect(self.inter_freq_analysis)
        self.widg.bttn_table_freq.clicked.connect(self.freq_dialog.inter_freq_table)

        # Interaction based analysis
        self.widg.calc_corr.clicked.connect(self.correlation_obj)

        # Analysis
        self.widg.resi_plot.clicked.connect(self.resi_plot_fn)
        self.widg.chain_graph.clicked.connect(self.chain_graph_fn)
        self.widg.show_inter_heatmap.clicked.connect(self.inter_heatmap)

        # Misc
        self.widg.timer = QtCore.QTimer()
        self.widg.timer.timeout.connect(self.refresh_sele)
        self.widg.timer.start(1500)

    def processEvents(self):
        self.app.processEvents()

    def log(self, s: str, timed=True, error=False, warning=False, process=True):
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

        self.widg.console_log.insertItem(0, out_s)
        if error:
            self.widg.console_log.item(0).setForeground(QColor(237, 67, 55))
        if warning:
            self.widg.console_log.item(0).setForeground(QColor(255, 146, 13))

        if process:
            self.processEvents()

    def disable_window(self):
        self.main.setEnabled(False)
        self.processEvents()

    def enable_window(self):
        self.main.setEnabled(True)
        self.processEvents()

    def slider_radius_change(self):
        value = self.widg.radius_value.value()
        cmd.set("cgo_line_width", value)

    def slider_transp_change(self):
        value = self.widg.transp_value.value() / 100
        cmd.set("cgo_transparency", value)

    def get_current_run_config(self):
        edge_policy = ""
        if self.widg.best_edge.isChecked():
            edge_policy = "--best_edge"
        if self.widg.multi_edge.isChecked():
            edge_policy = "--multi_edge"
        if self.widg.all_edge.isChecked():
            edge_policy = "--all_edges"

        seq_sep = str(self.widg.seq_separation.value())

        len_hbond = str(self.widg.len_hbond.value())
        len_pica = str(self.widg.len_pica.value())
        len_pipi = str(self.widg.len_pipi.value())
        len_salt = str(self.widg.len_salt.value())
        len_ss = str(self.widg.len_ss.value())
        len_vdw = str(self.widg.len_vdw.value())

        return {"-g"   : seq_sep,
                "-o"   : len_salt,
                "-s"   : len_ss,
                "-k"   : len_pipi,
                "-a"   : len_pica,
                "-b"   : len_hbond,
                "-w"   : len_vdw,
                "edges": edge_policy}

    def browse_ring_exe(self):
        filename = QtWidgets.QFileDialog.getOpenFileNames(self, "Select Ring executable")[0][0]

        if filename:
            self.widg.ring_path.setText(filename)

    def run(self):
        from pymol import stored
        import subprocess

        obj_name = self.widg.selections_list.currentText()

        if len(obj_name) == 0:
            self.log("No object selected", error=True)
            return

        # If it is a selection then try to just visualize
        if is_selection(obj_name):
            self.visualize(log_iter=True)
            return

        if not os.path.exists('/tmp/ring'):
            os.mkdir('/tmp/ring')

        stored.chains = ""
        cmd.iterate(obj_name, 'stored.chains += chain')
        if stored.chains == "":
            self.log("Pymol object does not contain chain name, please set it before running Ring. (Use alter)",
                     error=True)
            return

        file_pth = "/tmp/ring/" + obj_name + ".cif"

        stored.state = ''
        cmd.iterate_state(state=-1, selection=obj_name, expression='stored.state=state')

        current_run_config = self.get_current_run_config()

        # If Ring has already been run on that obj and the run config didn't changed then just visualize results
        if obj_name in self.prev_launch_config.keys() and self.prev_launch_config[obj_name] == current_run_config \
                and not self.widg.override_memory.isChecked():
            self.visualize(log_iter=True)
            return

        if len(self.widg.ring_path.text()) == 0:
            self.log("Ring path is not correct! Set it in the configuration tab", error=True)
            return

        self.disable_window()
        self.widg.visualize_btn.setText("Running ring...")

        self.log("Exporting pymol object {} in cif format ({})".format(obj_name, file_pth))

        cmd.save(filename=file_pth, selection=obj_name, state=0)
        self.log("Exporting done")

        try:
            p = subprocess.Popen(
                    [self.widg.ring_path.text(), "-i", file_pth, "--out_dir", "/tmp/ring/",
                     "-g", current_run_config["-g"],
                     "-o", current_run_config["-o"],
                     "-s", current_run_config["-s"],
                     "-k", current_run_config["-k"],
                     "-a", current_run_config["-a"],
                     "-b", current_run_config["-b"],
                     "-w", current_run_config["-w"],
                     "--all_chains", current_run_config["edges"], "--all_models"], stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE, universal_newlines=True)

            self.prev_launch_config[obj_name] = current_run_config
        except FileNotFoundError:
            self.log("Ring path is not correct! Set it in the configuration tab", error=True)
            self.widg.visualize_btn.setText("Show")
            self.enable_window()
            return

        self.log("Ring generation started")

        n_states = cmd.count_states(obj_name)
        while p.poll() is None:
            line = p.stderr.readline()
            if line != "":
                if "model" in line:
                    current_state = int(line.split("model ")[1].strip())
                    self.log("Running on model {} | {:.2%}".format(current_state, current_state / n_states))

        self.log("Ring generation finished")
        self.visualize(log_iter=True)

    def refresh_sele(self):
        self.widg.selections_list.blockSignals(True)
        present = set(self.widg.selections_list.itemText(i) for i in range(self.widg.selections_list.count()))
        selections = set(map(lambda x: "(" + x + ")", cmd.get_names('public_selections')))
        selections.update(list(filter(lambda x: x.split('_')[-1][-3:] != 'cgo',
                                      cmd.get_names('public_nongroup_objects'))))
        not_present_anymore = present - selections
        new_selections = selections - present

        for sele in not_present_anymore:
            for sele_in, idx in [(self.widg.selections_list.itemText(i), i) for i in
                                 range(self.widg.selections_list.count())]:
                if sele == sele_in:
                    self.widg.selections_list.removeItem(idx)
                    break

        for sele in new_selections:
            self.widg.selections_list.addItem(sele, self.widg.selections_list.count())

        self.widg.selections_list.blockSignals(False)

        obj = self.widg.selections_list.currentText()
        if len(obj) > 0 and obj[0] == "(" and obj[-1] == ")":
            self.widg.visualize_btn.setText("Show")
        elif obj in self.prev_launch_config.keys() and self.prev_launch_config[obj] == self.get_current_run_config() \
                and not self.widg.override_memory.isChecked():
            self.widg.visualize_btn.setText("Show")
        else:
            self.widg.visualize_btn.setText("Execute Ring")

    def visualize(self, selection=None, color=None, int_type=None, pair_set=None, block=True, log_iter=False):
        from pymol import stored

        if block:
            self.widg.main.setEnabled(False)
            self.processEvents()
        if selection:
            obj = selection
        else:
            obj = self.widg.selections_list.currentText()

        if obj == '':
            self.log("Please provide a selection", error=True)
            self.widg.main.setEnabled(True)
            return

        is_sele = is_selection(obj)

        obj = obj.lstrip('(').rstrip(')')

        cmd.delete(obj + "_edges")
        cmd.delete(obj + "_nodes")

        states = int(cmd.count_states(selection=obj))
        stored.model = ""
        cmd.iterate(obj, 'stored.model = model')
        file_pth = "/tmp/ring/" + stored.model + ".cif_ringEdges"

        if not os.path.exists(file_pth):
            self.log("Please run Ring on the selected object first!", error=True)
            return

        interactions_per_state = pd.read_csv(file_pth, sep='\t')

        for state in range(1, states + 1):
            df = interactions_per_state[interactions_per_state.Model == state]

            stored.chain_resi = set()

            conn_freq = get_freq(stored.model)

            if conn_freq is None:
                self.widg.main.setEnabled(True)
                self.processEvents()
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
                    if self.widg.interchain.isChecked() and chain1 == chain2:
                        continue
                    if self.widg.intrachain.isChecked() and chain1 != chain2:
                        continue
                tmp1 = ("{}/{}".format(chain1, pos1), "{}/{}".format(chain2, pos2))
                tmp2 = ("{}/{}".format(chain2, pos2), "{}/{}".format(chain1, pos1))
                if (chain1, pos1) in stored.chain_resi and (chain2, pos2) in stored.chain_resi \
                        and (self.widg.min_freq.value() <= freq <= self.widg.max_freq.value() or selection) \
                        and ((selection and (int_type == intType or int_type == "ALL") and
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
                                          color=intTypeMap[intType] if not color or int_type == "ALL" else color,
                                          coords=stored.coords,
                                          state=state)

            if log_iter:
                log_s = "Interactions state {}: ".format(state)
                other = ""
                for intType in sorted(interactions_per_type.keys()):
                    other += "{} {}, ".format(intType, len(interactions_per_type[intType]))
                log_s += other if len(other) > 1 else "No interaction for this state, maybe check the filters?"
                self.log(log_s.rstrip(', '), timed=False, process=False, warning=len(other) == 0)

            if not_present > 0:
                self.log("{} connections not displayed because atoms not present".format(not_present), warning=True,
                         process=False)

        cmd.hide(selection="*_edges")
        if not selection:
            self.create_node_edges_sele(stored.model, is_sele, obj)

        # Set transp and radius after updating the CGOs
        self.slider_radius_change()
        self.slider_transp_change()
        if block:
            self.enable_window()

    def create_node_edges_sele(self, model_name, is_sele, obj):
        members = ""
        for k in intTypeMap.keys():
            members += " {}_{}_cgo".format(obj, k)
        cmd.group(obj + "_edges", members=members)
        self.log("Created group {} for interaction edges".format(obj + "_edges"), timed=False)
        if not is_sele:
            members = ""
            for bond in intTypeMap.keys():
                sele = "{}_{}_resi".format(obj, bond)
                members += " {}".format(sele)

                freqs = get_freq_combined(model_name, bond, interchain=self.widg.interchain.isChecked(),
                                          intrachain=self.widg.intrachain.isChecked())

                for edge, freq in freqs.items():
                    if self.widg.min_freq.value() <= freq * 100 <= self.widg.max_freq.value():
                        cmd.select(sele,
                                   selection="chain {} and resi {}".format(edge.split(':')[0], edge.split(':')[1]),
                                   merge=1)
            cmd.group(obj + "_nodes", members=members)
            self.log("Created group {} for interaction nodes".format(obj + "_edges"), timed=False)

    def inter_freq_analysis(self):
        obj = self.widg.selections_list.currentText().lstrip('(').rstrip(')')
        if obj == '':
            self.log("Please provide a selection", error=True)
            return

        self.disable_window()

        inter = ""
        if self.widg.hbond.isChecked():
            inter = "HBOND"
        if self.widg.ionic.isChecked():
            inter = "IONIC"
        if self.widg.pipistack.isChecked():
            inter = "PIPISTACK"
        if self.widg.pication.isChecked():
            inter = "PICATION"
        if self.widg.vdw.isChecked():
            inter = "VDW"
        if self.widg.ssbond.isChecked():
            inter = "SSBOND"
        if self.widg.iac.isChecked():
            inter = "IAC"
        conn_freq = get_freq_combined(obj, inter)

        if conn_freq is not None:
            myspace = {'dict_freq': conn_freq}
            express = "b=dict_freq['{}:{}:{}'.format(chain,resi,resn)] if '{}:{}:{}'.format(chain,resi,resn) " \
                      "in dict_freq.keys() else 0.001"
            cmd.alter_state(-1, obj, expression=express, space=myspace)
            cmd.spectrum("b", "white yellow orange red", obj, minimum=0.001, maximum=1.0)

        self.log("Residues with interaction of type {} colored based on the frequency of contact".format(inter))
        self.enable_window()

    def correlation_obj(self):
        obj = self.widg.selections_list.currentText()

        if is_selection(obj):
            self.log("Correlation can only be computed on the whole object, please change selection!", error=True)
            return

        if obj == '':
            self.log("Please select an object first!", error=True)
            return

        self.disable_window()
        self.log("Calculation of correlated interactions started")

        coeff_thr = float(self.widg.coeff_thr.value())
        p_thr = float(self.widg.p_thr.value())

        max_presence = float(self.widg.max_presence.value()) / 100
        min_presence = float(self.widg.min_presence.value()) / 100

        states = cmd.count_states(obj)
        if states < 2:
            self.log("Correlation cannot be calculated if the number of states is {}".format(states), error=True)
            self.enable_window()
            return
        try:
            for inter in list(intTypeMap.keys()) + ["ALL"]:
                selections, corr_matr, p_matr = calculate_correlation(obj, states, int_type=inter,
                                                                      coeff_thresh=coeff_thr,
                                                                      p_thresh=p_thr, max_presence=max_presence,
                                                                      min_presence=min_presence)
                self.correlations.setdefault(obj, dict())
                self.correlations[obj][inter] = (selections, corr_matr, p_matr)

        except TypeError:
            self.log("Run ring on all the states first!", error=True)
            self.enable_window()
            return

        self.corr_dialog.create_table(obj)
        self.enable_window()

    def resi_plot_fn(self):
        from pymol import stored

        try:
            from matplotlib import pyplot as plt
        except ImportError:
            self.log("To run this feature you have to install matplotlib for the python version that PyMol is using",
                     error=True)
            return

        obj = self.widg.selections_list.currentText()
        if len(obj) == 0 or obj[0] != "(" or obj[-1] != ")":
            self.log("Please select a selection on the box above to use this feature", error=True)
            return

        states = int(cmd.count_states(selection=obj))
        if states == 1:
            self.log("You can use this feature only on multi-state objects", error=True)
            return

        stored.model = ""
        cmd.iterate(obj, 'stored.model = model')
        file_pth = "/tmp/ring/" + stored.model + ".cif_ringEdges"

        if not os.path.exists(file_pth):
            self.log(
                    "Before this you need to run Ring-md on the whole object first. Select it above and press the Show button",
                    error=True)
            return

        stored.chain_resi = set()
        cmd.iterate(obj, 'stored.chain_resi.add((chain, int(resi), resn))')
        if len(stored.chain_resi) != 2:
            self.log("You need to create a selection with exactly two residues to use this feature", error=True)
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
            plt.ylim(bottom=0, top=max(filter(lambda x: not np.isnan(x),
                                              [item for sublist in list(interaction_distance.values()) for item in
                                               sublist])) + 1)
            plt.xlim(left=0, right=states + 1)
            plt.xlabel("State")
            plt.ylabel("Distance (Å)")
            plt.legend()
            plt.tight_layout()
            plt.show()
        else:
            self.log("No interactions found between the two selected residues", warning=True)

    def chain_graph_fn(self):
        from pymol import stored
        try:
            from matplotlib import pyplot as plt
            import matplotlib.patches as mpatches
            import networkx as nx
        except ImportError:
            self.log("To run this feature you have to install matplotlib and networkx and for the python version "
                     "that PyMol is using", error=True)
            return

        obj = self.widg.selections_list.currentText()
        if len(obj) == 0:
            self.log("Please select an object to use this feature", error=True)
            return

        if obj[0] == "(" and obj[-1] == ")":
            self.log("Please select an object to use this feature", error=True)
            return

        plt.close()
        stored.model = ""
        cmd.iterate(obj, 'stored.model = model')
        file_pth = "/tmp/ring/" + stored.model + ".cif_ringEdges"
        if not os.path.exists(file_pth):
            self.log(
                    "Before this you need to run Ring-md on the object first. Select it above and press the Show button",
                    error=True)
            return

        if cmd.get_chains(stored.model) == 1:
            self.log("Only one chain is present in the selected object", warning=True)

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

    def inter_heatmap(self):
        from pymol import stored
        try:
            from matplotlib import pyplot as plt
            import seaborn as sn
        except ImportError:
            self.log("To run this feature you have to install matplotlib and seaborn and for the python version "
                     "that PyMol is using", error=True)
            return

        obj = self.widg.selections_list.currentText()
        if len(obj) == 0:
            self.log("Please select an object to use this feature", error=True)
            return

        if obj[0] == "(" and obj[-1] == ")":
            self.log("Please select an object to use this feature", error=True)
            return

        sele_inter = self.widg.interaction_sele.currentText()

        stored.model = ""
        cmd.iterate(obj, 'stored.model = model')
        if cmd.get_chains(stored.model) == 1:
            self.log("Only one chain is present in the selected object, cannot use this tool", error=True)
            return

        file_pth = "/tmp/ring/md/" + stored.model + ".gfreq_{}".format(sele_inter)
        if not os.path.exists(file_pth):
            self.log(
                    "Before this you need to run Ring-md on the object first. Select it in the View/Filter tab and press "
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
            self.log("No interaction of this type found", error=True)
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
