import datetime
import subprocess
from os import environ

from networkx import MultiGraph, draw_networkx_labels, draw_networkx_nodes, kamada_kawai_layout
from pandas import read_csv
from pymol import stored
from pymol.Qt import QtCore, QtWidgets
from pymol.Qt.utils import loadUi

from correlation_window import CorrelationDialog
from frequency_window import FreqDialog
from rmsd_clustering import *
from utilities import *


class MainDialog(QtWidgets.QDialog):
    def __init__(self, app=None, parent=None):
        super(MainDialog, self).__init__(parent)

        self.setWindowFlags(QtCore.Qt.WindowCloseButtonHint | QtCore.Qt.WindowMinimizeButtonHint)

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
        self.widg.SS_interaction.clicked.connect(self.secondary_structure_graph)

        # Clustering
        self.clustering_runned_ids = set()
        self.widg.hierarchy_plot.clicked.connect(self.hierarchy_plot_fn)
        self.widg.cluster_plot.clicked.connect(self.cluster_plot_fn)
        self.widg.create_obj.clicked.connect(self.create_cluster_obj)
        self.widg.rmsd_box.toggled.connect(lambda: self.checked(rmsd=True))
        self.widg.cluster_box.toggled.connect(lambda: self.checked(rmsd=False))

        # Interactions Color
        self.widg.b_color_h.clicked.connect(lambda: self.pick_color("HBOND"))
        self.widg.b_color_io.clicked.connect(lambda: self.pick_color("IONIC"))
        self.widg.b_color_ss.clicked.connect(lambda: self.pick_color("SSBOND"))
        self.widg.b_color_pipi.clicked.connect(lambda: self.pick_color("PIPISTACK"))
        self.widg.b_color_pica.clicked.connect(lambda: self.pick_color("PICATION"))
        self.widg.b_color_vdw.clicked.connect(lambda: self.pick_color("VDW"))
        self.widg.b_color_iac.clicked.connect(lambda: self.pick_color("IAC"))

        # Misc
        self.init_colors()
        self.widg.timer = QtCore.QTimer()
        self.widg.timer.timeout.connect(self.refresh_sele)
        self.widg.timer.start(1500)
        self.close_progress()
        self.center_qcombobox()

    # Helper functions
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
            self.widg.console_log.repaint()
            self.processEvents()

    def progress(self, p):
        if not self.widg.progress_bar.isVisible():
            self.widg.progress_bar.setVisible(True)
        self.widg.progress_bar.setValue(p)
        self.widg.progress_bar.setFormat("%.02f %%" % p)
        self.processEvents()

    def close_progress(self):
        self.widg.progress_bar.setVisible(False)
        self.processEvents()

    def disable_window(self):
        self.main.setEnabled(False)
        self.widg.visualize_btn.setEnabled(False)
        self.processEvents()

    def enable_window(self):
        self.main.setEnabled(True)
        self.widg.visualize_btn.setEnabled(True)
        self.processEvents()

    def checked(self, rmsd):
        self.widg.cluster_box.blockSignals(True)
        self.widg.rmsd_box.blockSignals(True)
        if rmsd:
            self.widg.cluster_box.setChecked(False)
            self.widg.rmsd_box.setChecked(True)
        else:
            self.widg.rmsd_box.setChecked(False)
            self.widg.cluster_box.setChecked(True)
        self.widg.cluster_box.blockSignals(False)
        self.widg.rmsd_box.blockSignals(False)

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
        filename = QtWidgets.QFileDialog.getOpenFileNames(self, "Select Ring executable")

        if len(filename[0]) > 0:
            self.widg.ring_path.setText(filename[0][0])

    def center_qcombobox(self):
        for item in [self.widg.interaction_sele, self.widg.clustering_method]:
            item.setEditable(True)

            # getting the line edit of combo box
            line_edit = item.lineEdit()

            # setting line edit alignment to the center
            line_edit.setAlignment(QtCore.Qt.AlignCenter)

            # setting line edit to read only
            line_edit.setReadOnly(True)

    # Ring related functions
    def run(self):
        obj_name = self.widg.selections_list.currentText()

        if len(obj_name) == 0:
            self.log("No object selected", error=True)
            return

        # If it is a selection then try to just visualize
        if is_selection(obj_name):
            self.visualize()
            return

        if not os.path.exists('/tmp/ring'):
            os.mkdir('/tmp/ring')

        stored.chains = ""
        cmd.iterate(obj_name, 'stored.chains += chain')
        if stored.chains == "":
            self.log("Pymol object does not contain chain name, please set it before running Ring. (Use alter)",
                     error=True)
            return

        stored.state = ''
        cmd.iterate_state(state=-1, selection=obj_name, expression='stored.state=state')

        current_run_config = self.get_current_run_config()

        # If Ring has already been run on that obj and the run config didn't changed then just visualize results
        if obj_name in self.prev_launch_config.keys() and self.prev_launch_config[obj_name] == current_run_config \
                and not self.widg.override_memory.isChecked():
            self.visualize()
            return

        if len(self.widg.ring_path.text()) == 0:
            self.log("Ring path is not correct! Set it in the configuration tab", error=True)
            return

        self.disable_window()

        file_pth = "/tmp/ring/" + obj_name + ".cif"

        self.log("Exporting pymol object {} in cif format ({})".format(obj_name, file_pth))

        cmd.save(filename=file_pth, selection=obj_name, state=0)
        self.log("Exporting done")

        if not os.path.exists(self.widg.ring_path.text()):
            self.log("Ring path is not correct! Set it in the configuration tab", error=True)
            self.widg.visualize_btn.setText("Show")
            self.enable_window()
            return

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
        self.log("Ring generation started")

        n_states = cmd.count_states(obj_name)
        prev_state = 0
        while p.poll() is None:
            line = p.stderr.readline()
            if line != "":
                if "model" in line:
                    current_state = int(line.split("model ")[1].strip())
                    if current_state > prev_state:
                        self.progress((current_state / n_states) * 100)
                        prev_state = current_state

        choice = QtWidgets.QMessageBox.question(self, 'Export',
                                                "Do you want to export the contact network in cytoscape format?\n"
                                                "It will be saved as:\n{}/{}.json".format(os.getcwd(), obj_name),
                                                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

        if choice == QtWidgets.QMessageBox.Yes:
            export_network_graph(obj_name)
            self.log("Cytoscape network format saved as {}/{}.json".format(os.getcwd(), obj_name))

        self.close_progress()
        self.log("Ring generation finished")
        self.enable_window()
        self.visualize()

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

    def visualize(self, selection=None, color=None, int_type=None):
        from pymol import stored

        if selection:
            obj = selection
        else:
            obj = self.widg.selections_list.currentText()

        if obj == '':
            self.log("Please provide a selection", error=True)
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

        self.log("Drawing started")

        stored.chain_resi = set()
        cmd.iterate(obj, 'stored.chain_resi.add((chain, int(resi)))')
        conn_freq = get_freq(stored.model)

        def draw():
            interactions_per_state = pd.read_csv(file_pth, sep='\t')

            if len(interactions_per_state) == 0:
                self.log("No interactions found in the object", warning=True)
                return

            for state in range(1, states + 1):
                df = interactions_per_state[interactions_per_state.Model == state]

                stored.coords = dict()
                stored.tmp = ""
                cmd.iterate_state(state=state, selection=obj,
                                  expression='stored.tmp = stored.coords.setdefault("{}/{}/{}"'
                                             '.format(chain, resi, name), [x,y,z])')

                interactions_per_type = dict()

                for (nodeId1, interaction, nodeId2, _, _, _, atom1, atom2, *_) in df.itertuples(index=False):
                    intType = interaction.split(":")[0]
                    node1 = Node(nodeId1)
                    node2 = Node(nodeId2)
                    edge = Edge(node1, node2)

                    # Global filters
                    if node1.id_tuple() not in stored.chain_resi or node2.id_tuple() not in stored.chain_resi:
                        continue

                    # Apply filters if selection is not present
                    if not selection:
                        try:
                            freq = conn_freq[intType][edge] * 100
                        except KeyError:
                            freq = 0.5 * 100
                        # Interchain
                        if self.widg.interchain.isChecked() and node1.chain == node2.chain:
                            continue
                        # Intrachain
                        if self.widg.intrachain.isChecked() and node1.chain != node2.chain:
                            continue
                        # Frequency
                        if not self.widg.min_freq.value() <= freq <= self.widg.max_freq.value():
                            continue
                    else:
                        # Apply filter if selection is present
                        if not (int_type == intType or int_type == "ALL"):
                            continue

                    interactions_per_type.setdefault(intType, [])
                    t = tuple()
                    if "," in atom1:
                        t += (atom1,)
                    else:
                        t += ("{}/{}/{}".format(node1.chain, str(node1.resi), atom1),)
                    if "," in atom2:
                        t += (atom2,)
                    else:
                        t += ("{}/{}/{}".format(node2.chain, str(node2.resi), atom2),)
                    interactions_per_type[intType].append(t)

                for intType, interactions in interactions_per_type.items():
                    draw_links(interactions,
                               object_name=obj + "_" + intType + "_cgo" if not selection else selection + "_cgo",
                               color=intTypeMap[intType],
                               coords=stored.coords,
                               state=state)

            cmd.hide(selection="*_edges")
            if not selection:
                self.create_node_edges_sele(stored.model, is_sele, obj)

            # Set transp and radius after updating the CGOs
            self.slider_radius_change()
            self.slider_transp_change()

        if selection is None:
            cmd.async_(draw)
            self.log("Created group {} for interaction edges".format(obj + "_edges"), timed=False)
            self.log("Created group {} for interaction nodes".format(obj + "_nodes"), timed=False)
        else:
            draw()

    def create_node_edges_sele(self, model_name, is_sele, obj):
        members = ""
        for k in intTypeMap.keys():
            members += " {}_{}_cgo".format(obj, k)
        cmd.group(obj + "_edges", members=members)
        if not is_sele:
            members = ""
            for bond in intTypeMap.keys():
                sele = "{}_{}_resi".format(obj, bond)
                members += " {}".format(sele)

                freqs = get_freq_combined(model_name, bond, interchain=self.widg.interchain.isChecked(),
                                          intrachain=self.widg.intrachain.isChecked())

                for node, freq in freqs.items():
                    if self.widg.min_freq.value() <= freq * 100 <= self.widg.max_freq.value():
                        cmd.select(sele,
                                   selection="chain {} and resi {}".format(node.chain, node.resi),
                                   merge=1)
            cmd.group(obj + "_nodes", members=members)

    # Ring Analysis functions
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

        conn_freq = get_freq_combined(obj, inter, interchain=self.widg.interchain.isChecked(),
                                      intrachain=self.widg.intrachain.isChecked(), key_string=True)

        if conn_freq is not None:
            myspace = {'dict_freq': conn_freq}
            express = "b=dict_freq['{}/{}/{}'.format(chain,resi,resn)] if '{}/{}/{}'.format(chain,resi,resn) " \
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
            for i, inter in enumerate(list(intTypeMap.keys()) + ["ALL"]):
                self.progress(i / (len(intTypeMap) + 1) * 100)
                edges, corr_matr, p_matr = calculate_correlation(obj, states, int_type=inter,
                                                                 coeff_thresh=coeff_thr,
                                                                 p_thresh=p_thr, max_presence=max_presence,
                                                                 min_presence=min_presence)
                self.correlations.setdefault(obj, dict())
                self.correlations[obj][inter] = (edges, corr_matr, p_matr)
            self.close_progress()
        except FileNotFoundError:
            self.log("Run ring on all the states first!", error=True)
            self.enable_window()
            return

        self.corr_dialog.create_table(obj)
        self.enable_window()

    def resi_plot_fn(self):
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

        node1, node2 = list(stored.chain_resi)
        node1 = Node(node1[0], node1[1], node1[2])
        node2 = Node(node2[0], node2[1], node2[2])
        edge1 = Edge(sorted([node1, node2]))
        interaction_distance = dict()

        interactions_per_state = read_csv(file_pth, sep='\t')
        for state in range(1, states + 1):
            df = interactions_per_state[interactions_per_state.Model == state]
            for (nodeId1, interaction, nodeId2, distance, _, _, atom1, atom2, *_) in df.itertuples(index=False):
                intType, intSubType = interaction.split(":")
                node1_1 = Node(nodeId1)
                node2_2 = Node(nodeId2)
                edge2 = Edge(sorted([node1_1, node2_2]))

                if edge1 == edge2:
                    interaction_distance.setdefault(intType, np.ones(states) * 999)
                    interaction_distance[intType][state - 1] = min(interaction_distance[intType][state - 1],
                                                                   float(distance))

        for inter in interaction_distance.keys():
            interaction_distance[inter] = list(map(lambda x: x if x != 999 else np.nan, interaction_distance[inter]))

        plt.figure()
        plt.style.use('default')
        something = False
        for inter in interaction_distance.keys():
            something = True
            plt.scatter(np.arange(1, states + 1), interaction_distance[inter], label=inter, color=intTypeMap[inter],
                        marker='.', s=100)
            plt.plot(np.arange(1, states + 1), interaction_distance[inter], color=intTypeMap[inter])

        if something:
            plt.title("{} - {}".format(node1, node2))
            plt.grid()
            plt.ylim(bottom=0, top=max(filter(lambda x: not np.isnan(x),
                                              [item for sublist in list(interaction_distance.values()) for item in
                                               sublist])) + 1)
            plt.xlim(left=0, right=states + 1)
            plt.xlabel("State")
            plt.ylabel("Distance (Å)")
            plt.legend()
            plt.tight_layout()
            plt.show(block=False)
        else:
            self.log("No interactions found between the two selected residues", warning=True)

    def inter_heatmap(self):
        try:
            obj, model, _ = self.get_values_from_input()
        except ValueError:
            return
        sele_inter = self.widg.interaction_sele.currentText()

        stored.model = ""
        cmd.iterate(obj, 'stored.model = model')
        if cmd.get_chains(stored.model) == 1:
            self.log("Only one chain is present in the selected object, cannot use this tool", error=True)
            return

        if sele_inter != "ALL":
            file_pth = "/tmp/ring/md/" + stored.model + ".gfreq_{}".format(sele_inter)
            if not os.path.exists(file_pth):
                self.log("Before this you need to run Ring-md on the object first!", error=True)
                return

            contact_freq = dict()
            order = []
            present = set()
            with open(file_pth, 'r') as f:
                for line in f:
                    node1, _, node2, perc = line.split('\t')
                    node1 = Node(node1)
                    node2 = Node(node2)
                    if node1.chain != node2.chain:
                        contact_freq.setdefault(Edge(node1, node2), perc)
                        if node1 not in present:
                            present.add(node1)
                            order.append(node1)

            if len(order) == 0:
                self.log("No interaction of this type found", error=True)
                return

        else:
            order = get_node_names_ordered(obj)
            contact_freq = get_freq_combined_all_interactions(obj, interchain=True)

            tmp = [x.node1 for x in contact_freq.keys()]
            to_remove = []
            for node in order:
                if node not in tmp:
                    to_remove.append(node)
            for item in to_remove:
                order.remove(item)

        matr = np.zeros((len(order), len(order))) * np.nan
        for i, node1 in enumerate(order):
            for j, node2 in enumerate(order):
                try:
                    matr[i, j] = contact_freq[Edge(node1, node2)]
                    matr[j, i] = contact_freq[Edge(node1, node2)]
                except KeyError:
                    pass
        plt.figure()
        plt.style.use('default')
        ax = plt.subplot()
        str_order = [str(x) for x in order]
        sn.set(font_scale=0.6)
        sn.heatmap(matr, square=False, vmin=0, vmax=1, xticklabels=str_order, yticklabels=str_order, cmap='OrRd',
                   ax=ax)
        change_chain = dict()
        for i, x in enumerate(order):
            change_chain.setdefault(x.chain, i)
        ax.hlines(list(change_chain.values())[1:], *ax.get_xlim(), colors=["k"])
        ax.vlines(list(change_chain.values())[1:], *ax.get_ylim(), colors=["k"])
        plt.title("{} interchain interactions".format(sele_inter))
        plt.tight_layout()
        plt.show(block=False)

    # Chain and SS interaction graphs
    def get_values_from_input(self):
        obj = self.widg.selections_list.currentText()
        if len(obj) == 0 or obj[0] == "(" and obj[-1] == ")":
            self.log("Please select an object to use this feature", error=True)
            raise ValueError

        stored.model = ""
        cmd.iterate(obj, 'stored.model = model')
        file_pth = "/tmp/ring/" + stored.model + ".cif_ringEdges"
        if not os.path.exists(file_pth):
            self.log(
                    "Before this you need to run Ring-md on the object first. Select it above and press the Show button",
                    error=True)
            raise ValueError

        return obj, stored.model, file_pth

    def chain_graph_fn(self):
        try:
            obj, model, file_pth = self.get_values_from_input()
        except ValueError:
            return

        if cmd.get_chains(model) == 1:
            self.log("Only one chain is present in the selected object", warning=True)

        G = MultiGraph()
        edges = dict()
        df = read_csv(file_pth, sep='\t')
        for (nodeId1, interaction, nodeId2, *_) in df.itertuples(index=False):
            chain1, *_ = nodeId1.split(":")
            chain2, *_ = nodeId2.split(":")
            intType, *_ = interaction.split(":")
            key = tuple(sorted([chain1, chain2]))
            edges.setdefault(key, set())
            edges[key].add(intType)

        for k, interactions in edges.items():
            for intType in interactions:
                G.add_node(k[0], chain=k[0])
                G.add_node(k[1], chain=k[1])
                G.add_edge(k[0], k[1], type=intType, freq=1.0)

        self.draw_multigraph(G, "Chain interaction graph for {}".format(model))

    def secondary_structure_graph(self):
        try:
            obj, model, file_pth = self.get_values_from_input()
        except ValueError:
            return

        stored.sec_struct = dict()
        ss_id = dict()
        stored.order = []
        cmd.iterate(obj, 'stored.sec_struct.setdefault((chain, int(resi)), set()).add(ss)')
        cmd.iterate(obj, 'stored.order.append((chain, int(resi)))')

        seen = set()
        stored.order = [x for x in stored.order if not (x in seen or seen.add(x))]

        for k in stored.sec_struct.keys():
            if '' in stored.sec_struct[k]:
                stored.sec_struct[k].remove('')

        n_alpha = dict()
        n_beta = dict()
        is_prev_ss = 0
        for i, chain_resi in enumerate(stored.order):
            chain = chain_resi[0]
            ss_char = stored.sec_struct[chain_resi].pop() if len(stored.sec_struct[chain_resi]) > 0 else ''
            is_alpha = ss_char == "H" or ss_char == "G" or ss_char == "I"
            is_beta = ss_char == "B" or ss_char == "E" or ss_char == "S"
            if is_alpha or is_beta:
                if is_prev_ss == 0 or is_prev_ss == 1 and is_beta or is_prev_ss == 2 and is_alpha:
                    if is_alpha:
                        n_alpha.setdefault(chain, 0)
                        n_alpha[chain] += 1
                        is_prev_ss = 1
                    else:
                        n_beta.setdefault(chain, 0)
                        n_beta[chain] += 1
                        is_prev_ss = 2

                ss_id[chain_resi] = "{} α{}".format(chain, n_alpha[chain]) \
                    if is_alpha else "{} β{}".format(chain, n_beta[chain])
            else:
                is_prev_ss = 0

        G = MultiGraph()
        interactions = read_csv(file_pth, sep='\t')
        frequency = dict()
        max_model = 0
        for (nodeId1, interaction, nodeId2, *_, n_model) in interactions.itertuples(index=False):
            node1 = Node(nodeId1)
            node2 = Node(nodeId2)
            intType, *_ = interaction.split(":")
            if n_model > max_model:
                max_model = n_model
                seen_this_model = set()

            if node1.id_tuple() in ss_id and node2.id_tuple() in ss_id:
                composite_key = (ss_id[node1.id_tuple()], ss_id[(node2.chain, node2.resi)], intType)
                if composite_key not in seen_this_model:
                    frequency.setdefault(composite_key, 0.0)
                    frequency[composite_key] += 1.0
                    G.add_node(ss_id[node1.id_tuple()], chain=node1.chain)
                    G.add_node(ss_id[node2.id_tuple()], chain=node2.chain)
                    seen_this_model.add(composite_key)

        for k, v in frequency.items():
            G.add_edge(k[0], k[1], type=k[2], freq=v / max_model)

        nchar = max([len(x) for x in ss_id.values()])
        if nchar == 3:
            size = 14
        elif nchar == 4:
            size = 12
        elif nchar == 5:
            size = 9
        else:
            size = 7
        self.draw_multigraph(G, "Secondary structure interaction graph for {}".format(model), shrink=16,
                             node_size=1000, text_size=size)

    @staticmethod
    def draw_multigraph(G, title, shrink=12, node_size=700, text_size=12):
        seen = dict()
        present_interaction = set()
        pos = kamada_kawai_layout(G, center=(20, 20))

        plt.figure()
        plt.style.use('default')
        ax = plt.gca()

        default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        chain_id = dict()
        cont = -1
        for n, chain in G.nodes(data="chain"):
            if chain not in chain_id:
                cont += 1
            chain_id.setdefault(chain, cont)

        draw_networkx_nodes(G, pos, nodelist=G.nodes,
                            node_color=[default_colors[chain_id[c] % len(default_colors)] for n, c in
                                        G.nodes(data='chain')],
                            node_size=node_size, ax=ax)
        draw_networkx_labels(G, pos, font_color='white', font_size=text_size, ax=ax)
        for e in G.edges(data=True):
            seen.setdefault((e[0], e[1]), 0)
            present_interaction.add(e[2]["type"])
            ax.annotate("",
                        xy=pos[e[0]], xycoords='data',
                        xytext=pos[e[1]], textcoords='data',
                        arrowprops=dict(arrowstyle="<->",
                                        color=intTypeMap[e[2]["type"]],
                                        shrinkA=shrink, shrinkB=shrink,
                                        patchA=None, patchB=None,
                                        lw=discrete_mapping(e[2]["freq"]),
                                        connectionstyle="arc3,rad={}".format(0.2 * seen[(e[0], e[1])]),
                                        ),
                        )
            seen[(e[0], e[1])] += 1
        handler_list = []
        for k, v in intTypeMap.items():
            if k in present_interaction:
                handler_list.append(mpatches.Patch(color=v, label=k))
        plt.legend(handles=handler_list, loc='best')
        plt.title(title)
        plt.tight_layout()
        plt.axis('off')
        plt.show(block=False)

    # Clustering related functions
    def calculate_clustering(self):
        obj = self.widg.selections_list.currentText()

        if len(obj) == 0 or obj[0] == "(" and obj[-1] == ")":
            self.log("Please select an object to use this feature", error=True)
            raise ValueError

        if not os.path.exists('/tmp/ring'):
            os.mkdir('/tmp/ring')

        if self.widg.CA_atoms.isChecked():
            sele = '{} and name CA'.format(obj)
            obj = '{}_ca'.format(obj)
            file_pth = "/tmp/ring/" + obj + ".xyz"
        else:
            sele = obj
            file_pth = "/tmp/ring/" + obj + ".xyz"

        self.log("Distance matrix for structure {} not found, creation started".format(obj), warning=True)
        self.log("Exporting pymol object {} in xyz format ({})".format(sele, file_pth))
        self.disable_window()
        cmd.save(filename=file_pth, selection=sele, state=0, format='xyz')
        self.log("Exporting done")
        compute_rmsd_dist_matrix(self, obj)
        self.enable_window()

    def init_clustering(self):
        obj = self.widg.selections_list.currentText()
        if len(obj) == 0 or obj[0] == "(" and obj[-1] == ")":
            self.log("Please select an object to use this feature", error=True)
            raise ValueError
        if self.widg.CA_atoms.isChecked():
            obj = '{}_ca'.format(obj)
        n_cluster = None
        rmsd_val = None
        if self.widg.cluster_box.isChecked():
            n_cluster = int(self.widg.n_cluster.value())
        else:
            rmsd_val = float(self.widg.rmsd_val.value())
        file_pth = "/tmp/ring/" + obj + ".npy"

        if not os.path.exists(file_pth) or obj not in self.clustering_runned_ids:
            self.calculate_clustering()
            self.clustering_runned_ids.add(obj)
        method = self.widg.clustering_method.currentText()
        return method, n_cluster, obj, rmsd_val

    def hierarchy_plot_fn(self):
        method, n_cluster, obj, rmsd_val = self.init_clustering()
        hierarchy_cut_plot(self, obj, method, rmsd_val, n_cluster)

    def cluster_plot_fn(self):
        method, n_cluster, obj, rmsd_val = self.init_clustering()
        cluster_distribution_heatmap(self, obj, method, rmsd_val, n_cluster)

    def create_cluster_obj(self):
        method, n_cluster, obj, rmsd_val = self.init_clustering()
        cluster_states_obj(self, obj, method, rmsd_val, n_cluster)

    def init_colors(self):
        for i in intTypeMap.keys():
            self.set_inter_colors(i)

    def pick_color(self, type):
        color = QtWidgets.QColorDialog.getColor()
        intTypeMap[type] = (float(color.red()) / 255.0, float(color.green()) / 255.0, float(color.blue()) / 255.0)
        self.set_inter_colors(type)

    def set_inter_colors(self, type):
        color = intTypeMap[type]
        color = color[0] * 255, color[1] * 255, color[2] * 255
        style_sheet = "QLabel{{background: rgb({},{},{}); border: 1.3px solid black; border-radius: 8px; " \
                      "min-height: 16px;min-width: 16px; max-height: 16px; max-width: 16px;}}".format(*color)
        if type == "HBOND":
            self.widg.color_h.setStyleSheet(style_sheet)
        if type == "IONIC":
            self.widg.color_io.setStyleSheet(style_sheet)
        if type == "PIPISTACK":
            self.widg.color_pipi.setStyleSheet(style_sheet)
        if type == "PICATION":
            self.widg.color_pica.setStyleSheet(style_sheet)
        if type == "SSBOND":
            self.widg.color_ss.setStyleSheet(style_sheet)
        if type == "VDW":
            self.widg.color_vdw.setStyleSheet(style_sheet)
        if type == "IAC":
            self.widg.color_iac.setStyleSheet(style_sheet)
