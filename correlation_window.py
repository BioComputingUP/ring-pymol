from os import path

from PyQt5 import QtCore, QtWidgets
from PyQt5.uic import loadUi

from utilities import *


class CorrelationDialog(QtWidgets.QDialog):
    def __init__(self, main_dialog, parent=None):
        super(CorrelationDialog, self).__init__(parent)

        self.setWindowFlags(self.windowFlags() & QtCore.Qt.WindowMinimizeButtonHint)

        # populate the Window from our *.ui file which was created with the Qt Designer
        uifile = path.join(path.dirname(__file__), 'GUIs/correlated.ui')
        loadUi(uifile, self)
        self.parent = main_dialog

        self.current_obj = None

        self.show_corr.clicked.connect(self.show_corr_fn)
        self.text_filter.textChanged.connect(self.filter_table)

    def create_table(self, obj):
        try:
            import numpy as np
            import pandas as pd
        except ImportError:
            self.parent.log("Please install numpy and pandas to use this plugin")
            return

        df = pd.DataFrame()

        freq_inter = get_freq(obj, self.parent.temp_dir.name)
        freq_combined = get_freq_combined_all_interactions(obj, self.parent.temp_dir.name)

        self.current_obj = obj

        tableWidget = self.corrTable
        for inter in list(intTypeMap.keys()) + ["ALL"]:
            edges, corr_matr, p_matr = self.parent.correlations[obj][inter]
            with np.errstate(divide='ignore', invalid='ignore'):
                indexes = np.argwhere(~np.isnan(corr_matr))

            edge1s = [edges[x] for x in [y[0] for y in indexes]]
            freqs1 = [freq_inter[inter][edge] if inter != 'ALL' else freq_combined[edge] for edge in edge1s]
            inter_labels = [inter for _ in edge1s]
            edge2s = [edges[x] for x in [y[1] for y in indexes]]
            freqs2 = [freq_inter[inter][edge] if inter != 'ALL' else freq_combined[edge] for edge in edge2s]
            corr_vals = [corr_matr[i, j] for (i, j) in indexes]
            p_vals = [p_matr[i, j] for (i, j) in indexes]

            tableWidget.setRowCount(0)
            tableWidget.setSortingEnabled(False)

            rowPosition = tableWidget.rowCount()  # necessary even when there are no rows in the table
            df_of_interacton = pd.DataFrame(
                    [edge1s, freqs1, inter_labels, edge2s, freqs2, corr_vals, p_vals]).transpose()
            df = df.append(df_of_interacton)

            self.parent.log("{} {} interactions correlates/anti-correlates".format(len(df_of_interacton), inter))

        df = df.sort_values([0, 3, 5], ascending=(True, True, False))

        prev_edge = None
        color = 2
        tableWidget.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        for _, row in df.iterrows():
            x, f1, p, y, f2, z, w, *_ = row.to_list()
            if x != prev_edge:
                bg_color, fg_color, color = get_bg_fg_colors(color)

            tableWidget.insertRow(rowPosition)

            tableWidget.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(str(x)))

            f1Item = QtWidgets.QTableWidgetItem()
            f1Item.setData(QtCore.Qt.EditRole, round(float(f1), 3))
            tableWidget.setItem(rowPosition, 1, f1Item)

            tableWidget.setItem(rowPosition, 2, QtWidgets.QTableWidgetItem(p))
            tableWidget.setItem(rowPosition, 3, QtWidgets.QTableWidgetItem(str(y)))

            f2Item = QtWidgets.QTableWidgetItem()
            f2Item.setData(QtCore.Qt.EditRole, round(float(f2), 3))
            tableWidget.setItem(rowPosition, 4, f2Item)

            zItem = QtWidgets.QTableWidgetItem()
            zItem.setData(QtCore.Qt.EditRole, round(float(z), 2))
            tableWidget.setItem(rowPosition, 5, zItem)

            wItem = QtWidgets.QTableWidgetItem()
            wItem.setData(QtCore.Qt.EditRole, float(w))
            tableWidget.setItem(rowPosition, 6, wItem)
            for i in range(7):
                tableWidget.item(rowPosition, i).setBackground(bg_color)
                tableWidget.item(rowPosition, i).setForeground(fg_color)
                tableWidget.item(rowPosition, i).setTextAlignment(QtCore.Qt.AlignCenter)

            if z > 0:
                tableWidget.item(rowPosition, 5).setForeground(QColor(0, 0, 255))
            else:
                tableWidget.item(rowPosition, 5).setForeground(QColor(255, 0, 0))

            prev_edge = x
            rowPosition += 1
        tableWidget.setSortingEnabled(True)
        tableWidget.viewport().update()
        self.show()

    def filter_table(self):
        filters = self.text_filter.toPlainText().strip().split(' ')
        tableWidget = self.corrTable

        if len(filters) > 0:
            to_filter_total = set()

            for f in filters:
                to_filter = set()
                items = tableWidget.findItems(f, QtCore.Qt.MatchContains)
                for item in items:
                    to_filter.add(item.row())
                if len(to_filter_total) == 0:
                    to_filter_total = to_filter
                else:
                    to_filter_total.intersection_update(to_filter)

            for row in range(tableWidget.rowCount()):
                tableWidget.setRowHidden(row, row not in to_filter_total)

    def show_corr_fn(self):
        self.parent.disable_window()

        cmd.delete("edge1 edge2 edge1_cgo edge2_cgo")

        table = self.corrTable
        selection = table.selectionModel().selectedRows()

        if len(selection) == 0:
            self.parent.log("Please select at least one row first!", error=True)
            return

        sele_set = set()
        corr_set = set()
        anti_set = set()
        for i in range(len(selection)):
            row = selection[i].row()
            edge1 = table.item(row, 0).text()
            inter = table.item(row, 2).text()
            edge2 = table.item(row, 3).text()
            corr_val = float(table.item(row, 5).text())

            edge1 = Edge([Node(node) for node in edge1.split(' - ')])
            edge2 = Edge([Node(node) for node in edge2.split(' - ')])

            print(edge1, edge2)

            cmd.select("edge1", "/{}//{}/{} or /{}//{}/{}".format(self.current_obj, edge1.node1.chain, edge1.node1.resi,
                                                                  self.current_obj, edge1.node2.chain,
                                                                  edge1.node2.resi), merge=1)
            sele_set.add(edge1)

            cmd.select("edge2", "/{}//{}/{} or /{}//{}/{}".format(self.current_obj, edge2.node1.chain, edge2.node1.resi,
                                                                  self.current_obj, edge2.node2.chain,
                                                                  edge2.node2.resi), merge=1)
            if corr_val > 0:
                corr_set.add(edge2)
            else:
                anti_set.add(edge2)

        self.parent.visualize(selection="edge1", int_type=inter)
        if len(corr_set) > 0:
            self.parent.visualize(selection="edge2", int_type=inter)
        if len(anti_set) > 0:
            self.parent.visualize(selection="edge2", int_type=inter)
        self.parent.log("Selection edge1 contains all the residues from the edges selected in the first column",
                        timed=False)
        self.parent.log("Selection edge2 contains all the residues from the edges selected in the second column",
                        timed=False)
        self.parent.log(
                "CGO objects edge1_cgo and edge2_cgo are the selected edges from the first and second column respectively",
                timed=False)

        self.parent.enable_window()
