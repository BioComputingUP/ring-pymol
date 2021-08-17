from pymol import cmd

from pymol.Qt import QtWidgets
from pymol.Qt import QtCore
from pymol.Qt.utils import loadUi

import os

from utilities import get_bg_fg_colors, get_freq_combined, intTypeMap


class FreqDialog(QtWidgets.QDialog):
    def __init__(self, main_dialog, parent=None):
        super(FreqDialog, self).__init__(parent)

        self.setWindowFlags(self.windowFlags() & QtCore.Qt.WindowMinimizeButtonHint)

        # populate the Window from our *.ui file which was created with the Qt Designer
        uifile = os.path.join(os.path.dirname(__file__), 'GUIs/frequency.ui')
        loadUi(uifile, self)
        self.parent = main_dialog

        self.freqTable.itemClicked.connect(self.sele_selected_resi_freq_table)

    def inter_freq_table(self):
        try:
            import numpy as np
            import pandas as pd
        except ImportError:
            self.parent.log("Please install numpy and pandas to use this plugin")
            return

        obj = self.parent.widg.selections_list.currentText()
        if obj == '' or (obj[0] == "(" and obj[-1] == ")"):
            self.parent.log("Please select an object to use this feature", error=True)
            return
        self.parent.disable_window()
        self.parent.log("Creation of residue interaction frequency table")
        freq_bond = dict()
        for bondType in intTypeMap.keys():
            try:
                freq_bond.setdefault(bondType,
                                     get_freq_combined(obj, bondType,
                                                       interchain=self.parent.widg.interchain.isChecked(),
                                                       intrachain=self.parent.widg.intrachain.isChecked()))
            except FileNotFoundError:
                self.parent.log("Run Ring on the selected object first!", error=True)
                self.parent.enable_window()
                return

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
        tableWidget = self.freqTable
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
                    bg_color, fg_color, color = get_bg_fg_colors(color)

                if i == 0:
                    tableWidget.setItem(rowPosition, i, QtWidgets.QTableWidgetItem(item))
                else:
                    wItem = QtWidgets.QTableWidgetItem()
                    wItem.setData(QtCore.Qt.DisplayRole, round(float(item), 2))
                    tableWidget.setItem(rowPosition, i, wItem)
                tableWidget.item(rowPosition, i).setBackground(bg_color)
                tableWidget.item(rowPosition, i).setForeground(fg_color)
                tableWidget.item(rowPosition, i).setTextAlignment(QtCore.Qt.AlignCenter)
                if i == 0:
                    prevEdge = item
        tableWidget.setSortingEnabled(True)
        tableWidget.viewport().update()
        self.show()
        self.parent.enable_window()

    def sele_selected_resi_freq_table(self):
        tableWidget = self.freqTable
        indexes = tableWidget.selectionModel().selectedRows()
        for index in sorted(indexes):
            chain, resi = tableWidget.item(index.row(), 0).text().split(":")[0:2]
            cmd.select("sele_row", selection="chain {} and resi {}".format(chain, resi), merge=1)
            self.parent.log("Updated selection sele_row with the residue selected in the frequency table")
