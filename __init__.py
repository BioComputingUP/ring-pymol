# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.
import sys

from pymol import cmd


def __init_plugin__(app=None):
    """
    Add an entry to the PyMOL "Plugin" menu
    """
    from pymol.plugins import addmenuitemqt
    from pathlib import Path

    sys.path.append(str(Path(__file__).absolute().parent))

    addmenuitemqt('Ring plugin', ring_plugin)


# global reference to avoid garbage collection of our dialog
dialog = None


@cmd.extend
def ring_plugin(test=False):
    from pymol.Qt import QtWidgets
    from main_window import MainDialog

    global dialog

    app = QtWidgets.QApplication([])

    if "Fusion" in QtWidgets.QStyleFactory.keys():
        app.setStyle('Fusion')

    dialog = MainDialog(app=app)
    dialog.show()
    if test:
        sys.exit(app.exec_())


@cmd.extend
def chain_label():
    cmd.alter("resi 526-542", "chain='B'")
    cmd.alter("not chain B", "chain='A'")


if __name__ == '__main__':
    import pymol

    pymol.finish_launching()
    cmd.set("defer_builds_mode", 3)
    cmd.fetch("2h9r")
    # cmd.load("/home/alessio/dynamics/trj.cif")
    cmd.dss()
    cmd.show_as("cartoon")
    cmd.util.cbc()
    cmd.do("ring_plugin")
