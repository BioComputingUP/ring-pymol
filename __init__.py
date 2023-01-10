# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.
import os
import sys
import time
from pathlib import Path

from pymol import cmd

BASE_DIR = Path(__file__).resolve().parent


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
    from pymol import Qt
    from qt_material import apply_stylesheet

    from main_window import MainDialog

    global dialog

    app = Qt.QtWidgets.QApplication([])

    if "Fusion" in Qt.QtWidgets.QStyleFactory.keys():
        app.setStyle('Fusion')

    extra = {
        # Density Scale
        'density_scale': '-1',

        'font_family': 'Roboto',

        # environ
        'pyside6': False,
        'linux': True,
    }

    dialog = MainDialog(app=app)

    apply_stylesheet(dialog, theme='light_blue.xml', extra=extra, invert_secondary=False)

    stylesheet = dialog.styleSheet()
    with open(os.path.join(BASE_DIR, "GUIs", "custom.scss")) as file:
        dialog.setStyleSheet(stylesheet + file.read().format(**os.environ))

    dialog.show()
    if test:
        sys.exit(app.exec_())


if __name__ == '__main__':
    import pymol

    pymol.finish_launching()
    cmd.set("defer_builds_mode", 3)
    cmd.fetch("2h9r")
    # cmd.load("/home/alessio/dynamics/trj.cif")
    # cmd.remove("chain D")
    # cmd.remove("chain C")
    # cmd.set_name("trj", "md-sim")
    cmd.dss()
    cmd.show_as("cartoon")
    cmd.util.cbc()
    cmd.do("ring_plugin")
