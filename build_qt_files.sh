#!/bin/sh
set -eu
pyside6-rcc resources.qrc -o resources_rc.py
pyside6-uic aurora_mainwindow.ui -o ui_aurora_mainwindow.py
python aurora_pyside6_designer.py
