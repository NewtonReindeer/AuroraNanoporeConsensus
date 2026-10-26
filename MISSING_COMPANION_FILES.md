# Companion files not present in the upload

The supplied source imports generated Qt files, but those files were not included
in the uploaded material. They therefore cannot be reproduced exactly from the
Python source alone:

- `ui_aurora_mainwindow.py`, generated from `aurora_mainwindow_fixed(1).ui`
- `resources_rc.py`, generated from `resources.qrc`
- the `.ui`, `.qrc`, images, logos, and optional metadata TSV files referenced by the application

Copy those existing project files into this package directory. The runtime error
in `gui.py` retains the original generation commands.
