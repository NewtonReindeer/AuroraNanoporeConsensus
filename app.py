"""Aurora application entry point."""
import os
import sys


def main():
    # Import lazily so non-GUI modules remain usable on systems without PySide6.
    from PySide6.QtGui import QIcon
    from PySide6.QtWidgets import QApplication
    from .gui import AuroraWindow
    from .filesystem import app_base_dir

    app = QApplication(sys.argv)
    app.setApplicationName("Aurora")
    app.setApplicationDisplayName("Aurora")

    icon_path = os.path.join(app_base_dir(), "howling.png")

    if os.path.exists(icon_path):
        app.setWindowIcon(QIcon(icon_path))

    window = AuroraWindow()

    if os.path.exists(icon_path):
        window.setWindowIcon(QIcon(icon_path))

    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())