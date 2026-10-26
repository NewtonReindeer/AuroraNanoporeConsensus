"""PySide6 user interface and background worker integration."""
import os
import sys
import threading

from PySide6.QtCore import QObject, QPropertyAnimation, QEasingCurve, QSettings, QThread, Qt, Signal, Slot, QRectF
from PySide6.QtGui import (
    QColor, QFont, QTextCursor, QPainter, QPainterPath,
    QLinearGradient, QPen, QPixmap,
)
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QFileDialog, QFormLayout,
    QFrame, QGraphicsDropShadowEffect, QHBoxLayout, QLabel, QLineEdit,
    QMainWindow, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton,
    QSizePolicy, QVBoxLayout, QWidget,
)

try:
    from . import resources_rc  # noqa: F401
    from .ui_aurora_mainwindow import Ui_MainWindow
except ImportError:
    try:
        import resources_rc  # type: ignore # noqa: F401
        from ui_aurora_mainwindow import Ui_MainWindow  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            'Missing generated Qt files. Run:\n'
            '  pyside6-rcc resources.qrc -o resources_rc.py\n'
            '  pyside6-uic "aurora_mainwindow_fixed(1).ui" -o ui_aurora_mainwindow.py'
        ) from exc

from .config import MEDAKA_MODEL_CHOICES
from .external_tools import check_required_tools, which_or_none
from .filesystem import app_base_dir
from .pipeline import run_full_pipeline
from .cancellation import PipelineCancelled, request_cancellation, reset_cancellation

try:
    import pysam
except ImportError:
    pysam = None


class ConsoleStream(QObject):
    """Small file-like object that forwards print output to the GUI."""

    text_written = Signal(str)

    def write(self, text):
        if text:
            self.text_written.emit(str(text))
        return len(text or "")

    def flush(self):
        pass

    def isatty(self):
        return False

class NativeOutputCapture:
    """Capture output written directly to file descriptors 1 and 2.

    Redirecting only ``sys.stdout`` does not catch output from programs such as
    minimap2, samtools, and Medaka. On POSIX systems this class temporarily
    redirects the process stdout/stderr file descriptors into a pipe drained by
    a small reader thread, then forwards that text to the Qt log widget.
    """

    def __init__(self, console_stream):
        self.console_stream = console_stream
        self._saved_stdout_fd = None
        self._saved_stderr_fd = None
        self._read_fd = None
        self._reader_thread = None
        self._active = False

    def start(self):
        if os.name != "posix" or self._active:
            return

        read_fd, write_fd = os.pipe()
        self._saved_stdout_fd = os.dup(1)
        self._saved_stderr_fd = os.dup(2)
        self._read_fd = read_fd

        os.dup2(write_fd, 1)
        os.dup2(write_fd, 2)
        os.close(write_fd)

        self._active = True
        self._reader_thread = threading.Thread(
            target=self._read_loop,
            name="AuroraNativeOutputReader",
            daemon=True,
        )
        self._reader_thread.start()

    def _read_loop(self):
        try:
            while True:
                chunk = os.read(self._read_fd, 4096)
                if not chunk:
                    break
                self.console_stream.write(chunk.decode("utf-8", errors="replace"))
        except OSError:
            pass
        finally:
            if self._read_fd is not None:
                try:
                    os.close(self._read_fd)
                except OSError:
                    pass
                self._read_fd = None

    def stop(self):
        if not self._active:
            return

        try:
            os.dup2(self._saved_stdout_fd, 1)
            os.dup2(self._saved_stderr_fd, 2)
        finally:
            for saved_fd in (self._saved_stdout_fd, self._saved_stderr_fd):
                if saved_fd is not None:
                    try:
                        os.close(saved_fd)
                    except OSError:
                        pass
            self._saved_stdout_fd = None
            self._saved_stderr_fd = None
            self._active = False

        if self._reader_thread is not None:
            self._reader_thread.join(timeout=2.0)
            self._reader_thread = None

class PipelineWorker(QObject):
    finished = Signal()
    failed = Signal(str)
    cancelled = Signal()

    def __init__(self, parameters, console_stream):
        super().__init__()
        self.parameters = parameters
        self.console_stream = console_stream

    @Slot()
    def run(self):
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        native_capture = NativeOutputCapture(self.console_stream)

        # Python print() output goes directly to the GUI. Native/subprocess
        # output is captured at the file-descriptor level below.
        sys.stdout = self.console_stream
        sys.stderr = self.console_stream
        native_capture.start()
        try:
            # Fail visibly in the GUI instead of letting run_full_pipeline()
            # print an error and then look like a successful completion.
            check_required_tools()
            run_full_pipeline(**self.parameters)
        except PipelineCancelled:
            self.cancelled.emit()
        except Exception:
            import traceback
            traceback.print_exc()
            self.failed.emit("The pipeline stopped because of an unexpected error.")
        else:
            self.finished.emit()
        finally:
            native_capture.stop()
            sys.stdout = old_stdout
            sys.stderr = old_stderr

class AuroraWindow(QMainWindow):
    """Controller for the window designed in Qt Widgets Designer."""

    def __init__(self):
        super().__init__()
        self.worker_thread = None
        self.worker = None

        # Remember the most recently browsed folder across all Browse buttons
        # and across Aurora sessions. This makes selecting sibling input, output,
        # and reference folders much quicker.
        self.settings = QSettings("Aurora", "Aurora")
        self.last_browse_directory = str(
            self.settings.value("last_browse_directory", app_base_dir()) or app_base_dir()
        )

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Respect the compact 573 x 795 geometry from Qt Designer instead of
        # overriding it with a wide desktop layout.
        self.resize(590, 820)
        self.setMinimumSize(540, 700)
        self.setMaximumWidth(680)

        # Short aliases keep the pipeline/controller code readable while using
        # the exact object names already present in the Designer file.
        self.virus_type_combo = self.ui.virus_combo
        self.input_folder_edit = self.ui.input_edit
        self.output_folder_edit = self.ui.output_edit
        self.reference_folder_edit = self.ui.reference_edit
        self.medaka_model_combo = self.ui.medaka_combo
        self.run_ncbi_checkbox = self.ui.ncbi_checkbox
        # Local BLAST now runs automatically; the legacy NCBI checkbox is obsolete.
        self.run_ncbi_checkbox.setChecked(False)
        self.run_ncbi_checkbox.hide()
        self.progress_bar = self.ui.progress
        self.console = self.ui.log
        self.run_button = self.ui.run_button
        self.clear_button = self.ui.clear_button
        self.check_button = self.ui.check_button
        self.cancel_button = QPushButton("Cancel run", self.ui.centralwidget)
        self.cancel_button.setObjectName("cancel_run_button")
        self.cancel_button.setEnabled(False)
        self.cancel_button.setToolTip("Stop the current Aurora run and terminate active analysis tools.")
        # Put Cancel in the same horizontal row as Check/Clear/Run.  Using the
        # parent widget's layout would add it to outerLayout instead, which can
        # place the button below the log and outside the visible window.
        self.ui.buttonRow.addWidget(self.cancel_button)

        # Let the form fields use the available compact width without imposing
        # the large minimum widths that caused the stretched window.
        self.ui.settingsForm.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )
        self.ui.settingsBox.setMinimumWidth(0)
        self.ui.settingsBox.setMaximumWidth(16777215)
        # Do not let the settings group absorb spare vertical space.  Keeping
        # it at its natural height removes the empty band below Medaka model
        # while preserving the form layout itself.
        self.ui.settingsBox.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Maximum,
        )

        for field_widget in (
            self.ui.inputFieldWidget,
            self.ui.outputFieldWidget,
            self.ui.referenceFieldWidget,
        ):
            field_widget.setMinimumWidth(0)
            field_widget.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Preferred,
            )

        for line_edit in (
            self.input_folder_edit,
            self.output_folder_edit,
            self.reference_folder_edit,
        ):
            line_edit.setMinimumWidth(150)
            line_edit.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Fixed,
            )

        for browse_button in (
            self.ui.input_browse_button,
            self.ui.output_browse_button,
            self.ui.reference_browse_button,
        ):
            browse_button.setFixedWidth(86)

        for button in (self.check_button, self.clear_button, self.run_button, self.cancel_button):
            button.setMinimumWidth(0)
            button.setMaximumWidth(180)

        # Keep the header compact. QLabel's default vertical size policy can
        # absorb spare window height, which makes the Aurora title/subtitle look
        # much farther apart than the layout spacing suggests. Make both labels
        # vertically fixed and use a modest global spacing.
        self.ui.outerLayout.setSpacing(6)
        for header_label in (self.ui.titleLabel, self.ui.subtitleLabel):
            header_label.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Fixed,
            )
            header_label.setContentsMargins(0, 0, 0, 0)
        self.ui.titleLabel.setFixedHeight(34)
        self.ui.subtitleLabel.setFixedHeight(20)

        # Keep the lower part compact while still leaving enough room for logs.
        self.console.setMinimumHeight(105)
        self.console.setMaximumHeight(150)
        # Add a little breathing room inside the log so console text does not
        # sit directly against the widget border.
        self.console.setStyleSheet(
            self.console.styleSheet() + "QPlainTextEdit#log { padding: 4px 7px; }"
        )
        self.progress_bar.setMaximumWidth(16777215)
        self.ui.buttonSpacer.changeSize(
            0,
            0,
            QSizePolicy.Policy.Minimum,
            QSizePolicy.Policy.Fixed,
        )
        self.ui.outerLayout.invalidate()

        # Keep the original full-width banner size and add only a soft glow.
        # The artwork remains controlled by the Designer label rather than being
        # redrawn into a smaller rounded card.
        self.ui.bannerLabel.setScaledContents(True)
        self.ui.bannerLabel.setMinimumHeight(150)
        self.ui.bannerLabel.setMaximumHeight(170)
        self.ui.bannerLabel.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self._create_banner_effect()

        self.status_label = QLabel("Ready", self)
        # Give status text the same small inset as the console instead of
        # letting it sit directly against the window edge.
        self.status_label.setContentsMargins(7, 0, 4, 0)
        self.ui.statusBar.addWidget(self.status_label, 1)

        self.console_stream = ConsoleStream(self)
        self.console_stream.text_written.connect(self.append_console_text)

        # Start with an explicit prompt instead of silently selecting a virus.
        virus_prompt = "Choose Virus..."
        virus_prompt_index = self.virus_type_combo.findText(virus_prompt)
        if virus_prompt_index < 0:
            self.virus_type_combo.insertItem(0, virus_prompt)
            virus_prompt_index = 0
        self.virus_type_combo.setCurrentIndex(virus_prompt_index)

        self.ui.input_browse_button.clicked.connect(
            lambda: self.choose_directory(self.input_folder_edit)
        )
        self.ui.output_browse_button.clicked.connect(
            lambda: self.choose_directory(self.output_folder_edit)
        )
        self.ui.reference_browse_button.clicked.connect(
            lambda: self.choose_directory(self.reference_folder_edit)
        )
        self.run_button.clicked.connect(self.start_pipeline)
        self.cancel_button.clicked.connect(self.cancel_pipeline)
        self.clear_button.clicked.connect(self.console.clear)
        self.check_button.clicked.connect(self.check_dependencies)

        # Keep the optional automatic model detection available even if the
        # current Designer combobox only lists concrete Medaka models.
        auto_model_text = "Auto-detect from FASTQ header"
        auto_model_index = self.medaka_model_combo.findText(auto_model_text)
        if auto_model_index < 0:
            self.medaka_model_combo.insertItem(0, auto_model_text)
            auto_model_index = 0
        self.medaka_model_combo.setCurrentIndex(auto_model_index)

        # The progress duration is unknown, so it is switched to busy mode
        # only while the pipeline is running.
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setVisible(False)

        self._create_progress_glow()

        banner_pixmap = self.ui.bannerLabel.pixmap()
        if banner_pixmap is None or banner_pixmap.isNull():
            self.append_console_text(
                "[GUI WARNING] The banner resource did not load. "
                "Rebuild resources_rc.py and verify the :/ resource path.\n"
            )

    def _create_banner_effect(self):
        """Apply a full-size, softly pulsing aurora glow to the banner."""
        # Keep a crisp teal edge around the original Designer-sized image.
        self.ui.bannerLabel.setStyleSheet(
            "QLabel#bannerLabel {"
            " background: #08162c;"
            " border: 1px solid #2dd4bf;"
            " padding: 0;"
            "}"
        )

        self.banner_glow = QGraphicsDropShadowEffect(self.ui.bannerLabel)
        self.banner_glow.setBlurRadius(24)
        self.banner_glow.setOffset(0, 0)
        self.banner_glow.setColor(QColor(34, 211, 238, 150))
        self.ui.bannerLabel.setGraphicsEffect(self.banner_glow)

        self.banner_glow_animation = QPropertyAnimation(
            self.banner_glow,
            b"blurRadius",
            self,
        )
        self.banner_glow_animation.setDuration(2800)
        self.banner_glow_animation.setKeyValueAt(0.0, 18.0)
        self.banner_glow_animation.setKeyValueAt(0.5, 32.0)
        self.banner_glow_animation.setKeyValueAt(1.0, 18.0)
        self.banner_glow_animation.setEasingCurve(QEasingCurve.Type.InOutSine)
        self.banner_glow_animation.setLoopCount(-1)
        self.banner_glow_animation.start()

    @Slot()
    def check_dependencies(self):
        missing = [
            tool
            for tool in ("minimap2", "samtools", "medaka_consensus")
            if which_or_none(tool) is None
        ]
        if missing:
            message = "Missing tools on PATH: " + ", ".join(missing)
            self.status_label.setText("Dependencies missing")
            self.append_console_text(f"[DEPENDENCY ERROR] {message}\n")
            QMessageBox.warning(self, "Missing dependencies", message)
            return

        message = "minimap2, samtools, and medaka_consensus were found on PATH."
        if pysam is None:
            message += "\n\npysam is not installed; HCV minor-variant and IUB features may be unavailable."
        self.status_label.setText("Dependencies checked")
        self.append_console_text(f"[DEPENDENCIES] {message}\n")
        QMessageBox.information(self, "Dependencies", message)

    @Slot()
    def choose_directory(self, target_edit):
        """Choose a folder, starting where the user last browsed.

        If this particular field already contains a valid directory, that path
        is preferred. Otherwise all Browse buttons share the most recently
        selected directory. The value is persisted with QSettings, so Aurora
        remembers it after the application is restarted.
        """
        current = target_edit.text().strip()
        if current and os.path.isdir(current):
            start = current
        elif self.last_browse_directory and os.path.isdir(self.last_browse_directory):
            start = self.last_browse_directory
        else:
            start = app_base_dir()

        selected = QFileDialog.getExistingDirectory(self, "Select folder", start)
        if selected:
            target_edit.setText(selected)
            self.last_browse_directory = selected
            self.settings.setValue("last_browse_directory", selected)

    def _create_progress_glow(self):
        self.progress_glow = QGraphicsDropShadowEffect(self.progress_bar)
        self.progress_glow.setBlurRadius(20)
        self.progress_glow.setOffset(0, 0)
        self.progress_glow.setColor(QColor(34, 211, 238, 190))
        self.progress_bar.setGraphicsEffect(self.progress_glow)

        self.glow_animation = QPropertyAnimation(
            self.progress_glow,
            b"blurRadius",
            self,
        )
        self.glow_animation.setDuration(1500)
        self.glow_animation.setKeyValueAt(0.0, 14.0)
        self.glow_animation.setKeyValueAt(0.5, 32.0)
        self.glow_animation.setKeyValueAt(1.0, 14.0)
        self.glow_animation.setEasingCurve(QEasingCurve.Type.InOutSine)
        self.glow_animation.setLoopCount(-1)

    @Slot(str)
    def append_console_text(self, text):
        self.console.moveCursor(QTextCursor.MoveOperation.End)
        self.console.insertPlainText(text)
        self.console.ensureCursorVisible()

    def _validated_parameters(self):
        input_folder = self.input_folder_edit.text().strip()
        output_folder = self.output_folder_edit.text().strip()
        reference_folder = self.reference_folder_edit.text().strip()

        selected_virus = self.virus_type_combo.currentText().strip()

        missing = []
        if not selected_virus or selected_virus == "Choose Virus...":
            missing.append("Virus type")
        if not input_folder:
            missing.append("Input folder")
        if not output_folder:
            missing.append("Output folder")
        if not reference_folder:
            missing.append("Reference folder")
        if missing:
            raise ValueError("Please select: " + ", ".join(missing))

        if not os.path.isdir(input_folder):
            raise ValueError("The input folder does not exist.")
        if not os.path.isdir(reference_folder):
            raise ValueError("The reference folder does not exist.")

        selected_model = self.medaka_model_combo.currentText()
        auto_detect_models = selected_model == "Auto-detect from FASTQ header"
        fallback_model = (
            "r1041_e82_400bps_sup_v4.2.0"
            if auto_detect_models
            else selected_model
        )

        return {
            "input_folder": input_folder,
            "output_folder": output_folder,
            "reference_folder": reference_folder,
            "medaka_model": fallback_model,
            "virus_type": selected_virus,
            "run_ncbi": False,
            "auto_detect_models": auto_detect_models,
        }

    def _set_input_controls_enabled(self, enabled):
        for widget in (
            self.virus_type_combo,
            self.ui.mode_combo,
            self.input_folder_edit,
            self.output_folder_edit,
            self.reference_folder_edit,
            self.medaka_model_combo,
            self.ui.input_browse_button,
            self.ui.output_browse_button,
            self.ui.reference_browse_button,
            self.check_button,
        ):
            widget.setEnabled(enabled)

    @Slot()
    def start_pipeline(self):
        if self.worker_thread is not None:
            return

        try:
            parameters = self._validated_parameters()
        except ValueError as exc:
            QMessageBox.warning(self, "Missing or invalid settings", str(exc))
            return

        reset_cancellation()
        self.console.clear()
        self._set_input_controls_enabled(False)
        self.run_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.clear_button.setEnabled(False)
        self.status_label.setText("Running pipeline...")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setVisible(True)
        self.glow_animation.start()

        self.worker_thread = QThread(self)
        self.worker = PipelineWorker(parameters, self.console_stream)
        self.worker.moveToThread(self.worker_thread)

        self.worker_thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.pipeline_finished)
        self.worker.failed.connect(self.pipeline_failed)
        self.worker.cancelled.connect(self.pipeline_cancelled)
        self.worker.finished.connect(self.worker_thread.quit)
        self.worker.failed.connect(self.worker_thread.quit)
        self.worker.cancelled.connect(self.worker_thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.worker.failed.connect(self.worker.deleteLater)
        self.worker.cancelled.connect(self.worker.deleteLater)
        self.worker_thread.finished.connect(self._thread_finished)
        self.worker_thread.finished.connect(self.worker_thread.deleteLater)
        self.worker_thread.start()

    @Slot()
    def cancel_pipeline(self):
        if self.worker_thread is None or not self.worker_thread.isRunning():
            return
        self.cancel_button.setEnabled(False)
        self.status_label.setText("Cancelling pipeline...")
        self.append_console_text("\n[GUI] Cancellation requested. Stopping active analysis tools...\n")
        request_cancellation()

    @Slot()
    def pipeline_cancelled(self):
        self.status_label.setText("Pipeline cancelled")
        self.append_console_text("\n[GUI] Pipeline cancelled by user. No further samples will be processed.\n")
        self._stop_busy_indicator()

    @Slot()
    def pipeline_finished(self):
        self.status_label.setText("Pipeline finished")
        self.append_console_text("\n[GUI] Pipeline finished.\n")
        self._stop_busy_indicator()

    @Slot(str)
    def pipeline_failed(self, message):
        self.status_label.setText("Pipeline failed")
        self.append_console_text(f"\n[GUI ERROR] {message}\n")
        self._stop_busy_indicator()
        print(f"[PIPELINE ERROR] {message}", flush=True)

    def _stop_busy_indicator(self):
        self.glow_animation.stop()
        self.progress_bar.setVisible(False)
        self.run_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.clear_button.setEnabled(True)
        self._set_input_controls_enabled(True)
        # Analysis mode is intentionally fixed in the Designer file.
        self.ui.mode_combo.setEnabled(False)

    @Slot()
    def _thread_finished(self):
        self.worker = None
        self.worker_thread = None

    def closeEvent(self, event):
        if self.worker_thread is not None and self.worker_thread.isRunning():
            QMessageBox.information(
                self,
                "Pipeline still running",
                "The pipeline is still running. Close the window after it finishes.",
            )
            event.ignore()
            return
        event.accept()