# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'aurora_mainwindow_fixed.ui'
##
## Created by: Qt User Interface Compiler version 6.11.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QFormLayout,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit,
    QMainWindow, QPlainTextEdit, QProgressBar, QPushButton,
    QSizePolicy, QSpacerItem, QStatusBar, QVBoxLayout,
    QWidget)
from . import resources_rc

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(573, 795)
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.outerLayout = QVBoxLayout(self.centralwidget)
        self.outerLayout.setSpacing(14)
        self.outerLayout.setObjectName(u"outerLayout")
        self.outerLayout.setContentsMargins(20, 20, 20, 20)
        self.bannerLabel = QLabel(self.centralwidget)
        self.bannerLabel.setObjectName(u"bannerLabel")
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.bannerLabel.sizePolicy().hasHeightForWidth())
        self.bannerLabel.setSizePolicy(sizePolicy)
        self.bannerLabel.setMinimumSize(QSize(0, 70))
        self.bannerLabel.setStyleSheet(u"QLabel {\n"
"    background-color: #2eb;\n"
"    color: white;\n"
"    font-size: 24px;\n"
"    font-weight: bold;\n"
"    padding: 2px;\n"
"}")
        self.bannerLabel.setPixmap(QPixmap(u":/Aurora/images:/croppeddeeraurora.png"))
        self.bannerLabel.setScaledContents(True)
        self.bannerLabel.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.outerLayout.addWidget(self.bannerLabel)

        self.titleLabel = QLabel(self.centralwidget)
        self.titleLabel.setObjectName(u"titleLabel")
        font = QFont()
        font.setPointSize(24)
        font.setBold(True)
        self.titleLabel.setFont(font)
        self.titleLabel.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.outerLayout.addWidget(self.titleLabel)

        self.subtitleLabel = QLabel(self.centralwidget)
        self.subtitleLabel.setObjectName(u"subtitleLabel")
        self.subtitleLabel.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.outerLayout.addWidget(self.subtitleLabel)

        self.settingsBox = QGroupBox(self.centralwidget)
        self.settingsBox.setObjectName(u"settingsBox")
        self.settingsForm = QFormLayout(self.settingsBox)
        self.settingsForm.setObjectName(u"settingsForm")
        self.virusLabel = QLabel(self.settingsBox)
        self.virusLabel.setObjectName(u"virusLabel")

        self.settingsForm.setWidget(0, QFormLayout.ItemRole.LabelRole, self.virusLabel)

        self.virus_combo = QComboBox(self.settingsBox)
        self.virus_combo.addItem("")
        self.virus_combo.addItem("")
        self.virus_combo.addItem("")
        self.virus_combo.addItem("")
        self.virus_combo.addItem("")
        self.virus_combo.setObjectName(u"virus_combo")

        self.settingsForm.setWidget(0, QFormLayout.ItemRole.FieldRole, self.virus_combo)

        self.modeLabel = QLabel(self.settingsBox)
        self.modeLabel.setObjectName(u"modeLabel")

        self.settingsForm.setWidget(1, QFormLayout.ItemRole.LabelRole, self.modeLabel)

        self.mode_combo = QComboBox(self.settingsBox)
        self.mode_combo.addItem("")
        self.mode_combo.setObjectName(u"mode_combo")
        self.mode_combo.setEnabled(False)

        self.settingsForm.setWidget(1, QFormLayout.ItemRole.FieldRole, self.mode_combo)

        self.inputLabel = QLabel(self.settingsBox)
        self.inputLabel.setObjectName(u"inputLabel")

        self.settingsForm.setWidget(2, QFormLayout.ItemRole.LabelRole, self.inputLabel)

        self.inputFieldWidget = QWidget(self.settingsBox)
        self.inputFieldWidget.setObjectName(u"inputFieldWidget")
        self.inputFieldLayout = QHBoxLayout(self.inputFieldWidget)
        self.inputFieldLayout.setObjectName(u"inputFieldLayout")
        self.inputFieldLayout.setContentsMargins(0, 0, 0, 0)
        self.input_edit = QLineEdit(self.inputFieldWidget)
        self.input_edit.setObjectName(u"input_edit")

        self.inputFieldLayout.addWidget(self.input_edit)

        self.input_browse_button = QPushButton(self.inputFieldWidget)
        self.input_browse_button.setObjectName(u"input_browse_button")

        self.inputFieldLayout.addWidget(self.input_browse_button)


        self.settingsForm.setWidget(2, QFormLayout.ItemRole.FieldRole, self.inputFieldWidget)

        self.outputLabel = QLabel(self.settingsBox)
        self.outputLabel.setObjectName(u"outputLabel")

        self.settingsForm.setWidget(3, QFormLayout.ItemRole.LabelRole, self.outputLabel)

        self.outputFieldWidget = QWidget(self.settingsBox)
        self.outputFieldWidget.setObjectName(u"outputFieldWidget")
        self.outputFieldLayout = QHBoxLayout(self.outputFieldWidget)
        self.outputFieldLayout.setObjectName(u"outputFieldLayout")
        self.outputFieldLayout.setContentsMargins(0, 0, 0, 0)
        self.output_edit = QLineEdit(self.outputFieldWidget)
        self.output_edit.setObjectName(u"output_edit")

        self.outputFieldLayout.addWidget(self.output_edit)

        self.output_browse_button = QPushButton(self.outputFieldWidget)
        self.output_browse_button.setObjectName(u"output_browse_button")

        self.outputFieldLayout.addWidget(self.output_browse_button)


        self.settingsForm.setWidget(3, QFormLayout.ItemRole.FieldRole, self.outputFieldWidget)

        self.referenceLabel = QLabel(self.settingsBox)
        self.referenceLabel.setObjectName(u"referenceLabel")

        self.settingsForm.setWidget(4, QFormLayout.ItemRole.LabelRole, self.referenceLabel)

        self.referenceFieldWidget = QWidget(self.settingsBox)
        self.referenceFieldWidget.setObjectName(u"referenceFieldWidget")
        self.referenceFieldLayout = QHBoxLayout(self.referenceFieldWidget)
        self.referenceFieldLayout.setObjectName(u"referenceFieldLayout")
        self.referenceFieldLayout.setContentsMargins(0, 0, 0, 0)
        self.reference_edit = QLineEdit(self.referenceFieldWidget)
        self.reference_edit.setObjectName(u"reference_edit")

        self.referenceFieldLayout.addWidget(self.reference_edit)

        self.reference_browse_button = QPushButton(self.referenceFieldWidget)
        self.reference_browse_button.setObjectName(u"reference_browse_button")

        self.referenceFieldLayout.addWidget(self.reference_browse_button)


        self.settingsForm.setWidget(4, QFormLayout.ItemRole.FieldRole, self.referenceFieldWidget)

        self.medakaLabel = QLabel(self.settingsBox)
        self.medakaLabel.setObjectName(u"medakaLabel")

        self.settingsForm.setWidget(5, QFormLayout.ItemRole.LabelRole, self.medakaLabel)

        self.medaka_combo = QComboBox(self.settingsBox)
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.addItem("")
        self.medaka_combo.setObjectName(u"medaka_combo")
        self.medaka_combo.setMaxVisibleItems(15)

        self.settingsForm.setWidget(5, QFormLayout.ItemRole.FieldRole, self.medaka_combo)

        self.ncbi_checkbox = QCheckBox(self.settingsBox)
        self.ncbi_checkbox.setObjectName(u"ncbi_checkbox")

        self.settingsForm.setWidget(6, QFormLayout.ItemRole.FieldRole, self.ncbi_checkbox)


        self.outerLayout.addWidget(self.settingsBox)

        self.buttonRow = QHBoxLayout()
        self.buttonRow.setObjectName(u"buttonRow")
        self.check_button = QPushButton(self.centralwidget)
        self.check_button.setObjectName(u"check_button")

        self.buttonRow.addWidget(self.check_button)

        self.clear_button = QPushButton(self.centralwidget)
        self.clear_button.setObjectName(u"clear_button")

        self.buttonRow.addWidget(self.clear_button)

        self.run_button = QPushButton(self.centralwidget)
        self.run_button.setObjectName(u"run_button")

        self.buttonRow.addWidget(self.run_button)


        self.outerLayout.addLayout(self.buttonRow)

        self.buttonSpacer = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.outerLayout.addItem(self.buttonSpacer)

        self.progress = QProgressBar(self.centralwidget)
        self.progress.setObjectName(u"progress")
        self.progress.setStyleSheet(u"QProgressBar {\n"
"    background-color: #10131c;\n"
"    border: 1px solid #394056;\n"
"    border-radius: 10px;\n"
"    text-align: center;\n"
"    color: white;\n"
"    font-weight: bold;\n"
"    min-height: 20px;\n"
"}\n"
"\n"
"QProgressBar::chunk {\n"
"    border-radius: 9px;\n"
"\n"
"    background: qlineargradient(\n"
"        x1: 0, y1: 0,\n"
"        x2: 1, y2: 0,\n"
"        stop: 0.00 #7c3aed,\n"
"        stop: 0.20 #ec4899,\n"
"        stop: 0.42 #22d3ee,\n"
"        stop: 0.65 #34d399,\n"
"        stop: 0.82 #38bdf8,\n"
"        stop: 1.00 #8b5cf6\n"
"    );\n"
"}")
        self.progress.setMinimum(0)
        self.progress.setMaximum(1)
        self.progress.setValue(0)

        self.outerLayout.addWidget(self.progress)

        self.log = QPlainTextEdit(self.centralwidget)
        self.log.setObjectName(u"log")
        font1 = QFont()
        font1.setFamilies([u"Menlo"])
        self.log.setFont(font1)
        self.log.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.log.setReadOnly(True)

        self.outerLayout.addWidget(self.log)

        MainWindow.setCentralWidget(self.centralwidget)
        self.statusBar = QStatusBar(MainWindow)
        self.statusBar.setObjectName(u"statusBar")
        MainWindow.setStatusBar(self.statusBar)

        self.retranslateUi(MainWindow)

        self.run_button.setDefault(True)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"Aurora - Virus Analysis", None))
        self.bannerLabel.setText("")
        self.titleLabel.setText(QCoreApplication.translate("MainWindow", u"Aurora", None))
        self.subtitleLabel.setText(QCoreApplication.translate("MainWindow", u"Universal Bioinformatic Tool For Nanopore Analysis", None))
        self.settingsBox.setTitle(QCoreApplication.translate("MainWindow", u"Pipeline settings", None))
        self.virusLabel.setText(QCoreApplication.translate("MainWindow", u"Virus type:", None))
        self.virus_combo.setItemText(0, QCoreApplication.translate("MainWindow", u"Adenovirus", None))
        self.virus_combo.setItemText(1, QCoreApplication.translate("MainWindow", u"Hepatitis B Virus", None))
        self.virus_combo.setItemText(2, QCoreApplication.translate("MainWindow", u"Hepatitis C Virus", None))
        self.virus_combo.setItemText(3, QCoreApplication.translate("MainWindow", u"Cytomegalovirus (CMV)", None))
        self.virus_combo.setItemText(4, QCoreApplication.translate("MainWindow", u"Other", None))

        self.modeLabel.setText(QCoreApplication.translate("MainWindow", u"Analysis mode:", None))
        self.mode_combo.setItemText(0, QCoreApplication.translate("MainWindow", u"Full Pipeline", None))

        self.inputLabel.setText(QCoreApplication.translate("MainWindow", u"Input folder:", None))
        self.input_browse_button.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.outputLabel.setText(QCoreApplication.translate("MainWindow", u"Output folder:", None))
        self.output_browse_button.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.referenceLabel.setText(QCoreApplication.translate("MainWindow", u"Reference folder:", None))
        self.reference_browse_button.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.medakaLabel.setText(QCoreApplication.translate("MainWindow", u"Medaka model:", None))
        self.medaka_combo.setItemText(0, QCoreApplication.translate("MainWindow", u"r1041_e82_400bps_sup_v5.2.0", None))
        self.medaka_combo.setItemText(1, QCoreApplication.translate("MainWindow", u"r1041_e82_400bps_hac_v5.0.0", None))
        self.medaka_combo.setItemText(2, QCoreApplication.translate("MainWindow", u"r1041_e82_400bps_sup_v4.3.0", None))
        self.medaka_combo.setItemText(3, QCoreApplication.translate("MainWindow", u"r1041_e82_400bps_sup_v5.0.0", None))
        self.medaka_combo.setItemText(4, QCoreApplication.translate("MainWindow", u"r1041_e82_400bps_sup_v4.2.0", None))
        self.medaka_combo.setItemText(5, QCoreApplication.translate("MainWindow", u"r1041_e82_400bps_sup_v3.5.1", None))
        self.medaka_combo.setItemText(6, QCoreApplication.translate("MainWindow", u"r1041_e82_400bps_sup_variant_v4.2.0", None))
        self.medaka_combo.setItemText(7, QCoreApplication.translate("MainWindow", u"r1041_e82_400bps_hac_v4.2.0", None))
        self.medaka_combo.setItemText(8, QCoreApplication.translate("MainWindow", u"r1041_e82_400bps_hac_variant_v4.2.0", None))
        self.medaka_combo.setItemText(9, QCoreApplication.translate("MainWindow", u"r1041_e82_400bps_fast_g632", None))
        self.medaka_combo.setItemText(10, QCoreApplication.translate("MainWindow", u"r1041_e82_400bps_fast_variant_g632", None))
        self.medaka_combo.setItemText(11, QCoreApplication.translate("MainWindow", u"r1041_e82_260bps_sup_v4.1.0", None))
        self.medaka_combo.setItemText(12, QCoreApplication.translate("MainWindow", u"r1041_e82_260bps_sup_variant_v4.1.0", None))
        self.medaka_combo.setItemText(13, QCoreApplication.translate("MainWindow", u"r1041_e82_260bps_hac_v4.1.0", None))
        self.medaka_combo.setItemText(14, QCoreApplication.translate("MainWindow", u"r1041_e82_260bps_hac_variant_v4.1.0", None))
        self.medaka_combo.setItemText(15, QCoreApplication.translate("MainWindow", u"r1041_e82_260bps_fast_g632", None))
        self.medaka_combo.setItemText(16, QCoreApplication.translate("MainWindow", u"r1041_e82_260bps_fast_variant_g632", None))
        self.medaka_combo.setItemText(17, QCoreApplication.translate("MainWindow", u"r941_min_fast_g507", None))
        self.medaka_combo.setItemText(18, QCoreApplication.translate("MainWindow", u"r941_min_fast_variant_g507", None))
        self.medaka_combo.setItemText(19, QCoreApplication.translate("MainWindow", u"r941_min_hac_g507", None))
        self.medaka_combo.setItemText(20, QCoreApplication.translate("MainWindow", u"r941_min_hac_variant_g507", None))
        self.medaka_combo.setItemText(21, QCoreApplication.translate("MainWindow", u"r941_min_sup_g507", None))
        self.medaka_combo.setItemText(22, QCoreApplication.translate("MainWindow", u"r941_min_sup_variant_g507", None))
        self.medaka_combo.setItemText(23, QCoreApplication.translate("MainWindow", u"r941_prom_fast_g507", None))
        self.medaka_combo.setItemText(24, QCoreApplication.translate("MainWindow", u"r941_prom_fast_variant_g507", None))
        self.medaka_combo.setItemText(25, QCoreApplication.translate("MainWindow", u"r941_prom_hac_g507", None))
        self.medaka_combo.setItemText(26, QCoreApplication.translate("MainWindow", u"r941_prom_hac_variant_g507", None))
        self.medaka_combo.setItemText(27, QCoreApplication.translate("MainWindow", u"r941_prom_sup_g507", None))
        self.medaka_combo.setItemText(28, QCoreApplication.translate("MainWindow", u"r941_prom_sup_variant_g507", None))

        self.ncbi_checkbox.setText(QCoreApplication.translate("MainWindow", u"Run NCBI pipeline after consensus (Under Construction)", None))
        self.check_button.setText(QCoreApplication.translate("MainWindow", u"Check dependencies", None))
        self.clear_button.setText(QCoreApplication.translate("MainWindow", u"Clear log", None))
        self.run_button.setText(QCoreApplication.translate("MainWindow", u"Run pipeline", None))
        self.progress.setFormat(QCoreApplication.translate("MainWindow", u"Ready", None))
    # retranslateUi

