"""Static assay, reporting, and GUI configuration copied from the original application."""


HBV_EXPECTED_GENOME_LENGTH = 3200

HBV_SYSTEM1_FULL_REGIONS = [(251, 1801)]

HBV_SYSTEM2_FULL_REGIONS = [(1689, HBV_EXPECTED_GENOME_LENGTH), (1, 262)]

HBV_SYSTEM1_UNIQUE_REGIONS = [(263, 1688)]

HBV_SYSTEM2_UNIQUE_REGIONS = [(1802, HBV_EXPECTED_GENOME_LENGTH), (1, 250)]

HBV_SYSTEM_OVERLAP_REGIONS = [(251, 262), (1689, 1801)]

MEDAKA_MODEL_HEADER_KEYS = (
    "model_version_id", "basecall_model", "basecaller_model", "model",
)

HCV_MIN_COVERAGE = 100

HCV_MIN_VARIANT_FREQUENCY = 0.15  # 15 percent reporting threshold; lower-frequency variants are hidden from the main report.

HCV_SCRIPT_VERSION = "Aurora3.1.0.2026_07_16"

HCV_NS5A_PANEL_BY_SUBTYPE = {
    "1A": {
        28: {"label": "M28", "reportable_aas": {"A", "G", "T", "V"}},
        30: {"label": "Q30", "reportable_aas": {"H", "R"}},
        31: {"label": "L31", "reportable_aas": {"M", "V"}},
        58: {"label": "H58", "reportable_aas": {"D"}},
        93: {"label": "Y93", "reportable_aas": {"H", "N"}},
    },
    "1B": {
        31: {"label": "L31", "reportable_aas": {"F", "I", "M", "V"}},
        93: {"label": "Y93", "reportable_aas": {"H", "N"}},
    },
    "3A": {
        30: {"label": "A30", "reportable_aas": {"K"}},
        93: {"label": "Y93", "reportable_aas": {"H"}},
    },
}

HCV_NS5A_POSITIONS = tuple(sorted({
    aa_pos
    for subtype_panel in HCV_NS5A_PANEL_BY_SUBTYPE.values()
    for aa_pos in subtype_panel.keys()
}))

HCV_NS5A_REFERENCE_MOTIF_BY_SUBTYPE = {
    "1A": {28: "M", 30: "Q", 31: "L", 58: "H", 93: "Y"},
    "1B": {28: "M", 30: "R", 31: "L", 58: "P", 93: "Y"},
    "3A": {28: "M", 30: "A", 31: "L", 58: "P", 93: "Y"},
}

HCV_MINOR_RAS_MAX_FREQUENCY = 0.50

HCV_REFERENCE_METADATA = {
    # Built-in fallbacks for common HCV references seen in this project.
    # Local hcv_reference_metadata.tsv and FASTA-header metadata are still supported,
    # but these fallbacks prevent common accession rows from being missed.
    "NC_004102.1": {"subtype": "1A", "ns5a_codon1_start": 6258},
    "NC_004102": {"subtype": "1A", "ns5a_codon1_start": 6258},
    "AJ238799.1": {"subtype": "1B", "ns5a_codon1_start": 6258},
    "AJ238799": {"subtype": "1B", "ns5a_codon1_start": 6258},
    "M58335.1": {"subtype": "1B", "ns5a_codon1_start": 6258},
    "M58335": {"subtype": "1B", "ns5a_codon1_start": 6258},
    "AF011753.1": {"subtype": "1A", "ns5a_codon1_start": 6258},
    "AF011753": {"subtype": "1A", "ns5a_codon1_start": 6258},
    "D90208.1": {"subtype": "1B", "ns5a_codon1_start": 6258},
    "D90208": {"subtype": "1B", "ns5a_codon1_start": 6258},
    "NC_009824.1": {"subtype": "3A", "ns5a_codon1_start": 6258},
    "NC_009824": {"subtype": "3A", "ns5a_codon1_start": 6258},
    "AB030907.1": {"subtype": "1B", "ns5a_codon1_start": 6258},
    "AB030907": {"subtype": "1B", "ns5a_codon1_start": 6258},
    "JQ065709.1": {"subtype": "3B", "ns5a_codon1_start": 6258},
    "JQ065709": {"subtype": "3B", "ns5a_codon1_start": 6258},
    "FJ462435.1": {"subtype": "4B", "ns5a_codon1_start": 6258},
    "FJ462435": {"subtype": "4B", "ns5a_codon1_start": 6258},
    "DQ988074.1": {"subtype": "4A", "ns5a_codon1_start": 6258},
    "DQ988074": {"subtype": "4A", "ns5a_codon1_start": 6258},
    "HQ639944.1": {"subtype": "2A", "ns5a_codon1_start": 6258},
    "HQ639944": {"subtype": "2A", "ns5a_codon1_start": 6258},
}

HCV_HTML_SCRIPT_VERSION = "HCV_NS5A_HTML_REPORT_2026_07_20"

HCV_NS5A_REPRESENTATIVE_MEDICATIONS = (
    "Daclatasvir",
    "Ledipasvir",
    "Elbasvir",
    "Ombitasvir",
    "Velpatasvir",
    "Pibrentasvir",
)

HBV_MUTATION_SCRIPT_VERSION = "Aurora3.1.0.HBV_LOCAL_DEPTH_SCOPE_2026_10_05"

HBV_MIN_COVERAGE = 100
# Findings at 30-99 complete codon-spanning reads may be reported with an explicit low-depth disclaimer.
# This is HBV-specific and does not alter HCV or other pipeline thresholds.
HBV_MIN_REPORTABLE_COVERAGE = 30

HBV_MIN_VARIANT_FREQUENCY = 0.15

HBV_RT_PANEL = {
    # Local priority markers: prominently flagged in the report rather than
    # downgraded to a generic compensatory warning.
    80:  {"label": "rtL80",  "aas": {"I", "V"}, "category": "priority resistance-associated marker", "flag": "PRIORITY", "drugs": "lamivudine/telbivudine-associated pattern"},
    173: {"label": "rtV173", "aas": {"L"}, "category": "priority resistance-associated marker", "flag": "PRIORITY", "drugs": "lamivudine-associated pattern"},

    # Well-established primary or key resistance-pattern substitutions.
    180: {"label": "rtL180", "aas": {"M"}, "category": "key resistance-pattern substitution", "flag": "ESTABLISHED", "drugs": "lamivudine; entecavir resistance background"},
    181: {"label": "rtA181", "aas": {"T", "V"}, "category": "primary resistance-associated", "flag": "ESTABLISHED", "drugs": "adefovir; reduced-susceptibility patterns"},
    204: {"label": "rtM204", "aas": {"I", "S", "V"}, "category": "primary resistance-associated", "flag": "ESTABLISHED", "drugs": "lamivudine/telbivudine"},
    236: {"label": "rtN236", "aas": {"T"}, "category": "primary resistance-associated", "flag": "ESTABLISHED", "drugs": "adefovir"},

    # These substitutions are biologically relevant but should not be interpreted
    # as standalone entecavir resistance without the appropriate background.
    169: {"label": "rtI169", "aas": {"T"}, "category": "context-dependent entecavir-associated", "flag": "CAUTION", "drugs": "interpret with rtL180M + rtM204V/I"},
    184: {"label": "rtT184", "aas": set("ACFGILMS"), "category": "context-dependent entecavir-associated", "flag": "CAUTION", "drugs": "interpret with rtL180M + rtM204V/I"},
    202: {"label": "rtS202", "aas": {"C", "G", "I"}, "category": "context-dependent entecavir-associated", "flag": "CAUTION", "drugs": "interpret with rtL180M + rtM204V/I"},
    250: {"label": "rtM250", "aas": {"I", "L", "V"}, "category": "context-dependent entecavir-associated", "flag": "CAUTION", "drugs": "interpret with rtL180M + rtM204V/I"},
}

HBV_MEDICATIONS = ("Lamivudine", "Adefovir", "Entecavir", "Tenofovir")

HBV_MUTATION_MEDICATION_MAP = {
    "L80I": {"Lamivudine": "associated", "Entecavir": "associated/background"},
    "L80V": {"Lamivudine": "associated", "Entecavir": "associated/background"},
    "I169T": {"Lamivudine": "associated", "Entecavir": "context-dependent"},
    "V173L": {"Lamivudine": "associated", "Entecavir": "associated/background"},
    "L180M": {"Lamivudine": "established pattern", "Entecavir": "resistance background"},
    "A181T": {"Lamivudine": "associated", "Adefovir": "associated"},
    "A181V": {"Lamivudine": "associated", "Adefovir": "associated"},
    "T184A": {"Entecavir": "context-dependent"},
    "T184G": {"Lamivudine": "associated", "Entecavir": "context-dependent"},
    "T184I": {"Lamivudine": "associated", "Entecavir": "context-dependent"},
    "T184C": {"Lamivudine": "associated", "Entecavir": "context-dependent"},
    "T184M": {"Lamivudine": "associated", "Entecavir": "context-dependent"},
    "T184S": {"Lamivudine": "associated", "Entecavir": "context-dependent"},
    "T184L": {"Lamivudine": "associated", "Tenofovir": "candidate association - caution"},
    "S202G": {"Lamivudine": "associated", "Entecavir": "context-dependent"},
    "S202I": {"Lamivudine": "associated", "Entecavir": "context-dependent"},
    "M204I": {"Lamivudine": "established", "Entecavir": "resistance background"},
    "M204V": {"Lamivudine": "established", "Entecavir": "resistance background"},
    "N236T": {"Adefovir": "established", "Tenofovir": "reported association - caution"},
    "M250V": {"Entecavir": "context-dependent"},
    "M250I": {"Entecavir": "context-dependent"},
    "M250L": {"Entecavir": "context-dependent"},
}

HBV_S_ESCAPE_PANEL = {
    120: {"label": "sP120", "aas": {"L", "Q", "S", "T"}, "category": "escape-associated"},
    123: {"label": "sT123", "aas": {"A", "N"}, "category": "escape-associated"},
    124: {"label": "sC124", "aas": {"R", "Y"}, "category": "escape-associated"},
    126: {"label": "sT126", "aas": {"A", "I", "N", "S"}, "category": "escape-associated"},
    127: {"label": "sP127", "aas": {"L", "S"}, "category": "escape-associated"},
    128: {"label": "sI128", "aas": {"V"}, "category": "escape-associated"},
    129: {"label": "sQ129", "aas": {"H", "N", "R"}, "category": "escape-associated"},
    130: {"label": "sG130", "aas": {"N", "R"}, "category": "escape-associated"},
    131: {"label": "sN131", "aas": {"I"}, "category": "escape-associated"},
    133: {"label": "sM133", "aas": {"I", "L", "T"}, "category": "escape-associated"},
    134: {"label": "sF134", "aas": {"H", "L", "S", "Y"}, "category": "escape-associated"},
    141: {"label": "sK141", "aas": {"E", "I"}, "category": "escape-associated"},
    142: {"label": "sP142", "aas": {"L", "S"}, "category": "escape-associated"},
    143: {"label": "sS143", "aas": {"L", "T"}, "category": "escape-associated", "flag": "HIGH_IMPACT"},
    144: {"label": "sD144", "aas": {"A", "E", "G"}, "category": "escape-associated", "flag": "HIGH_IMPACT"},
    145: {"label": "sG145", "aas": {"A", "E", "K", "R"}, "category": "escape-associated", "flag": "HIGH_IMPACT"},

    # Treatment-selected HBsAg changes caused by the overlapping polymerase/S reading frames.
    # These are retained as detectable surface alterations, but are not by themselves
    # classified as immune/vaccine/HBIG/diagnostic escape.
    164: {"label": "sE164", "aas": {"D"}, "category": "therapy-associated surface alteration", "flag": "THERAPY_SURFACE"},
    172: {"label": "sW172", "aas": {"L", "*"}, "category": "therapy-associated surface alteration", "flag": "THERAPY_SURFACE"},
    195: {"label": "sI195", "aas": {"M"}, "category": "therapy-associated surface alteration", "flag": "THERAPY_SURFACE"},
    196: {"label": "sW196", "aas": {"S", "L", "*"}, "category": "therapy-associated surface alteration", "flag": "THERAPY_SURFACE"},
}

# Mutation-specific HBsAg escape associations used by the HTML report.
# Labels are deliberately conservative: a mutation is assigned only to categories for
# which the curated literature supports that specific substitution. 'Terapiassocierad'
# means selected/created through nucleos(t)ide-analogue resistance in the overlapping
# polymerase/S reading frames; it is not a separate antiviral-resistance call.
HBV_S_ESCAPE_ASSOCIATIONS = {
    "sP120L": {
        "associations": ("Immunescape",),
        "interpretation": "Associerad med förändrad antikroppsigenkänning.",
    },
    "sP120S": {
        "associations": ("Immunescape", "Vaccin-escape", "Diagnostikescape"),
        "interpretation": "Rapporterad i vaccin-/immunescape-sammanhang och kan påverka HBsAg-detektion.",
    },
    "sP120T": {
        "associations": ("Immunescape", "Vaccin-escape", "HBIG-escape", "Diagnostikescape"),
        "interpretation": "Väl beskriven escape-substitution som kan påverka antikroppsigenkänning, HBIG-svar och HBsAg-detektion.",
    },

    "sT123N": {
        "associations": ("Immunescape", "HBIG-escape"),
        "interpretation": "Rapporterad under antikropps-/HBIG-selektion.",
    },

    "sT126A": {
        "associations": ("Immunescape", "Vaccin-escape"),
        "interpretation": "Rapporterad vaccin-/immunescape-substitution i HBsAg.",
    },
    "sT126I": {
        "associations": ("Immunescape", "Vaccin-escape"),
        "interpretation": "Rapporterad vaccin-/immunescape-substitution i HBsAg.",
    },
    "sT126N": {
        "associations": ("Immunescape", "Vaccin-escape"),
        "interpretation": "Rapporterad vaccin-/immunescape-substitution i HBsAg.",
    },
    "sT126S": {
        "associations": ("Immunescape", "Vaccin-escape", "HBIG-escape", "Diagnostikescape"),
        "interpretation": "Rapporterad vid vaccin-/immunescape samt i samband med HBIG och nedsatt HBsAg-detektion.",
    },

    "sP127L": {
        "associations": ("Immunescape",),
        "interpretation": "Rapporterad immunescape-associerad substitution i major hydrophilic region.",
    },

    # Position 128 is genotype-/referenceberoende och bör tolkas försiktigt.
    "sI128V": {
        "associations": ("Vaccin-escape",),
        "interpretation": "Rapporterad vaccin-escape-association; tolkning bör beakta genotyp och referenssekvens.",
    },

    "sQ129H": {
        "associations": ("Immunescape", "Vaccin-escape", "HBIG-escape", "Diagnostikescape"),
        "interpretation": "Rapporterad escape-substitution med association till vaccin/HBIG och möjlig påverkan på HBsAg-detektion.",
    },
    "sQ129R": {
        "associations": ("Immunescape", "Vaccin-escape", "Diagnostikescape"),
        "interpretation": "Rapporterad vaccin-/immunescape-substitution med möjlig påverkan på HBsAg-detektion.",
    },

    "sG130N": {
        "associations": ("Immunescape", "Diagnostikescape"),
        "interpretation": "Rapporterad immunescape-substitution med möjlig HBsAg-detekteringssvikt.",
    },

    # Position 131 kan ha olika wild-type-rest mellan genotyper/referenser.
    "sN131I": {
        "associations": ("Vaccin-escape", "Diagnostikescape"),
        "interpretation": "Rapporterad vaccin-/diagnostikescape-association; tolkning bör beakta genotyp och referenssekvens.",
    },

    "sM133I": {
        "associations": ("Immunescape", "Vaccin-escape", "HBIG-escape", "Diagnostikescape"),
        "interpretation": "Rapporterad i vaccin-, immunoglobulin- och diagnostikescape-sammanhang.",
    },
    "sM133L": {
        "associations": ("Immunescape", "Vaccin-escape", "HBIG-escape", "Diagnostikescape"),
        "interpretation": "Rapporterad i vaccin-, immunoglobulin- och diagnostikescape-sammanhang.",
    },
    "sM133T": {
        "associations": ("Immunescape", "Vaccin-escape", "HBIG-escape", "Diagnostikescape"),
        "interpretation": "Rapporterad i immun-, vaccin-, HBIG- och diagnostikescape-sammanhang.",
    },

    "sK141E": {
        "associations": ("Vaccin-escape",),
        "interpretation": "Rapporterad vaccin-escape-substitution.",
    },

    "sP142S": {
        "associations": ("Immunescape", "Vaccin-escape"),
        "interpretation": "Rapporterad vaccin-/immunescape-substitution.",
    },

    "sS143L": {
        "associations": ("Immunescape", "Vaccin-escape", "Diagnostikescape"),
        "interpretation": "Escape-associerad substitution som kan påverka antikroppsigenkänning och HBsAg-detektion.",
    },
    "sS143T": {
        "associations": ("Immunescape", "Vaccin-escape"),
        "interpretation": "Rapporterad vaccin-/immunescape-substitution.",
    },

    "sD144A": {
        "associations": ("Immunescape", "Vaccin-escape", "HBIG-escape", "Diagnostikescape"),
        "interpretation": "Väl beskriven escape-substitution med association till vaccin/HBIG och möjlig påverkan på HBsAg-detektion.",
    },
    "sD144E": {
        "associations": ("Immunescape", "Vaccin-escape", "HBIG-escape"),
        "interpretation": "Rapporterad vaccin-/HBIG-/immunescape-substitution.",
    },

    "sG145A": {
        "associations": ("Immunescape", "Vaccin-escape", "HBIG-escape", "Diagnostikescape"),
        "interpretation": "Högpåverkande escape-substitution i a-determinanten med möjlig påverkan på HBsAg-detektion.",
    },
    "sG145E": {
        "associations": ("Immunescape", "HBIG-escape"),
        "interpretation": "Rapporterad HBIG-/immunescape-substitution.",
    },
    "sG145K": {
        "associations": ("Immunescape", "Vaccin-escape"),
        "interpretation": "Rapporterad vaccin-/immunescape-substitution.",
    },
    "sG145R": {
        "associations": ("Immunescape", "Vaccin-escape", "HBIG-escape", "Diagnostikescape"),
        "interpretation": "Klassisk vaccin-escape-substitution; kan även påverka HBIG-svar och HBsAg-detektion.",
    },

    # Terapiassocierade ytantigenförändringar p.g.a. överlappande polymerase/S-ORF.
    "sE164D": {
        "associations": ("Terapiassocierad", "Immunescape"),
        "interpretation": "Kopplad till rtV173L vid antiviral selektion och associerad med reducerad anti-HBs-igenkänning.",
    },
    "sW172L": {
        "associations": ("Terapiassocierad",),
        "interpretation": "Ytantigenförändring kopplad till rtA181T; bedöms tillsammans med RT-fynd.",
    },
    "sW172*": {
        "associations": ("Terapiassocierad",),
        "interpretation": "Stopkodon kopplat till rtA181T och associerat med störd HBsAg/virion-sekretion.",
    },

    "sI195M": {
        "associations": ("Terapiassocierad", "Immunescape", "Vaccin-escape", "Diagnostikescape"),
        "interpretation": "Kopplad till rtM204V; kan minska anti-HBs-bindning och påverka vaccin-/immunescape samt HBsAg-detektion.",
    },

    "sW196S": {
        "associations": ("Terapiassocierad", "Immunescape"),
        "interpretation": "Kopplad till rtM204I och associerad med förändrad HBsAg-antigenicitet.",
    },
    "sW196L": {
        "associations": ("Terapiassocierad", "Immunescape"),
        "interpretation": "Kopplad till rtM204I och associerad med förändrad HBsAg-antigenicitet.",
    },
    "sW196*": {
        "associations": ("Terapiassocierad",),
        "interpretation": "Stopkodon kopplat till rtM204I; kan ge defekt HBsAg-sekretion.",
    },
}
REPORT_LOGO_FILENAMES = (
    "report_logo.png", "aurora_logo.png", "logo.png", "lab_logo.png", "cvr_logo.png"
)

MEDAKA_MODEL_CHOICES = [
    "Auto-detect from FASTQ header",
    "r1041_e82_400bps_sup_v5.2.0",
    "r1041_e82_400bps_hac_v5.0.0",
    "r1041_e82_400bps_sup_v4.3.0",
    "r1041_e82_400bps_sup_v5.0.0",
    "r1041_e82_400bps_sup_v4.2.0",
    "r1041_e82_400bps_sup_v3.5.1",
    "r1041_e82_400bps_sup_variant_v4.2.0",
    "r1041_e82_400bps_hac_v4.2.0",
    "r1041_e82_400bps_hac_variant_v4.2.0",
    "r1041_e82_400bps_fast_g632",
    "r1041_e82_400bps_fast_variant_g632",
    "r1041_e82_260bps_sup_v4.1.0",
    "r1041_e82_260bps_sup_variant_v4.1.0",
    "r1041_e82_260bps_hac_v4.1.0",
    "r1041_e82_260bps_hac_variant_v4.1.0",
    "r1041_e82_260bps_fast_g632",
    "r1041_e82_260bps_fast_variant_g632",
    "r941_min_fast_g507",
    "r941_min_fast_variant_g507",
    "r941_min_hac_g507",
    "r941_min_hac_variant_g507",
    "r941_min_sup_g507",
    "r941_min_sup_variant_g507",
    "r941_prom_fast_g507",
    "r941_prom_fast_variant_g507",
    "r941_prom_hac_g507",
    "r941_prom_hac_variant_g507",
    "r941_prom_sup_g507",
    "r941_prom_sup_variant_g507",
]
