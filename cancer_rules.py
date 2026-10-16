# SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

"""
Hard-coded rules and regex for classifying sample metadata of Methbase
methylomes according to cancer status; designed for human samples.
"""

import re

cancer_non_words = [
    "adjacent",
    "para-?cancer",
    "non-?cancer",
    "cancer-?free",
    "non-?tumou?r",
    "non-?malignan",
    "non-?neoplas",
    "matched normal",
    r"\bnormal\b",
    "healthy",
    "unaffected",
    "disease-?free",
    "no known",
]

cancer_non_re = re.compile("(" + "|".join(cancer_non_words) + ")", re.I)

cancer_control_words = [
    r"\bcontrols?\b",
    r"\bctrl\b",
    r"\bctr\b",
    "_control",
    "control_",
]

cancer_control_re = re.compile("(" + "|".join(cancer_control_words) + ")", re.I)

cancer_strong_words = [
    "carcinom",
    "cancer",
    "malignan",
    "metasta",
    "sarcoma",
    "leuk[a]?emi",
    "lymphoma",
    "myeloma",
    "melanoma",
    "blastoma",
    "glioma",
    "mesothelioma",
    "seminoma",
    r"\bhcc\b",
    r"\bnsclc\b",
    r"\bpdac\b",
    r"\bescc\b",
    r"\boscc\b",
    r"\bcrc\b",
    r"\baml\b",
    r"\bcafs?\b",
    r"\bcll\b",
    r"\bcml\b",
    r"\bb-?all\b",
]

cancer_strong_re = re.compile("(" + "|".join(cancer_strong_words) + ")", re.I)

cancer_growth_words = [
    r"\badenoma",
    "polyp",
    "hyperplas",
    "dysplas",
    "metaplas",
    "panin",
    "leiomyoma",
    "fibroid"
    r"\bnevus\b"
    r"\bnevi\b"
    "barrett",
    "precancer",
    "pre-?malignan",
    "non-?invasive neoplas",
    "benign tumou?r",
    "benign neoplas",
    "benign growth",
]

cancer_growth_re = re.compile("(" + "|".join(cancer_growth_words) + ")", re.I)

cancer_benign_words = [
    r"\bbenign\b",
]

cancer_benign_re = re.compile("|".join(cancer_benign_words), re.I)

cancer_weak_words = [
    "tumou?r",
    "neoplas",
    "xenograft",
    r"\bpdx\b",
]

cancer_weak_re = re.compile("(" + "|".join(cancer_weak_words) + ")", re.I)

# sample-level fields (describe the sample itself)
sample_level_fields = [
    "cancer_status",
    "sampling_site",
    "sample_site",
    "biopsy_site",
    "tissue",
    "tissue_type",
    "sample_type",
    "source_name",
    "sample_condition",
    "tumor_stage",
    "sample_name",
    "cell_type",
    "cell_line",
    "sample_title",
    "title",
    "description",
    "sample_comment",
    "source",
    "BIOMATERIAL_TYPE",
    "sample_group",
]

# patient/donor-level fields
patient_donor_fields = [
    "disease",
    "diagnosis",
    "patient_diagnosis",
    "health_state",
    "DONOR_HEALTH_STATUS",
    "subject_status",
    "phenotype",
    "condition",
    "disease_type",
    "pathological_type",
    "group",
    "case_control",
    "status_of_low-risk_mds",
]

# tumor fields
tumor_fields = [
    "tumor_stage",
    "tnm_stage",
    "tumor_grading",
    "tumor_pathological_stage",
    "tumor_source_tissue",
    "disease_staging",
    "developmental/tumor_stage",
]

cancer_lines_words = [
    "hct-?116",
    "ht-?29",
    "rko",
    "caco-?2",
    "sw-?480",
    "sw-?620",
    "lovo",
    "mcf-?7",
    "mda-?mb-?\d+",
    "t-?47d",
    "hcc-?1954",
    "hela",
    "a549",
    "h1299",
    "h1437",
    "h1703",
    "h2126",
    "h2228",
    "h2347",
    "h23",
    "h460",
    "hepg2",
    "huh-?7",
    "snu-?\d+",
    "sw-?13",
    "k-?562",
    "jurkat",
    "hl-?60",
    "thp-?1",
    "u-?937",
    "molm-?13",
    "nalm-?6",
    "reh",
    "molt-?4",
    "kyse-?\d+",
    "sh-?sy5y",
    "sk-?n-?[a-z]{2}",
    "kns-?42",
    "lncap",
    "pc-?3",
    "du-?145",
    "vcap",
    "22rv1",
    "gp5d",
    "caki-?[12]",
    "786-?o",
    "769-?p",
    "detroit-?562",
    "cal-?27",
    "fadu",
    "scc-?4",
    "ovcar-?\d",
    "a2780",
    "cov434",
    "panc-?1",
    "mia ?paca",
    "u-?2 ?os",
    "saos-?2",
    "mg-?63",
    "a375",
    "u-?87",
    "ln-?229",
    "aml-?3",
    "bt-?245",
    "su-?dipg\w*",
    "oci-?ly\d+",
    "blaer\d?",
    "mel-?1",
    "akata",
    "gm-?ct",
]

cancer_lines_re = re.compile(
    r"(?<![a-z0-9])"
    + "("
    + "|".join(cancer_lines_words)
    + ")"
    + "(?![a-z0-9])"
    + "|"
    + "|".join(
        [
            "b-?all cell line",
            "leukemia cell line",
            "cancer cell line",
        ]
    ),
    re.I,
)

noncancer_lines_words = [
    "imr-?90",
    "wi-?38",
    "mrc-?5",
    "hek-?293t?",
    "293t?",
    "huvecs?",
    "hacat",
    "gm\d{5}",
    "na\d{5}",
    "hg\d{5}",
    "lcl",
    "hdfa",
    "bj fibroblast",
    "hbec3-kt",
    "saec",
    "luhmes",
    "lhcn-m2",
    "lhsar",
    "hpr1(-ar)?",
    "rpe-?1",
    "wtc-?11",
    "pgp1",
    "h1",
    "h9",
    "hues ?\d+",
    "ucla\d",
    "wibr\d",
    "rues\d",
    "khes\d",
    "hnes\d",
]

noncancer_lines_trailing_words = [
    "lymphoblastoid",
    "embryonic stem",
    r"\bh?escs?\b",
    "hipsc",
    "ipsc",
    r"\bips cells?\b",
    "keratinocyte",
    "fibroblast",
]

noncancer_lines_re = re.compile(
    r"(?<![a-z0-9])" + "("
    + "|".join(noncancer_lines_words)
    + ")"
    + "(?![a-z0-9])"
    + "|"
    + "("
    +"|".join(noncancer_lines_trailing_words)
    + ")",
    re.I,
)

cancer_neg_paren_words = [
    "?:no",
    "without",
    "negative for",
    "free of",
]

cancer_neg_paren_re = re.compile(
    r"\(" + "(" + "|".join(cancer_neg_paren_words) + ")" + r"\b[^)]*" + "\)",
    re.I,
)

cancer_neg_clause_prefix_words = [
    "?:no",
    "without",
    "negative for",
    "free of",
    "absence of"
]

cancer_neg_clause_suffix_words = [
    "?:cancer",
    "tumou?r",
    "malignan\w*",
    "carcinom\w*",
]

cancer_neg_clause_re = re.compile(
    "("
    + "|".join(cancer_neg_clause_prefix_words)
    + ")"
    + r"\b[^,;.)]*?"  # symbols
    + "("
    + "|".join(cancer_neg_clause_suffix_words)
    + ")"
    + "[^,;.)]*",  # symbols
    re.I,
)
