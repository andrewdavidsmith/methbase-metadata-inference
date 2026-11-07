# SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

"""
Hard-coded rules and regex for classifying sample metadata of Methbase
methylomes according to cancer status; designed for human samples.
"""

import re

cancer_non_words = [
    r"adjacent",
    r"para-?cancer",
    r"non-?cancer",
    r"cancer-?free",
    r"non-?tumou?r",
    r"non-?malignan",
    r"non-?neoplas",
    r"matched normal",
    r"\bnormal\b",
    r"healthy",
    r"unaffected",
    r"disease-?free",
    r"no known",
]

cancer_non_re = re.compile("(" + "|".join(cancer_non_words) + ")", re.I)

cancer_control_words = [
    r"\bcontrols?\b",
    r"\bctrl\b",
    r"\bctr\b",
    r"_control",
    r"control_",
]

cancer_control_re = re.compile("(" + "|".join(cancer_control_words) + ")", re.I)

cancer_strong_words = [
    r"carcinom",
    r"cancer",
    r"malignan",
    r"metasta",
    r"sarcoma",
    r"leuk[a]?emi",
    r"lymphoma",
    r"myeloma",
    r"melanoma",
    r"blastoma",
    r"glioma",
    r"mesothelioma",
    r"seminoma",
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
    r"polyp",
    r"hyperplas",
    r"dysplas",
    r"metaplas",
    r"panin",
    r"leiomyoma",
    r"fibroid"
    r"\bnevus\b"
    r"\bnevi\b"
    r"barrett",
    r"precancer",
    r"pre-?malignan",
    r"non-?invasive neoplas",
    r"benign tumou?r",
    r"benign neoplas",
    r"benign growth",
]

cancer_growth_re = re.compile("(" + "|".join(cancer_growth_words) + ")", re.I)

cancer_benign_words = [
    r"\bbenign\b",
]

cancer_benign_re = re.compile("|".join(cancer_benign_words), re.I)

cancer_weak_words = [
    r"tumou?r",
    r"neoplas",
    r"xenograft",
    r"\bpdx\b",
]

cancer_weak_re = re.compile("(" + "|".join(cancer_weak_words) + ")", re.I)

# sample-level fields (describe the sample itself)
sample_level_fields = [
    r"cancer_status",
    r"sampling_site",
    r"sample_site",
    r"biopsy_site",
    r"tissue",
    r"tissue_type",
    r"sample_type",
    r"source_name",
    r"sample_condition",
    r"tumor_stage",
    r"sample_name",
    r"cell_type",
    r"cell_line",
    r"sample_title",
    r"title",
    r"description",
    r"sample_comment",
    r"source",
    r"BIOMATERIAL_TYPE",
    r"sample_group",
]

# patient/donor-level fields
patient_donor_fields = [
    r"disease",
    r"diagnosis",
    r"patient_diagnosis",
    r"health_state",
    r"DONOR_HEALTH_STATUS",
    r"subject_status",
    r"phenotype",
    r"condition",
    r"disease_type",
    r"pathological_type",
    r"group",
    r"case_control",
    r"status_of_low-risk_mds",
]

# tumor fields
tumor_fields = [
    r"tumor_stage",
    r"tnm_stage",
    r"tumor_grading",
    r"tumor_pathological_stage",
    r"tumor_source_tissue",
    r"disease_staging",
    r"developmental/tumor_stage",
]

cancer_lines_words = [
    r"hct-?116",
    r"ht-?29",
    r"rko",
    r"caco-?2",
    r"sw-?480",
    r"sw-?620",
    r"lovo",
    r"mcf-?7",
    r"mda-?mb-?\d+",
    r"t-?47d",
    r"hcc-?1954",
    r"hela",
    r"a549",
    r"h1299",
    r"h1437",
    r"h1703",
    r"h2126",
    r"h2228",
    r"h2347",
    r"h23",
    r"h460",
    r"hepg2",
    r"huh-?7",
    r"snu-?\d+",
    r"sw-?13",
    r"k-?562",
    r"jurkat",
    r"hl-?60",
    r"thp-?1",
    r"u-?937",
    r"molm-?13",
    r"nalm-?6",
    r"reh",
    r"molt-?4",
    r"kyse-?\d+",
    r"sh-?sy5y",
    r"sk-?n-?[a-z]{2}",
    r"kns-?42",
    r"lncap",
    r"pc-?3",
    r"du-?145",
    r"vcap",
    r"22rv1",
    r"gp5d",
    r"caki-?[12]",
    r"786-?o",
    r"769-?p",
    r"detroit-?562",
    r"cal-?27",
    r"fadu",
    r"scc-?4",
    r"ovcar-?\d",
    r"a2780",
    r"cov434",
    r"panc-?1",
    r"mia ?paca",
    r"u-?2 ?os",
    r"saos-?2",
    r"mg-?63",
    r"a375",
    r"u-?87",
    r"ln-?229",
    r"aml-?3",
    r"bt-?245",
    r"su-?dipg\w*",
    r"oci-?ly\d+",
    r"blaer\d?",
    r"mel-?1",
    r"akata",
    r"gm-?ct",
]

cancer_lines_re = re.compile(
    r"(?<![a-z0-9])"
    + "("
    + "|".join(cancer_lines_words)
    + ")"
    + r"(?![a-z0-9])"
    + "|"
    + "|".join(
        [
            r"b-?all cell line",
            r"leukemia cell line",
            r"cancer cell line",
        ]
    ),
    re.I,
)

mouse_noncancer_lines_words_1 = [
    r"nih-?3t3",
    r"3t3(-?l1)?",
    r"c2c12",
    r"mefs?",
    r"e14(tg2a)?",
    r"r1 (mesc|es\w*|129)",
    r"j1",
    r"cgr8",
    r"v6\.?5",
    r"hm1",
    r"kh2",
    r"mesc",
    r"l929",
    r"nrk",
    r"mc3t3",
]

mouse_noncancer_lines_words_2 = [
    r"\bm?escs?\b",
    r"haesc",
    r"ipsc",
    r"embryonic stem",
    r"\bm?ipscs?\b",
    r"mouse embryonic fibroblast",
]

noncancer_lines_re = re.compile(
    r"(?<![a-z0-9])(" +
    "|".join(mouse_noncancer_lines_words_1) +
    r")(?![a-z0-9])" + "|" +
    "|".join(mouse_noncancer_lines_words_2),
    re.I
)

cancer_neg_paren_words = [
    r"?:no",
    r"without",
    r"negative for",
    r"free of",
]

cancer_neg_paren_re = re.compile(
    r"\(" + "(" + "|".join(cancer_neg_paren_words) + ")" + r"\b[^)]*" + r"\)",
    re.I,
)

cancer_neg_clause_prefix_words = [
    r"?:no",
    r"without",
    r"negative for",
    r"free of",
    r"absence of"
]

cancer_neg_clause_suffix_words = [
    r"?:cancer",
    r"tumou?r",
    r"malignan\w*",
    r"carcinom\w*",
]

cancer_neg_clause_re = re.compile(
    r"("
    + "|".join(cancer_neg_clause_prefix_words)
    + ")"
    + r"\b[^,;.)]*?"  # symbols
    + "("
    + "|".join(cancer_neg_clause_suffix_words)
    + ")"
    + r"[^,;.)]*",  # symbols
    re.I,
)

cancer_lines_words_1 = [
    r"4t1",
    r"b16(-?f\d+)?",
    r"mc-?38",
    r"ct-?26",
    r"llc",
    r"lewis lung",
    r"el-?4",
    r"hepa ?1-?6",
    r"neuro-?2a",
    r"n2a",
    r"raw ?264\.?7",
    r"j774",
    r"p388",
    r"l1210",
    r"wehi-?\d+",
    r"a20",
    r"tramp-?c\d",
    r"min6",
    r"mel",
    r"kpc\w*",
    r"e0771",
    r"renca",
    r"py8119",
    r"b16-?f10",
    r"eo771",
    r"c1498",
    r"ba/f3",
    r"pan02",
    r"colon-?26",
    r"hepa1c1c7",
]

cancer_lines_words_2 = [
    r"cancer cell line",
    r"tumou?r cell line",
    r"leukemia cell line",
]

cancer_lines_re = re.compile(
    r"(?<![a-z0-9])(" +
    "|".join(cancer_lines_words_1) +
    r")(?![a-z0-9])" + "|" +
    "|".join(cancer_lines_words_2),
    re.I
)
