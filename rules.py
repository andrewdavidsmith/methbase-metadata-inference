# SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

"""
Hard-coded rules and regex for classifying sample metadata of Methbase
methylomes according to tissue of origin; designed for human samples.  The
"tissues" aren't necessarily tissues.

Also includes a mapping of tissues to UBERON terms.
"""

import re

# (label, regex): ordered, most specific first
rules = [
 ("cord blood",        r"cord ?blood|\bcbmc"),
 ("umbilical cord",    r"umbilical cord"),
 ("placenta",          r"placenta|planceta|trophoblast|chorion|villous|decidua"),
 ("bone marrow",       r"bone marrow|\bbm\b|marrow"),
 ("leukemia/lymphoma", r"leukemi|leukaemi|lymphoma|myeloma|\bk-?562|jurkat|hl-?60|thp-?1|\bu937|molt-?4|\bcll\b|\baml\d?\b|\breh\b|nalm-?6|b-?all\b|blaer"),
 ("blood",             r"(?<![a-z])blood|pbmc|pbl\b|\bwbc|white blood|leukocyte|leucocyte|lymphocyte|\bplasma\b|serum|buffy|leukapheresis|monocyte|neutrophil|granulocyte|\bt[- ]?cells?|car-?t\b|nkt|hspc|macrophage|\bgm12878|\bb[- ]?cell|b-lymphocyte|nk cell|cd[348]\+?|eosinophil|platelet|erythro|hematopoietic|haematopoietic|hsc\b|lcl\b|lymphoblastoid"),
 ("brain",             r"luhmes|spinal cord|sh-?sy5y|sk-?n-?[a-z]{2}|neuroblastoma|brain|cortex|hippocamp|amygdala|cerebell|cerebr|nucleus accumbens|nacc|striatum|thalamus|hypothalam|neuron|neural|astrocyte|microglia|glia|frontal|temporal lobe|parietal|occipital|\bba\d+|substantia|midbrain|glioma|glioblastoma|dlpfc|prefrontal|cingulate|pons|medulla|putamen|caudate|olfactory bulb"),
 ("stool",             r"stool|fec(es|al)|faeces|feces|meconium"),
 ("csf",               r"\bcsf\b|cerebrospinal"),
 ("sperm",             r"sperm|semen|spermatozoa"),
 ("testis",            r"testis|testes|testicular"),
 ("ovary",             r"ovary|ovarian|oocyte|granulosa|follic"),
 ("fallopian tube",    r"fallopian|salpinx|salpingectomy|oviduct"),
 ("uterus/cervix",     r"endometri|uterus|uterine|cervi|myometri|\bhela\b"),
 ("liver",             r"liver|\bhcc\b(?!-?\d)|hepat|hepg2|huh-?7|cholangio|biliary"),
 ("lung",              r"\bat2\b|\blung|pulmonary|bronch|alveol|\ba549|h1299|h1437|h1703|h2126|h2228|h2347|\bh\d{3,4}\b lung|nsclc|imr-?90|wi-?38|mrc-?5|airway"),
 ("small intestine",   r"small intestine|ileum|jejunum|duoden|\bileal"),
 ("colon",             r"colon|\brko\b|gp5d|colorectal|rectal|rectum|caecum|cecum|intestin|sigmoid|hct-?116|caco-?2|ht-?29|sw480|sw620|lovo|crc\b"),
 ("esophagus",         r"esophag|oesophag|escc|kyse-?\d+"),
 ("stomach",           r"stomach|gastric|\bmkn|ags\b"),
 ("pancreas",          r"pancrea|panc-?1|mia ?paca|pdac|islet"),
 ("kidney",            r"kidney|renal|nephr|caki-?[12]|hek-?293|293t?\b|\b786-?o"),
 ("bladder",           r"bladder|urothel|urinary|\bt24\b|urine"),
 ("prostate",          r"vcap|lhsar|prostat|lncap|\bpc-?3\b|du145|22rv1"),
 ("breast",            r"breast|mammary|hcc-?1954|mcf-?7|mda-?mb|\bt47d|\bbt-?474|sk-?br-?3"),
 ("skin",              r"skin|dermis|dermal|epiderm|keratinocyte|fibroblast|melanocyte|melanoma|foreskin|\bhacat|epidermal|scalp|hair"),
 ("heart",             r"heart|cardiac|cardiomyocyte|ventricle|atrium|atrial|myocard"),
 ("muscle",            r"muscle|myoblast|myotube|myofiber|skeletal|\bc2c12|vastus|soleus|psoas"),
 ("adipose",           r"adipos|\bfat\b|adipocyte|omentum"),
 ("nerve",             r"nerve|neural crest ganglia|ganglion"),
 ("thyroid",           r"thyroid"),
 ("adrenal gland",     r"adrenal"),
 ("spleen",            r"spleen|splenic"),
 ("lymph node",        r"lymph node|lymphoid"),
 ("thymus",            r"thymus|thymic"),
 ("synovium",          r"synovi"),
 ("bone",              r"\bbone\b|tibia|osteo|chondro|cartilage|\bsaos|u-?2 ?os|\bmg-?63|periosteum|skull|calvari"),
 ("oral cavity",       r"oral|mouth|tongue|saliva|buccal|gingiv|lip tissue|\blip\b|oscc|head and neck|hnscc|tonsil|pharyn|larynx|laryng"),
 ("nasal",             r"nasal|nasopharyn|nose|sinus"),
 ("eye",               r"\beye|retina|retinal|cornea|lens|ocular|rpe\b|uveal"),
 ("ear",               r"cochlea|inner ear|\bear\b"),
 ("stem/embryonic cells", r"\bh?escs?\b|hescs?|ipsc|\bips\b|ips cell|sees\d|tslprimed|imelc|\bte\d+|trophectoderm|\bh?ipscs?\b|\bpgp1\b|\begc?s?\b|embryonic germ|\bicm\b|inner cell mass|reprogramming|oogonia|embryonic stem|pluripotent|embryo|blastocyst|\bhues\d|ucla\d|epilc|pgc|primordial germ|germ cell|epiblast|gastrulo|naive|hes\d|hksr|msc|mesenchymal|stem cell|organoid|\bh1\b|\bh9\b"),
 ("gonad",             r"gonad"),
 ("vasculature",       r"endotheli|\bhuvecs?|\bh?ecs?\b|vsmc|smooth muscle|huvec|aorta|carotid|artery|arterial|vein|venous|vascular|umbilical vein"),
 # ("tumor (unspecified)", r"xenograft|\bcafs?\b|tumou?r|carcinoma|cancer|neoplas|sarcoma|metasta"),
]
rules_re = [(l, re.compile(r, re.I)) for l, r in rules]

# Field priority: most direct tissue evidence first
evidence_priority = [
    r"tissue",
    r"tissue_type",
    r"tissue_source",
    r"isolation_source",
    r"body_site",
    r"organ",
    r"sample_type",
    r"source_name",
    r"cell_type",
    r"cell_line",
    r"sample_comment",
    r"sample_name",
    r"diagnosis",
    r"disease",
    r"brain_region",
    r"common_name",
    r"sample_title",
    r"title",
    r"description",
    r"phenotype",
    r"BIOMATERIAL_TYPE",
    r"dev_stage",
    r"health_state",
    r"cancer_status",
    r"group",
    r"study",
    r"sample_set",
]

nullish_words = [
    r"na",
    r"n/a",
    r"n\.a\.",
    r"none",
    r"null",
    r"missing",
    r"not applicable",
    r"not collected",
    r"not available",
    r"unknown",
    r"unspecified",
    r"-",
    r"\?",
]
nullish_re = re.compile(r"^(" + "|".join(nullish_words) + r")$", re.I)

uninformative_words = [
    r"na",
    r"n/a",
    r"none",
    r"null",
    r"missing",
    r"not applicable",
    r"not collected",
    r"not available",
    r"unknown",
    r"unspecified",
    r"-",
    r"\?",
    r"primary",
    r"cell line",
    r"culture cells?",
    r"cells?",
    r"tissue",
    r"normal",
    r"control",
    r"n\.a\.",
    r"ns",
    r"dna",
    r"whole organism",
    r"ln",
    r"human somatic cells",
    r"matched normal",
    r"adjacent normal",
    r"normal adjacent tissue",
    r"adjacent normal tissue",
]

uninformative_re = re.compile(r"^(" + "|".join(uninformative_words) + r")$", re.I)

administrative_words = [
    r"provider",
    r"center",
    r"centre",
    r"^lab$",
    r"submitter",
    r"institut",
    r"geo_loc",
    r"collected_by",
    r"^INSDC",
    r"^ENA",
    r"external_id",
    r"^study$",
    r"batch",
    r"sample_set",
]
admin_re = re.compile("(" + "|".join(administrative_words) + ")", re.I)

tech_words = [
    r"bs-?seq",
    r"oxbs",
    r"tab-?seq",
    r"methylation",
    r"spike",
    r"lambda",
    r"unmethylated",
    r"methylated",
    r"conversion",
    r"input",
    r"mock",
    r"no[- ]template",
    r"water",
]

control_words = r"[ _-]*(control|ctrl)|control[ _-]*(dna|oligo)"

tech_re = re.compile("(" + "|".join(tech_words) + ")" + control_words, re.I)

culture_hint_words = [
    r"cell line",
    r"cultured",
    r"culture",
    r"organoid",
    r"ipsc",
    r"\bips\b",
    r"\bh?escs?\b",
    r"hescs",
    r"xenograft",
]

culture_hint_re = re.compile("(" + "|".join(culture_hint_words) + ")", re.I)

cell_line_hint_words = [
    r"cell line",
    r"\bline\b",
]

cell_line_hint_re = re.compile("(" + "|".join(cell_line_hint_words) + ")", re.I)

# Stored UBERON mapping (tissue label -> UBERON id, name) Labels that are not
# anatomical structures (stem cells, tumors, leukemia) do not have a
# corresponding UBERON term.
uberon_term_map = {
    r"blood": ("UBERON:0000178", "blood"),
    r"cord blood": ("UBERON:0012168", "umbilical cord blood"),
    r"umbilical cord": ("UBERON:0002331", "umbilical cord"),
    r"placenta": ("UBERON:0001987", "placenta"),
    r"bone marrow": ("UBERON:0002371", "bone marrow"),
    r"brain": ("UBERON:0000955", "brain"),
    r"csf": ("UBERON:0001359", "cerebrospinal fluid"),
    r"sperm": ("UBERON:0006530", "seminal fluid"),
    r"testis": ("UBERON:0000473", "testis"),
    r"ovary": ("UBERON:0000992", "ovary"),
    r"fallopian tube": ("UBERON:0003889", "fallopian tube"),
    r"uterus/cervix": ("UBERON:0000995", "uterus"),
    r"liver": ("UBERON:0002107", "liver"),
    r"lung": ("UBERON:0002048", "lung"),
    r"small intestine": ("UBERON:0002108", "small intestine"),
    r"colon": ("UBERON:0001155", "colon"),
    r"esophagus": ("UBERON:0001043", "esophagus"),
    r"stomach": ("UBERON:0000945", "stomach"),
    r"pancreas": ("UBERON:0001264", "pancreas"),
    r"kidney": ("UBERON:0002113", "kidney"),
    r"bladder": ("UBERON:0001255", "urinary bladder"),
    r"prostate": ("UBERON:0002367", "prostate gland"),
    r"breast": ("UBERON:0000310", "breast"),
    r"skin": ("UBERON:0002097", "skin of body"),
    r"heart": ("UBERON:0000948", "heart"),
    r"muscle": ("UBERON:0001630", "muscle organ"),
    r"adipose": ("UBERON:0001013", "adipose tissue"),
    r"nerve": ("UBERON:0001021", "nerve"),
    r"thyroid": ("UBERON:0002046", "thyroid gland"),
    r"adrenal gland": ("UBERON:0002369", "adrenal gland"),
    r"spleen": ("UBERON:0002106", "spleen"),
    r"lymph node": ("UBERON:0000029", "lymph node"),
    r"thymus": ("UBERON:0002370", "thymus"),
    r"synovium": ("UBERON:0002018", "synovial membrane of synovial joint"),
    r"bone": ("UBERON:0001474", "bone element"),
    r"oral cavity": ("UBERON:0000167", "oral cavity"),
    r"nasal": ("UBERON:0000004", "nose"),
    r"eye": ("UBERON:0000970", "eye"),
    r"ear": ("UBERON:0001690", "ear"),
    r"gonad": ("UBERON:0000991", "gonad"),
    r"vasculature": ("UBERON:0002049", "vasculature"),
    r"stool": ("UBERON:0001988", "feces"),
}
