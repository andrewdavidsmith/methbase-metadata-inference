# SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

"""
Hard-coded rules and regex for classifying sample metadata of Methbase
methylomes according to tissue of origin; designed for MOUSE samples.  The
"tissues" aren't necessarily tissues.

Also includes a mapping of tissues to UBERON terms.
"""

import re

# checked after organ rules, before anything else
mouse_rules_embryo_words = [
    r"embryo",
    r"zygote",
    r"\b[1-9]-?cell\b",
    r"\b(two|four|eight)[- ]cell",
    r"\b[248]c\b",
    r"morula",
    r"blastula",
    r"blastocyst",
    r"\bicm\b",
    r"inner cell mass",
    r"trophectoderm",
    r"\bte\b",
    r"epiblast",
    r"ectoderm",
    r"endoderm",
    r"mesoderm",
    r"gastrula",
    r"preimplantation",
    r"pre-implantation",
    r"\be\d+(\.\d+)?\b",
    r"embryonic day",
    r"\bscnt\b",
    r"nuclear transfer",
    r"yolk sac",
    r"extraembryonic",
    r"\bexe\b",
]

mouse_rules_first = [   # checked before the core rules
    ("kidney",              r"kidney cortex|renal cortex"),
    ("adrenal gland",       r"adrenal cortex"),
    ("muscle",              r"tibialis|gastrocnemius|soleus|quadriceps|plantaris|extensor digitorum|\bedl\b"),
    ("testis",              r"sertoli|leydig|spermatogon|spermatocyte|spermatid|seminiferous"),
    ("ovary",               r"cumulus|granulosa"),
    ("small intestine",     r"small intest|duoden|jejun|ile(um|al)\b|enterocyte|villi|villus"),
    ("colon",               r"colon|colorectal|rectum|caecum|cecum|large intestine|(?<![a-z0-9])(ct-?26|mc-?38|colon-?26)(?![a-z0-9])"),
    ("intestine (unspecified)", r"intestin|crypt|paneth|goblet"),
    ("nasal",               r"olfactory epithel"),
    ("ear",                 r"organ of corti|hair cells?"),
    ("oocyte",              r"oocyte|\bmii\b|\bgv\b|germinal vesicle|\bova\b|\bovum"),
    ("stem/embryonic cells", r"\bm?escs?\b|\bes cells?\b|embryonic stem|\bm?ipscs?\b|\bips\b|\bepilcs?\b|\bepisc|pgclc|primordial germ|\bpgcs?\b|reprogramming|pluripotent|\bnaive\b|\bprimed\b|\btscs?\b|trophoblast stem|\bxen\b|\be14\b|\br1 129|\bcgr8|\bj1\b|\bv6\.?5|\bhm1\b|\bkh2\b|\bp19\b|\bf9\b|embryonic germ|\begc|haesc|ipsc|(?<![a-z0-9])mes(?![a-z0-9])"),
    ("germ cells",          r"germ cells?"),
    ("embryonic fibroblast", r"(?<![a-z])mefs?(?![a-z])|embryonic fibroblast|nih-?3t3|\b3t3\b|3t3-?l1"),
    ("neural tube",         r"neural tube|neuroepithel|headfold|head fold"),
    ("limb",                r"\blimb|forelimb|hindlimb|limb bud"),
    ("face",                r"facial prominence|craniofacial|palate|mandib|maxill|branchial arch|pharyngeal arch"),
    ("retina/photoreceptor", r"retina|retinal|photoreceptor|horizontal cell|bipolar cell|amacrine|ganglion cell|rod\b|\bcone\b|m[uü]ller"),
    ("bone marrow",         r"\bwbm\b|whole[ _]bone[ _]marrow|\bksl\b|\blsk\b|\bhsc|\bhspc|\bmpp|\bclp\b|\bcmp\b|\bgmp\b|\bmep\b|lin-?\s*sca|hematopoietic progenitor|haematopoietic progenitor|myeloid progenitor|lymphoid progenitor"),
    ("lymph node",          r"(?<![a-z])ln[_ -]|lymph node"),
    ("liver",               r"hepa ?1-?6|hepatocyte|hepatoma|kupffer"),
    ("brain",               r"neuro-?2a|(?<![a-z])n2a(?![a-z])|neuroblastoma|granule cell|purkinje|oligodendro|\bopc\b|dentate gyrus|subventricular|habenula|\bsvz\b"),
    ("blood",               r"tconv|treg|(?<![a-z])t[- ]?cells?|(?<![a-z])b[- ]?cells?|nk cells?|macrophage|monocyte|neutrophil|dendritic|raw ?264|j774|myeloid|splenocyte(?!.*spleen)"),
    ("whisker/tail",  r"\btail\b|whisker"),
]

mouse_rules_embryo = [
    ("embryo", "|".join(mouse_rules_embryo_words)),
]

# (label, regex) -- ordered, most specific first
rules = mouse_rules_first + [
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
    ("gonad",             r"gonad"),
    ("vasculature",       r"endotheli|\bhuvecs?|\bh?ecs?\b|vsmc|smooth muscle|huvec|aorta|carotid|artery|arterial|vein|venous|vascular|umbilical vein"),
    ("embryo", "|".join(mouse_rules_embryo_words)),
    ("tumor (unspecified)", r"xenograft|\bcafs?\b|tumou?r|carcinoma|cancer|neoplas|sarcoma|metasta"),
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
    r"breast": ("UBERON:0001911", "mammary gland"),
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
    r"embryo": ("UBERON:0000922", "embryo"),
    r"neural tube": ("UBERON:0001049", "neural tube"),
    r"limb": ("UBERON:0002101", "limb"),
    r"face": ("UBERON:0001456", "face"),
    r"retina/photoreceptor": ("UBERON:0000966", "retina"),
    r"whisker/tail": ("", ""),
    r"intestine (unspecified)": ("UBERON:0000160", "intestine"),
}

# genotype/strain/treatment/etc. hold gene names, diets and strains -> never
# tissue evidence
exclude_words = [
    "genotype",
    "strain",
    "background",
    "background_strain",
    "strain_background",
    "breed",
    "ecotype",
    "cultivar",
    "treatment",
    "diet",
    "sex",
    "age",
    "donor_age",
    "lab",
    "litter",
    "replicate",
    "biological_replicate",
    "time",
    "timepoint",
    "isolate",
    "subject_id",
    "donor_ID",
    "growth_protocol",
    "collection_date",
    "dev_stage_note",
]

exclude_re = re.compile(r"^(" + "|".join(exclude_words) + r")$", re.I)

strain_rules = [
    ("unannotated", r"(missing|unknown)?"),
    # ("Cell line", r"(cell|HT22|E14\w*|J1|R1|NJ2008)"),
    # ("B6D2F1", r"B6D2(F1)?"),
    ("B6C3F1", r"(B6C3F1(/J)?|C3B6F\d|B6C3\(Cg\))"),
    ("B6x129 hybrid", r"(V6\.?5|KH2|C57BL6?J?/129|C57BL/6J-129Sv|FVB/129)"),
    (
        "F1/mixed", # " (other)",
        r".*(\bF1\b|hybrid|mixed|offspring|\)\*|\(MAT\)|J\*DBA|"
        r"C57BL/6xCBA|Spretus|CZECHII|WSB|CB6F1|CastF1|BALB/cfC3H).*|"
        r"M[BP][SP](x\w+)?",
    ),
    ("C3H", r".*(C3H|^C3\.Cg).*"),
    (
        "transgenic/KI/KO",
        r".*(tm\d|Tg\(|cre|GFP|Gfp|mRFP|mCherry|rtTA|Trap2|OG2|"
        r"Tet1|Pvalb|Dnmt3a|Y\^|ob/ob|Dp\(\d+\)).*",
    ),
    ("129", r"129.*"),
    # ("DBA/2", r"D(BA|AB)/2J?"),
    ("BALB/c", r"BALB/c"),
    ("FVB", r"FVB(/NJ?)?"),
    ("NOD", r"NOD"),
    ("A/J", r"A/J"),
    ("NZO", r"NZO.*"),
    # ("SAM", r"SAM[PR]\d*"),
    ("CD-1/ICR", r"(CD-?1|ICR)(-ICR)?"),
    ("CAST", r"CAST/EiJ"),
    # ("PWK", r"PWK/PHJ"),
    # ("SPRET", r"SPRET.*"),
    ("C57BL/6", r"(B6(n_\d)?|C5[67].*|CD57BL6|Bl/?6)"),
]

strain_rules_re = [
    (label, re.compile(pattern, re.I)) for label, pattern in strain_rules
]

include_fields_for_strain = [
    "strain",
    "background_strain",
    "hybrid_strain",
    "genotype",
    "genetic_background",
    "background",
    "breed",
    "litter",
    ## ADS: source_name includes the wrong things too often
    # "source_name",
]

def get_strain(metadata: dict) -> str:
    """
    Use metadata to identify strain for mouse
    """
    for field in include_fields_for_strain:
        if field not in metadata:
            continue
        value = metadata[field]
        if not isinstance(value, str):
            continue
        for label, rx in strain_rules_re:
            if rx.fullmatch(value):
                return label, f"{field}: {value}"
    return "other", ""
