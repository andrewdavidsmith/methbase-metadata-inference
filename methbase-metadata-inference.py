#!/usr/bin/env python3
# SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith
"""
methbase-metadata-inference.py: rule-based tissue and cancer-status labeling
of SRA sample metadata for Methbase.

Input:  JSON of the form {"SRX...": {"metadata": {field: value, ...}}, ...}
Output: <out>/tissue_labels.csv        one row per accession (with uberon info)
        <out>/uberon_mapping.csv       label -> UBERON id/name table with counts
        <out>/labeled.json             same as input with label fields added
        <out>/review_unresolved.json   accessions where tissue = other
        <out>/<input>_tissue.json      identified tissue names with counts
Console: summary counts and fields NOT used by the cancer classifier.

Usage: methbase-metadata-inference.py --out results_dir input.json

Rules used for classification are hard-coded in source files.
"""

import argparse
import collections
import csv
import json
import os
import re
import sys

from rules import (
    rules_re,
    evidence_priority,
    uninformative_re,
    admin_re,
    tech_re,
    control_words,
    nullish_re,
    culture_hint_re,
    cell_line_hint_re,
    uberon_term_map,
)

from cancer_rules import (
    cancer_non_re,
    cancer_control_re,
    cancer_strong_re,
    cancer_growth_re,
    cancer_benign_re,
    cancer_weak_re,
    patient_donor_fields,
    sample_level_fields,
    tumor_fields,
    cancer_lines_re,
    noncancer_lines_re,
    cancer_neg_paren_re,
    cancer_neg_clause_re,
)


def count_by_col(rows, col):
    """
    Count the appearances of each value in a given col for a set of rows
    """
    return collections.Counter(row[col] for row in rows)


def strip_all(obj):
    """
    Recursively strips all 'values' and converts all to string if they are None.
    """
    if isinstance(obj, str):
        return "" if obj is None else obj.strip()
    if isinstance(obj, dict):
        return {k.strip(): strip_all(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [strip_all(v) for v in obj]
    return obj


def rules_match(text):
    """
    Match text to rules using 'rules_re' and return a 'label' based on the first
    matched rule or None. Note: rules must remain ordered.
    """
    for label, rx in rules_re:
        if rx.search(text):
            return label
    return None


def get_tissue_label(metadata):
    """
    Returns a 4-tuple of (label, field used, confidence, evidence)
    """
    if not metadata:
        return "other", None, "none", ""
    for meta in metadata.values():
        if not type(meta) is str:
            raise TypeError("Expected string: " + str(meta))
    # prioritized fields
    fields = [field for field in evidence_priority if field in metadata] + [
        field
        for field in metadata
        if field not in evidence_priority and not admin_re.search(field)
    ]
    for field in fields:
        curr_meta = metadata[field]
        if curr_meta is None:
            continue
        if not curr_meta or uninformative_re.match(curr_meta):
            continue
        if label := rules_match(curr_meta):
            ## validate positions in evidence priority
            confidence = (
                "high"
                if field in evidence_priority[:6]
                else "medium"
                if field in evidence_priority[:14]
                else "low"
            )
            return label, field, confidence, curr_meta
    return "other", None, "none", ""


def evidence_from_field(text, allow_control_sample=True):
    """
    Classify a text value (metadata for a field) and return (classification,
    matching_text) pair or (None, None). The first member of the pair is the
    classification {non-cancer, cancer, benign/premalignant, control} and the
    second is the matching text that formed the evidence.
    """
    if not text or nullish_re.match(text):
        return None
    t = cancer_neg_clause_re.sub(
        "", cancer_neg_paren_re.sub("", tech_re.sub("", text))
    )
    if cancer_non_re.search(t):
        return "non-cancer", cancer_non_re.search(t)[0]
    if cancer_strong_re.search(t):
        return "cancer", cancer_strong_re.search(t)[0]
    if cancer_growth_re.search(t):
        return "benign/premalignant", cancer_growth_re.search(t)[0]
    if cancer_benign_re.search(t):
        return "non-cancer", "benign"
    if cancer_weak_re.search(t):
        return "cancer", cancer_weak_re.search(t)[0]
    if allow_control_sample and cancer_control_re.search(t):
        return "control", cancer_control_re.search(t)[0]
    return (None, None)


def get_cancer_evidence(metadata, fields, allow_control_sample):
    """
    Iterate over fields and construct a list of all (field, metadata, outcome)
    for cancer classification; the combination of (field, metadata) is the
    evidence for the outcome.
    """
    hits = []
    for field in fields:
        if field in metadata:
            if r := evidence_from_field(metadata[field], allow_control_sample):
                hits.append((field, metadata[field], r[0]))
    return hits


def add_cancer_label(hits):
    """
    Precedence when adding labels:
    non-cancer (adjacent/normal/healthy) > cancer > benign growth > control.
    """
    cancer_status_labels = (
        "non-cancer",
        "cancer",
        "benign/premalignant",
        "control",
    )
    for label in cancer_status_labels:
        for field, value, cancer_status in hits:
            if cancer_status == label:
                return label, f"{field}: {value[:80]}"
    return None


def get_blob(metadata):
    """
    Project the metadata onto fields relevant for classifying samples by cancer
    status
    """
    relevant_fields = (
        "cell_line",
        "source_name",
        "cell_type",
        "sample_name",
        "tissue",
        "sample_type",
        "sample_comment",
        "sample_title",
        "title",
        "description",
    )
    return " ".join(metadata.get(field, "") for field in relevant_fields)


def get_is_cell_line(metadata, blob):
    """
    Determine if the metadata says this is a cell line
    """
    return bool(
        metadata.get("cell_line")
        and not nullish_re.match(metadata.get("cell_line"))
    ) or bool(cell_line_hint_re.match(blob))


def get_is_cultured(metadata, blob, cancer_status):
    """
    Determine if the metadata says the sample is cultured cells
    """
    if "cell_line" in metadata:
        return "yes", "metadata key: cell_line"
    if (evidence := culture_hint_re.search(blob)):
        return "yes", evidence[0]
    if cancer_status == "cancer cell line":
        return "yes", "cancer status: cancer cell line"
    return "no", "none"


def cancer_classification(metadata):
    """
    Returns a triple with (classification, confidence, evidence)
    """
    if not metadata:
        return "no evidence", "none", ""

    blob = get_blob(metadata)

    cancer_line = cancer_lines_re.search(blob)
    normal_line = noncancer_lines_re.search(blob)
    is_cell_line = get_is_cell_line(metadata, blob) or cancer_line

    # Tumor-specific fiels (currently all are staging) imply a tumor sample
    tumor_specific_field = next(
        (
            field
            for field in tumor_fields
            if field in metadata
            and not nullish_re.match(metadata[field])
            and not re.match(r"normal|non", metadata[field], re.I)
        ),
        None,
    )

    # Sample-level evidence
    sample_level_evidence = get_cancer_evidence(
        metadata, sample_level_fields, allow_control_sample=not is_cell_line
    )
    res = add_cancer_label(sample_level_evidence)
    if not res and tumor_specific_field:
        res = (
            "cancer",
            f"{tumor_specific_field}: {metadata[tumor_specific_field][:80]}",
        )
    if res:
        classification, evidence = res
        if classification == "control":
            return "non-cancer", "medium", evidence
        if classification == "cancer" and is_cell_line:
            return "cancer cell line", "high", evidence
        return classification, "high", evidence

    # Known cell lines (from hard-coded names)
    if cancer_line:
        evidence = f"known cancer-derived line: {cancer_line[0]}"
        return "cancer cell line", "medium", evidence
    if normal_line:
        evidence = f"known non-cancer line/cell type: {normal_line[0]}"
        return "non-cancer", "medium", evidence

    # Patient/donor-level evidence
    patient_level_evidence = get_cancer_evidence(
        metadata, patient_donor_fields, allow_control_sample=True
    )
    if res := add_cancer_label(patient_level_evidence):
        classification, evidence = res
        if classification == "control":
            return "non-cancer", "medium", evidence
        if classification == "cancer" and is_cell_line:
            return "cancer cell line", "medium", evidence
        return classification, "medium", evidence

    # If we see "disease=none" then we conclude non-cancer
    if metadata.get("disease", "").lower() == "none":
        return "non-cancer", "low", "disease: none"

    # If we still have no classification, report no evidence
    return "no evidence", "none", ""


def tissue_to_uberon(tissue_label):
    """
    Get the UBERON info for a given tissue label
    """
    return uberon_term_map.get(tissue_label, ("", ""))


def write_tissue_set(rows, outdir, assembly):
    """
    Write a JSON file containing the set of tissues and their counts
    """
    tissue_col = 1
    outpath = os.path.join(outdir, f"{assembly}_tissues.json")
    with open(outpath, "w", encoding="utf-8") as out:
        print(json.dumps(count_by_col(rows, tissue_col), indent=4), file=out)


def report_unused_for_cancer(data, output_labels):
    """
    Report fields the cancer classifier never attempts to use
    """
    read = (
        set(sample_level_fields) | set(patient_donor_fields) | set(tumor_fields)
    )
    field_count = collections.Counter(
        field
        for rec in data.values()
        if isinstance(rec, dict)
        for field in rec.get("metadata", {})
        if field not in read and not admin_re.search(field)
        and field not in output_labels
    )
    print()
    print("Fields NOT used by the cancer classifier (top 25):")
    for field, count in field_count.most_common(25):
        print(f"{count:8d} {field}")


def report_summary_stats(rows, unresolved):
    """
    Report summary stats to the terminal
    """
    n_tissues_to_report = 12
    tissue_col = 1
    cancer_col = 7
    uberon_id_col = 5
    tissue_confid_col = 2
    print(f"{len(rows)} accessions")
    tissues = count_by_col(rows, tissue_col).most_common(n_tissues_to_report)
    print("tissue: ", tissues, "...")
    print("cancer: ", dict(count_by_col(rows, cancer_col)))
    with_uberon_id = sum(1 for row in rows if row[uberon_id_col])
    print("with UBERON id:", with_uberon_id, "of", len(rows))
    tissue_confid = dict(count_by_col(rows, tissue_confid_col))
    print("tissue confidence:", tissue_confid)
    n_unknown = len(unresolved)
    print(f"unknown tissue w metadata (see review_unresolved.csv): {n_unknown}")


def write_tissue_labels_csv(outdir, rows, assembly):
    """
    Write a CSV file of tissue labels
    """
    tissue_labels_csv_header = [
        "accession",
        "tissue",
        "tissue_confidence",
        "tissue_source_field",
        "tissue_evidence",
        "uberon_id",
        "uberon_name",
        "cancer_status",
        "cancer_confidence",
        "cancer_evidence",
        "cultured",
        "cultured_evidence",
    ]
    tissue_labels_csv_out = os.path.join(
        outdir, f"{assembly}_tissue_labels.csv"
    )
    with open(tissue_labels_csv_out, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(tissue_labels_csv_header)
        writer.writerows(rows)


def main():
    """
    main
    """
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "infile",
        help="input JSON file",
    )
    ap.add_argument(
        "--assembly",
        required=True,
        help="name of the assembly"
    )
    ap.add_argument(
        "--out",
        required=True,
        help="output directory (default: current dir)"
    )
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    with open(args.infile, "r", encoding="utf-8") as infile:
        data = json.loads(infile.read())
        if not isinstance(data, dict):
            sys.exit("Expected a JSON object keyed by accession.")
        data = strip_all(data)

    output_labels = (
        "tissue",
        "tissue_confidence",
        "tissue_source_field",
        "tissue_evidence",
        "uberon_id",
        "uberon_name",
        "cancer_status",
        "cancer_confidence",
        "cancer_evidence",
        "cultured",
        "cultured_evidence",
    )

    blob_fields = (
        "cell_line",
        "cell_type",
        "source_name",
        "sample_type",
        "tissue",
    )

    results = {}
    rows = []
    unresolved = {}
    for accession, rec in data.items():
        if not isinstance(rec, dict):
            continue
        metadata = rec.get("metadata", {})
        (
            tissue_label,
            tissue_field,
            tissue_confidence,
            tissue_evidence,
        ) = get_tissue_label(metadata)
        (
            cancer_status,
            cancer_confidence,
            cancer_evidence,
        ) = cancer_classification(metadata)
        uberon_id, uberon_name = tissue_to_uberon(tissue_label)
        blob = " ".join(metadata.get(field, "") for field in blob_fields)
        is_cultured, cultured_evidence = get_is_cultured(metadata, blob, cancer_status)
        rec["metadata"] |= {
            "tissue": tissue_label,
            "tissue_confidence": tissue_confidence,
            "tissue_source_field": tissue_field,
            "tissue_evidence": tissue_evidence,
            "uberon_id": uberon_id,
            "uberon_name": uberon_name,
            "cancer_status": cancer_status,
            "cancer_confidence": cancer_confidence,
            "cancer_evidence": cancer_evidence,
            "cultured": is_cultured,
            "cultured_evidence": cultured_evidence,
        }
        results[accession] = {
            "tissue": tissue_label,
            "tissue_confidence": tissue_confidence,
            "tissue_source_field": tissue_field,
            "tissue_evidence": tissue_evidence,
            "uberon_id": uberon_id,
            "uberon_name": uberon_name,
            "cancer_status": cancer_status,
            "cancer_confidence": cancer_confidence,
            "cancer_evidence": cancer_evidence,
            "cultured": is_cultured,
            "cultured_evidence": cultured_evidence,
        }
        rows.append(
            [
                accession,
                tissue_label,
                tissue_confidence,
                tissue_field,
                tissue_evidence,
                uberon_id,
                uberon_name,
                cancer_status,
                cancer_confidence,
                cancer_evidence,
                is_cultured,
                cultured_evidence,
            ]
        )
        if tissue_label == "other" and metadata:
            unresolved[accession] = metadata

    write_tissue_labels_csv(args.out, rows, args.assembly)

    outfile = os.path.join(args.out, f"{args.assembly}_labeled.json")
    with open(outfile, "w", encoding="utf-8") as out:
        json.dump({args.assembly: results}, out, indent=4, ensure_ascii=False)

    review_file = os.path.join(
        args.out, f"{args.assembly}_review_unresolved.json"
    )
    with open(review_file, "w", newline="", encoding="utf-8") as out:
        json.dump(unresolved, out, indent=4, ensure_ascii=False)

    uberon_mapping_out = os.path.join(
        args.out, f"{args.assembly}_uberon_mapping.csv"
    )
    with open(uberon_mapping_out, "w", newline="", encoding="utf-8") as out:
        tissue_col = 1
        headings = ("tissue", "uberon_id", "uberon_name", "n_accessions")
        writer = csv.writer(out)
        writer.writerow(headings)
        counted = count_by_col(rows, tissue_col)
        for label in sorted(counted, key=lambda x: -counted[x]):
            writer.writerow([label, *tissue_to_uberon(label), counted[label]])

    write_tissue_set(rows, args.out, args.assembly)
    report_summary_stats(rows, unresolved)
    report_unused_for_cancer(data, output_labels)


if __name__ == "__main__":
    main()
