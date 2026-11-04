"""Acquire paired CPTAC ITGB1 RNA and mass-spectrometry protein data."""

from __future__ import annotations

import json

import pandas as pd
import requests

from scripts.core.cli import parse_project_root


API = "https://www.cbioportal.org/api"
ENTREZ_ITGB1 = 3688
STUDIES = [
    "paad_cptac_2021", "luad_cptac_2020", "lusc_cptac_2021",
    "coad_cptac_2019", "ucec_cptac_2020", "gbm_cptac_2021",
    "brca_cptac_2020",
]


def get_json(path, params=None):
    response = requests.get(API + path, params=params, timeout=180)
    response.raise_for_status()
    return response.json()


def select_profiles(study):
    profiles = get_json(f"/studies/{study}/molecular-profiles")
    protein = [p for p in profiles if p.get("molecularAlterationType") == "PROTEIN_LEVEL"
               and p.get("datatype") != "Z-SCORE"]
    rna = [p for p in profiles if p.get("molecularAlterationType") == "MRNA_EXPRESSION"
           and p.get("datatype") == "CONTINUOUS" and "mirna" not in p.get("name", "").lower()]
    if not protein or not rna:
        return None, None
    rna.sort(key=lambda p: ("mrna" not in p.get("name", "").lower(),
                            "rsem" not in p.get("description", "").lower()))
    return protein[0], rna[0]


def fetch_gene(profile_id, study):
    return get_json(f"/molecular-profiles/{profile_id}/molecular-data", {
        "entrezGeneId": ENTREZ_ITGB1, "sampleListId": f"{study}_all",
    })


def main():
    root, _ = parse_project_root()
    raw = root / "data/raw/CPTAC_ITGB1"
    processed = root / "data/processed/CPTAC_ITGB1"
    raw.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)
    tables, provenance = [], []
    for study in STUDIES:
        protein_profile, rna_profile = select_profiles(study)
        if not protein_profile or not rna_profile:
            provenance.append({"study": study, "status": "no paired profiles"})
            continue
        protein = fetch_gene(protein_profile["molecularProfileId"], study)
        rna = fetch_gene(rna_profile["molecularProfileId"], study)
        (raw / f"{study}_protein.json").write_text(json.dumps(protein), encoding="utf-8")
        (raw / f"{study}_rna.json").write_text(json.dumps(rna), encoding="utf-8")
        protein_table = pd.DataFrame(protein)[["sampleId", "value"]].rename(
            columns={"value": "ITGB1_protein"})
        rna_table = pd.DataFrame(rna)[["sampleId", "value"]].rename(
            columns={"value": "ITGB1_RNA"})
        merged = rna_table.merge(protein_table, on="sampleId", how="inner")
        merged.insert(0, "study", study)
        merged["rna_profile"] = rna_profile["molecularProfileId"]
        merged["protein_profile"] = protein_profile["molecularProfileId"]
        tables.append(merged)
        provenance.append({
            "study": study, "status": "included", "n": len(merged),
            "rna_profile": rna_profile["molecularProfileId"],
            "protein_profile": protein_profile["molecularProfileId"],
        })
    if not tables:
        raise RuntimeError("No CPTAC study provided paired ITGB1 RNA/protein data")
    pd.concat(tables, ignore_index=True).to_csv(
        processed / "CPTAC_ITGB1_RNA_protein.tsv", sep="\t", index=False)
    pd.DataFrame(provenance).to_csv(
        processed / "CPTAC_ITGB1_profile_selection.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
