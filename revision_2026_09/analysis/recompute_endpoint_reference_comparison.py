#!/usr/bin/env python3
"""Compare corrected targeted-metabolite contrasts using 0 h and late controls."""

from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


WORKBOOK = Path("/Users/mingfeichen/Downloads/New_metabolites_corrected.xlsx")
SHEET = "Supplementary Table X Necromas "
OUTDIR = Path("/Users/mingfeichen/Manuscript/outputs/manual-20260921-a1/endpoint_reference_reanalysis")
OUTDIR.mkdir(parents=True, exist_ok=True)


def parse_sample(sample: str) -> tuple[str, str]:
    species = "Bacillus" if "Bacillus" in sample else "Rhodanobacter"
    if "-0hr_" in sample:
        return species, "0h"
    if species == "Bacillus":
        if "-8hr_" in sample:
            return species, "mid"
        if "-24hr-1mMAlCl3" in sample:
            return species, "Al"
        if "-24hr-20ulmLphage" in sample:
            return species, "Phage"
        if "-24hr-500ugmLKana" in sample:
            return species, "Kana"
        if "-24hr_" in sample:
            return species, "late"
    else:
        if "-48hr-1mMAlCl3" in sample:
            return species, "Al"
        if "-48hr-40ulmLphage" in sample:
            return species, "Phage"
        if "-48hr-500ugmLKana" in sample:
            return species, "Kana"
        if "-48hr_" in sample:
            return species, "late"
        if "-24hr_" in sample:
            return species, "mid"
    return species, "unknown"


def compare_feature(treatment: pd.Series, reference: pd.Series) -> tuple[float, float, int, int]:
    treatment = pd.to_numeric(treatment, errors="coerce").dropna()
    reference = pd.to_numeric(reference, errors="coerce").dropna()
    if len(treatment) == 0 or len(reference) == 0:
        return np.nan, np.nan, len(reference), len(treatment)
    mean_treat = treatment.mean()
    mean_ref = reference.mean()
    log2fc = np.log2((mean_treat + 1.0) / (mean_ref + 1.0))
    if len(treatment) < 2 or len(reference) < 2:
        pvalue = np.nan
    else:
        pvalue = ttest_ind(
            np.log2(treatment + 1.0),
            np.log2(reference + 1.0),
            equal_var=False,
            nan_policy="omit",
        ).pvalue
    return log2fc, pvalue, len(reference), len(treatment)


def add_status(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["status"] = np.select(
        [df["padj"].lt(0.05) & df["log2FC"].gt(0), df["padj"].lt(0.05) & df["log2FC"].lt(0)],
        ["Enriched", "Depleted"],
        default="No change",
    )
    df.loc[df["pvalue"].isna(), "status"] = "Not tested"
    return df


def main() -> None:
    raw = pd.read_excel(WORKBOOK, sheet_name=SHEET)
    feature_col = raw.columns[0]
    samples = list(raw.columns[1:])
    metadata = pd.DataFrame([parse_sample(s) for s in samples], columns=["species", "condition"], index=samples)

    # Keep the biological comparison scheme explicit: time-course controls use 0 h;
    # endpoint stress treatments use the late untreated control at the same endpoint.
    comparisons = [
        ("Bacillus", "mid", "0h", "B_mid", "0h baseline"),
        ("Bacillus", "late", "0h", "B_late", "0h baseline"),
        ("Bacillus", "Al", "late", "B_Al", "late untreated"),
        ("Bacillus", "Kana", "late", "B_K", "late untreated"),
        ("Bacillus", "Phage", "late", "B_P", "late untreated"),
        ("Rhodanobacter", "mid", "0h", "R_mid", "0h baseline"),
        ("Rhodanobacter", "late", "0h", "R_late", "0h baseline"),
        ("Rhodanobacter", "Al", "late", "R_Al", "late untreated"),
        ("Rhodanobacter", "Kana", "late", "R_K", "late untreated"),
        ("Rhodanobacter", "Phage", "late", "R_P", "late untreated"),
    ]

    records = []
    for species, treatment, reference, contrast, reference_label in comparisons:
        treat_samples = metadata.index[(metadata.species == species) & (metadata.condition == treatment)]
        ref_samples = metadata.index[(metadata.species == species) & (metadata.condition == reference)]
        for _, row in raw.iterrows():
            log2fc, pvalue, n_ref, n_treat = compare_feature(row[treat_samples], row[ref_samples])
            records.append(
                {
                    "feature": row[feature_col],
                    "species": species,
                    "treatment": treatment,
                    "reference": reference,
                    "reference_label": reference_label,
                    "contrast": contrast,
                    "n_reference": n_ref,
                    "n_treatment": n_treat,
                    "log2FC": log2fc,
                    "pvalue": pvalue,
                }
            )

    result = pd.DataFrame(records)
    result["padj"] = np.nan
    for contrast, idx in result.groupby("contrast").groups.items():
        p = result.loc[idx, "pvalue"].to_numpy(dtype=float)
        valid = np.isfinite(p)
        adj = np.full(len(p), np.nan)
        if valid.any():
            adj[valid] = multipletests(p[valid], method="fdr_bh")[1]
        result.loc[idx, "padj"] = adj
    result = add_status(result)
    result.to_csv(OUTDIR / "corrected_targeted_endpoint_reference_results.csv", index=False)

    summary = (
        result.groupby(["species", "contrast", "treatment", "reference", "reference_label", "status"], dropna=False)
        .size()
        .rename("n")
        .reset_index()
    )
    summary.to_csv(OUTDIR / "corrected_targeted_endpoint_reference_status_summary.csv", index=False)

    # Directly quantify how changing the reference changes the endpoint stress results.
    comparison_records = []
    for species, treatment, _, contrast, _ in comparisons:
        if treatment not in {"Al", "Kana", "Phage"}:
            continue
        treat_samples = metadata.index[(metadata.species == species) & (metadata.condition == treatment)]
        late_samples = metadata.index[(metadata.species == species) & (metadata.condition == "late")]
        zero_samples = metadata.index[(metadata.species == species) & (metadata.condition == "0h")]
        for _, row in raw.iterrows():
            late_fc, late_p, _, _ = compare_feature(row[treat_samples], row[late_samples])
            zero_fc, zero_p, _, _ = compare_feature(row[treat_samples], row[zero_samples])
            comparison_records.append(
                {
                    "feature": row[feature_col],
                    "species": species,
                    "treatment": treatment,
                    "contrast": contrast,
                    "log2FC_vs_late": late_fc,
                    "pvalue_vs_late": late_p,
                    "log2FC_vs_0h": zero_fc,
                    "pvalue_vs_0h": zero_p,
                    "delta_log2FC_late_minus_0h": late_fc - zero_fc,
                }
            )
    endpoint_compare = pd.DataFrame(comparison_records)
    for prefix in ["late", "0h"]:
        endpoint_compare[f"padj_vs_{prefix}"] = np.nan
        for contrast, idx in endpoint_compare.groupby("contrast").groups.items():
            p = endpoint_compare.loc[idx, f"pvalue_vs_{prefix}"].to_numpy(dtype=float)
            valid = np.isfinite(p)
            adj = np.full(len(p), np.nan)
            if valid.any():
                adj[valid] = multipletests(p[valid], method="fdr_bh")[1]
            endpoint_compare.loc[idx, f"padj_vs_{prefix}"] = adj
    endpoint_compare["status_vs_late"] = np.select(
        [endpoint_compare["padj_vs_late"].lt(0.05) & endpoint_compare["log2FC_vs_late"].gt(0),
         endpoint_compare["padj_vs_late"].lt(0.05) & endpoint_compare["log2FC_vs_late"].lt(0)],
        ["Enriched", "Depleted"], default="No change"
    )
    endpoint_compare["status_vs_0h"] = np.select(
        [endpoint_compare["padj_vs_0h"].lt(0.05) & endpoint_compare["log2FC_vs_0h"].gt(0),
         endpoint_compare["padj_vs_0h"].lt(0.05) & endpoint_compare["log2FC_vs_0h"].lt(0)],
        ["Enriched", "Depleted"], default="No change"
    )
    endpoint_compare.to_csv(OUTDIR / "endpoint_stress_reference_comparison.csv", index=False)

    status_wide = (
        endpoint_compare.groupby(["species", "treatment", "contrast", "status_vs_late", "status_vs_0h"])
        .size().rename("n_features").reset_index()
    )
    status_wide.to_csv(OUTDIR / "endpoint_status_transition_summary.csv", index=False)

    print("Endpoint-reference status counts")
    print(summary.to_string(index=False))
    print("\nEndpoint status transitions: late reference vs 0 h reference")
    print(status_wide.to_string(index=False))
    print("\nLargest absolute log2FC changes after switching reference")
    print(
        endpoint_compare.assign(abs_delta=lambda x: x.delta_log2FC_late_minus_0h.abs())
        .sort_values("abs_delta", ascending=False)
        .head(20)
        [["species", "treatment", "feature", "log2FC_vs_0h", "log2FC_vs_late", "delta_log2FC_late_minus_0h", "status_vs_0h", "status_vs_late"]]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
