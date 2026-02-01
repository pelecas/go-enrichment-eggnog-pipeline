#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 20:34:27 2025

@author: Eliel
"""
import pandas as pd
import numpy as np
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
import matplotlib.pyplot as plt
import chardet
import streamlit as st
import tempfile
import os

# ===============================
# Streamlit App Configuration
# ===============================
st.set_page_config(layout="wide")
st.title("GO Enrichment Analysis and Visualization App")

st.write("""
Upload:
- Annotation file (CSV / TSV)
- Study gene list (CSV)
- GO OBO file

The app performs GO enrichment analysis and generates a
publication-quality dot plot with SVG export.
""")

# ===============================
# Robust file upload
# ===============================
def robust_file_upload(label, file_types):
    upload_file = st.file_uploader(label, type=file_types)
    if upload_file:
        try:
            raw = upload_file.read()
            enc = chardet.detect(raw).get("encoding", "utf-8")
            upload_file.seek(0)
            df = pd.read_csv(upload_file, encoding=enc)
            st.success(f"Loaded: {upload_file.name}")
            st.dataframe(df.head())
            return df
        except Exception as e:
            st.error(f"File error: {e}")
            return None
    return None

# ===============================
# Upload data
# ===============================
anno_file = robust_file_upload("Upload Annotation File", ["csv", "tsv", "txt"])
if anno_file is None:
    st.stop()

study_file = robust_file_upload("Upload Study Gene List", ["csv", "tsv", "txt"])
if study_file is None:
    st.stop()

go_obo_file = st.file_uploader("Upload GO OBO File (.obo)", type=["obo"])
if go_obo_file is None:
    st.stop()

namespace_to_run = st.selectbox("GO Namespace", ["BP", "MF", "CC"])
alpha = st.slider("FDR threshold", 0.01, 0.1, 0.05, 0.01)

# ===============================
# Data processing
# ===============================
try:
    anno_df = anno_file.iloc[:, [0, 9]].copy()
    anno_df.columns = ["gene_id", "Gos"]
    anno_df.replace(["-", ""], pd.NA, inplace=True)
    anno_df = anno_df.dropna()

    go_pairs = []
    for g, gos in anno_df.values:
        for go in str(gos).split(","):
            if go.startswith("GO:"):
                go_pairs.append((g, go.strip()))

    go_df = pd.DataFrame(go_pairs, columns=["gene_id", "go_id"])
    st.write(f"Gene–GO mappings: {len(go_df)}")

    study_genes = set(study_file.iloc[:, 0].astype(str))
    study_genes = list(study_genes & set(go_df["gene_id"]))
    st.write(f"Study genes after filtering: {len(study_genes)}")

    # ---- Save OBO file to disk (required by GOATOOLS)
    with tempfile.NamedTemporaryFile(delete=False, suffix=".obo") as tmp:
        tmp.write(go_obo_file.getvalue())
        obo_path = tmp.name

    obodag = GODag(obo_path)
    st.write(f"Loaded GO DAG terms: {len(obodag)}")

    def split_namespace(go_df, obodag, ns):
        def match(go):
            if go in obodag:
                n = obodag[go].namespace
                return (
                    (n == "biological_process" and ns == "BP") or
                    (n == "molecular_function" and ns == "MF") or
                    (n == "cellular_component" and ns == "CC")
                )
            return False

        return go_df[go_df["go_id"].map(match)]

    filtered_go_df = split_namespace(go_df, obodag, namespace_to_run)

    gene2go = (
        filtered_go_df
        .groupby("gene_id")["go_id"]
        .apply(set)   # MUST be set
        .to_dict()
    )

    background_genes = list(gene2go.keys())
    study_ns_genes = [g for g in study_genes if g in background_genes]

    st.write(
        f"{namespace_to_run}: "
        f"Background = {len(background_genes)}, "
        f"Study = {len(study_ns_genes)}"
    )

except Exception as e:
    st.error(f"Processing error: {e}")
    st.stop()

# ===============================
# Enrichment analysis
# ===============================
try:
    goea = GOEnrichmentStudy(
        background_genes,
        gene2go,
        obodag,
        alpha=alpha,
        methods=["fdr_bh"],
        propagate_counts=True
    )

    results = goea.run_study(study_ns_genes)

    df_out = pd.DataFrame([
        {
            "GO_ID": r.GO,
            "Name": r.name,
            "Namespace": r.NS,
            "p_fdr_bh": r.p_fdr_bh,
            "Genes": ";".join(r.study_items)
        }
        for r in results if r.p_fdr_bh is not None
    ])

    df_out["study_count"] = df_out["Genes"].str.count(";") + 1
    df_out = df_out.sort_values("p_fdr_bh")

    st.subheader("Enrichment Results")
    st.dataframe(df_out.head(20))

    st.download_button(
        "Download enrichment table (CSV)",
        df_out.to_csv(index=False),
        "go_enrichment_results.csv"
    )

except Exception as e:
    st.error(f"Enrichment error: {e}")
    st.stop()

# ===============================
# Dot plot with legends + SVG
# ===============================
def create_bubble_plot(df, namespace):
    df = df[df["p_fdr_bh"] > 0].copy()
    df["-log10_fdr"] = -np.log10(df["p_fdr_bh"])
    df = df.sort_values("-log10_fdr", ascending=False).head(20)

    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(10, len(df) * 0.5))

    sc = ax.scatter(
        df["-log10_fdr"],
        y,
        s=df["study_count"] * 15,
        c=df["-log10_fdr"],
        cmap="viridis",
        edgecolors="black",
        alpha=0.85
    )

    ax.set_yticks(y)
    ax.set_yticklabels(df["Name"], fontsize=10)
    ax.set_xlabel(r"-log$_{10}$(FDR)")
    ax.set_title(f"Top 20 {namespace} GO Terms")
    ax.invert_yaxis()

    cbar = plt.colorbar(sc, ax=ax, shrink=0.45, pad=0.04, aspect=20)
    cbar.set_label(r"-log$_{10}$ FDR")

    # Dot size legend
    sizes = [5, 10, 20, 40]
    handles = [
        plt.scatter([], [], s=s * 15, color="gray", edgecolors="black", alpha=0.6)
        for s in sizes
    ]

    ax.legend(
        handles,
        [f"{s} genes" for s in sizes],
        title="Gene count",
        loc="center left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=False
    )

    plt.tight_layout()
    st.pyplot(fig)

    # Save SVG
    with tempfile.NamedTemporaryFile(delete=False, suffix=".svg") as tmp:
        fig.savefig(tmp.name, format="svg")
        with open(tmp.name, "rb") as f:
            st.download_button(
                "Download dot plot (SVG)",
                f,
                f"GO_{namespace}_dotplot.svg",
                "image/svg+xml"
            )
        os.remove(tmp.name)

if not df_out.empty:
    create_bubble_plot(df_out, namespace_to_run)

# ===============================
# Cleanup
# ===============================
try:
    os.remove(obo_path)
except Exception:
    pass
