#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 20:34:27 2025

@author: Eliel
"""

import streamlit as st
import pandas as pd
import numpy as np
import chardet
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
import matplotlib.pyplot as plt

st.title("GO Enrichment & Redundancy Reduction Tool")

st.markdown("""
Upload your annotation file (eggNOG output), your gene list, and a GO OBO file.<br>
Dotplots display the top GO enrichment hits:<br>
- **Dot size** = Number of genes with that GO term
- **Dot color** = Significance (-log10 p-value) gradient<br>
- **Bubble legend** = Gene count<br>
- **Colorbar** = P-value scale<br>
""", unsafe_allow_html=True)

def robust_annotation_upload(label, file_types=['csv', 'txt', 'tsv']):
    upload_file = st.file_uploader(label, type=file_types)
    if upload_file:
        raw_bytes = upload_file.read()
        upload_file.seek(0)
        encoding_guess = chardet.detect(raw_bytes)['encoding']
        st.info(f"Detected encoding: {encoding_guess}")
        try:
            df = pd.read_csv(
                upload_file,
                header=0,
                encoding=encoding_guess,
                engine="python",
                sep=None,
                on_bad_lines="skip"
            )
            st.success(f"Successfully loaded file. Columns detected: {list(df.columns)}")
        except Exception as e:
            st.error(f"Could not read file: {e}")
            return None

        st.write("Preview:", df.head())
        col_gene = st.selectbox(f"Pick the gene ID column for [{label}]:", df.columns, key=label+"-col-gene")
        col_go = st.selectbox(f"Pick the GO/annotation column for [{label}]:", df.columns, key=label+"-col-go")
        df_sub = df[[col_gene, col_go]].copy()
        df_sub.columns = ['gene_id', 'Gos']
        st.write("Mapped Annotation DataFrame:", df_sub.head())
        return df_sub
    return None

def robust_gene_list_upload(label, file_types=['csv', 'txt', 'tsv']):
    upload_file = st.file_uploader(label, type=file_types)
    if upload_file:
        raw_bytes = upload_file.read()
        upload_file.seek(0)
        encoding_guess = chardet.detect(raw_bytes)['encoding']
        st.info(f"Detected encoding: {encoding_guess}")
        try:
            df = pd.read_csv(
                upload_file,
                header=0,
                encoding=encoding_guess,
                engine="python",
                sep=None,
                on_bad_lines="skip"
            )
            st.success(f"Successfully loaded file. Columns detected: {list(df.columns)}")
        except Exception as e:
            st.error(f"Could not read file: {e}")
            return None

        st.write("Preview:", df.head())
        col_gene = st.selectbox(f"Pick the gene ID column for [{label}]:", df.columns, key=label+"-col-gene")
        df_sub = df[[col_gene]].copy()
        df_sub.columns = ['gene_id']
        st.write("Mapped Gene List DataFrame:", df_sub.head())
        return df_sub
    return None

anno_df = robust_annotation_upload("1. Annotation file (CSV)")
study_df = robust_gene_list_upload("2. Study gene list (CSV)")
go_obo_file = st.file_uploader("3. GO OBO file (obo)", type=['obo'])

namespace_to_run = st.selectbox("Choose GO namespace:", ["BP", "MF", "CC"])
alpha = st.number_input("Significance threshold (FDR)", min_value=0.001, max_value=0.5, value=0.05)
run_btn = st.button("Run GO Enrichment")

def split_namespace(go_df, obodag, namespace):
    def ns_for_go(go_id):
        if go_id in obodag:
            term = obodag[go_id]
            if hasattr(term, 'namespace'):
                if term.namespace.startswith('biological_process') and namespace == "BP":
                    return "BP"
                elif term.namespace.startswith('molecular_function') and namespace == "MF":
                    return "MF"
                elif term.namespace.startswith('cellular_component') and namespace == "CC":
                    return "CC"
        return None
    return go_df[go_df["go_id"].map(ns_for_go) == namespace].copy()

def semantic_reduction(df_out, obodag, gene2go):
    pval_col = None
    for col in ['p_uncorrected', 'p_fdr_bh', 'pval_uncorrected', 'p_value']:
        if col in df_out.columns:
            pval_col = col
            break
    if pval_col is None or df_out.empty:
        st.warning("No p-value column found for reduction. Returning all terms.")
        return df_out
    sig_terms = df_out[df_out[pval_col] < 0.05]['GO_ID'].tolist()
    terms_for_sim = [go for go in sig_terms if go in obodag]
    n_terms = len(terms_for_sim)
    if n_terms == 0:
        return df_out
    try:
        from goatools.semantic import TermCounts, resnik_sim
        from sklearn.cluster import AgglomerativeClustering
        termcounts = TermCounts(obodag, gene2go)
        dist_matrix = np.ones((n_terms, n_terms)) * 10
        for i, go1 in enumerate(terms_for_sim):
            for j, go2 in enumerate(terms_for_sim):
                if i < j and obodag[go1].namespace == obodag[go2].namespace:
                    sim = resnik_sim(go1, go2, obodag, termcounts) or 0
                    dist = 10 - sim
                    dist_matrix[i, j] = dist_matrix[j, i] = dist
        np.fill_diagonal(dist_matrix, 0)
        if n_terms > 1:
            clustering = AgglomerativeClustering(
                n_clusters=None, distance_threshold=6.0, metric='precomputed', linkage='average')
            labels = clustering.fit_predict(dist_matrix)
            rep_terms = []
            for cid in np.unique(labels):
                members = [terms_for_sim[i] for i in range(n_terms) if labels[i] == cid]
                best_idx = df_out[df_out['GO_ID'].isin(members)][pval_col].idxmin()
                best = df_out.loc[best_idx, 'GO_ID']
                rep_terms.append(best)
        else:
            rep_terms = terms_for_sim
        return df_out[df_out['GO_ID'].isin(rep_terms)]
    except Exception as e:
        st.warning(f"Reduction failed: {str(e)}")
        return df_out

if (anno_df is not None) and (study_df is not None) and go_obo_file and run_btn:
    with st.spinner("Processing..."):
        df_go = anno_df[anno_df['Gos'].notna()]
        go_map = []
        for pid, goterms in df_go.values:
            for go in str(goterms).split(','):
                go = go.strip()
                if go and go.startswith("GO:"):
                    go_map.append((pid, go))
        go_df = pd.DataFrame(go_map, columns=['gene_id', 'go_id'])
        go_df = go_df[go_df['go_id'] != 'go_id']
        go_df = go_df[go_df['go_id'].str.startswith("GO:")]

        st.write("Loaded gene-GO mapping:")
        st.write(go_df.head())
        st.write(f"Unique GO terms: {len(go_df['go_id'].unique())}")

        study_genes = study_df['gene_id'].astype(str).tolist()
        study_genes = [g for g in study_genes if g in go_df['gene_id'].unique()]
        st.write(f"Total genes in study set after mapping to annotation: {len(study_genes)}")

        go_obo_bytes = go_obo_file.read()
        with open("go-basic.obo", "wb") as f:
            f.write(go_obo_bytes)
        obodag = GODag("go-basic.obo")
        st.write(f"OBO has {len(obodag)} terms.")

        st.write(f"Filtering annotation to namespace: {namespace_to_run}")
        go_ns_df = split_namespace(go_df, obodag, namespace_to_run)
        st.write(f"{namespace_to_run}: {len(go_ns_df)} cleaned gene-GO pairs")

        gene2go = go_ns_df.groupby('gene_id')['go_id'].apply(list).to_dict()
        background_genes = list(gene2go.keys())
        study_ns_genes = [g for g in study_genes if g in background_genes]
        st.write("Background genes:", len(background_genes), "Study genes (namespace):", len(study_ns_genes))

        study_assigned_gos = []
        for g in study_ns_genes:
            terms = gene2go.get(g, [])
            study_assigned_gos.extend(terms)
        st.write(f"Total GO assignments for study genes in namespace: {len(study_assigned_gos)}")
        st.write(f"Unique GO terms assigned to study genes: {len(set(study_assigned_gos))}")
        st.write("Sample assigned GO terms:", list(set(study_assigned_gos))[:20])

        enrichment = GOEnrichmentStudy(
            background_genes, gene2go, obodag, alpha=alpha, methods=['fdr_bh'], propagate_counts=False
        )
        results = enrichment.run_study(study_ns_genes)
        st.write("Number of GOATOOLS results:", len(results))

        df_out = pd.DataFrame([
            {
                "GO_ID": getattr(r, "GO", None),
                "Name": getattr(r, "name", None),
                "Namespace": getattr(r, "NS", None),
                "p_uncorrected": getattr(r, "p_uncorrected", None),
                "p_fdr_bh": getattr(r, "p_fdr_bh", None),
                "Genes": ";".join(getattr(r, "study_items", []))
            }
            for r in results
        ], columns=["GO_ID", "Name", "Namespace", "p_uncorrected", "p_fdr_bh", "Genes"])

        st.write(df_out.head())
        st.write(df_out.describe())

        st.download_button("Download raw enrichment results (CSV)", df_out.to_csv(index=False), "go_enrichment_results.csv", "text/csv")

        reduced_df = semantic_reduction(df_out, obodag, gene2go)
        st.download_button("Download reduced GO terms (CSV)", reduced_df.to_csv(index=False), "go_enrichment_results_reduced.csv", "text/csv")

        # --- Dual Dotplot visualization (uncorrected AND FDR, matching bubble legend, wider figures) ---
        N = min(20, reduced_df.shape[0])
        if N > 0:
            df_plot = reduced_df.nsmallest(N, 'p_fdr_bh').copy()
            df_plot['-log10_p_uncorrected'] = -np.log10(df_plot['p_uncorrected'].replace(0, 1e-15))
            df_plot['-log10_p_fdr_bh'] = -np.log10(df_plot['p_fdr_bh'].replace(0, 1e-15))
            df_plot['gene_count'] = df_plot['Genes'].apply(lambda s: len(str(s).split(';')))
            sizes = df_plot['gene_count'] * 15

            min_size = df_plot['gene_count'].min()
            max_size = df_plot['gene_count'].max()
            legend_example_counts = np.linspace(min_size, max_size, 4).astype(int)
            legend_sizes = legend_example_counts * 15

            col1, col2 = st.columns(2)

            # Plot for p_uncorrected
            with col1:
                fig1, ax1 = plt.subplots(figsize=(9, 1 + 0.7*N))
                x1 = df_plot['-log10_p_uncorrected']
                y = np.arange(len(df_plot))

                scatter1 = ax1.scatter(
                    x1, y, s=sizes, c=x1, cmap='viridis', alpha=0.85, edgecolor='black'
                )
                ax1.set_yticks(y)
                ax1.set_yticklabels(df_plot['Name'])
                ax1.set_xlabel(r'$-\log_{10}$ p_uncorrected')
                ax1.set_title('Raw p-value Dotplot')
                ax1.invert_yaxis()
                plt.tight_layout()

                # Bubble legend (matches dot size)
                bubble_handles = [
                    plt.scatter([], [], s=s, c='gray', edgecolors='black', label=f'{int(count)} gene{"s" if int(count)>1 else ""}')
                    for count, s in zip(legend_example_counts, legend_sizes)
                ]
                bubble_legend = ax1.legend(
                    handles=bubble_handles,
                    title="Gene Count",
                    handletextpad=2,
                    loc='upper left',
                    bbox_to_anchor=(1.01, 0.95),
                    frameon=True
                )
                ax1.add_artist(bubble_legend)

                # Colorbar, thin and shrunk
                norm1 = plt.Normalize(x1.min(), x1.max())
                sm1 = plt.cm.ScalarMappable(cmap="viridis", norm=norm1)
                sm1.set_array([])
                cbar1 = plt.colorbar(sm1, ax=ax1, pad=0.03, aspect=25, shrink=0.4)
                cbar1.set_label(r'$-\log_{10}$ p_uncorrected', fontsize=12)

                st.pyplot(fig1)
                # Export if needed:
                # fig1.savefig("go_dotplot_p_uncorrected.png", dpi=300, bbox_inches='tight')

            # Plot for p_fdr_bh
            with col2:
                fig2, ax2 = plt.subplots(figsize=(9, 1 + 0.7*N))
                x2 = df_plot['-log10_p_fdr_bh']

                scatter2 = ax2.scatter(
                    x2, y, s=sizes, c=x2, cmap='viridis', alpha=0.85, edgecolor='black'
                )
                ax2.set_yticks(y)
                ax2.set_yticklabels(df_plot['Name'])
                ax2.set_xlabel(r'$-\log_{10}$ FDR')
                ax2.set_title('FDR-corrected p-value Dotplot')
                ax2.invert_yaxis()
                plt.tight_layout()

                # Bubble legend (matches dot size)
                bubble_handles2 = [
                    plt.scatter([], [], s=s, c='gray', edgecolors='black', label=f'{int(count)} gene{"s" if int(count)>1 else ""}')
                    for count, s in zip(legend_example_counts, legend_sizes)
                ]
                bubble_legend2 = ax2.legend(
                    handles=bubble_handles2,
                    title="Gene Count",
                    handletextpad=2,
                    loc='upper left',
                    bbox_to_anchor=(1.01, 0.95),
                    frameon=True
                )
                ax2.add_artist(bubble_legend2)

                # Colorbar, thin and shrunk
                norm2 = plt.Normalize(x2.min(), x2.max())
                sm2 = plt.cm.ScalarMappable(cmap="viridis", norm=norm2)
                sm2.set_array([])
                cbar2 = plt.colorbar(sm2, ax=ax2, pad=0.03, aspect=25, shrink=0.4)
                cbar2.set_label(r'$-\log_{10}$ FDR', fontsize=12)

                st.pyplot(fig2)
                # Export if needed:
                # fig2.savefig("go_dotplot_p_fdr_bh.png", dpi=300, bbox_inches='tight')

        else:
            st.warning("No significant/reduced terms found for plot.")