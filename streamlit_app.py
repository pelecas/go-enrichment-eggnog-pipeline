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
import filetype
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

def verify_obo_file(filebytes):
    """Verify if the uploaded file is a valid OBO file."""
    try:
        kind = filetype.guess(filebytes)
        if kind:
            if kind.extension == "obo":
                return True
            else:
                st.warning(f"Uploaded file is detected as {kind.extension} but expected OBO format.")
                return False
        else:
            st.warning("Could not determine file type. Ensure you upload a valid OBO file.")
            return False
    except Exception as e:
        st.error(f"Error verifying OBO file: {e}")
        return False

anno_df = robust_annotation_upload("1. Annotation file (CSV)")
study_df = robust_gene_list_upload("2. Study gene list (CSV)")
go_obo_file = st.file_uploader("3. GO OBO file (obo)", type=['obo'])

namespace_to_run = st.selectbox("Choose GO namespace:", ["BP", "MF", "CC"])
alpha = st.number_input("Significance threshold (FDR)", min_value=0.001, max_value=0.5, value=0.05)
run_btn = st.button("Run GO Enrichment")

if go_obo_file and run_btn:
    filebytes = go_obo_file.read()
    if not verify_obo_file(filebytes):
        st.stop()
    # Save the file locally for the GODag parser
    with open("go-basic.obo", "wb") as f:
        f.write(filebytes)
else:
    st.stop()

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

# remaining code logic is similar to the previous module passed.