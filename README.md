# go-enrichment-eggnog-pipeline

We provide a workflow to extract and process Gene Ontology (GO) annotations from eggNOG-mapper outputs, enabling GO enrichment analysis, redundancy reduction, and data visualization for non-model plant organisms.

## Easy web interface (Streamlit)

### 1. **Install required packages:**
```bash
pip install streamlit pandas goatools matplotlib scikit-learn numpy
```
*(Or use `pip install -r requirements.txt` for all dependencies)*

---

2. **Make sure your files are in your working directory:**
    - `streamlit_app.py` (the app script)
    - Your annotation file (eggNOG-mapper CSV)
    - Your gene list (CSV or TXT)
    - Your `go-basic.obo` file (download from [Gene Ontology](http://purl.obolibrary.org/obo/go/go-basic.obo))

3. **Launch the app:**
   ```
   streamlit run streamlit_app_1.py
   ```

4. **Follow the prompts in your browser!**
    - Upload your input files and run analysis.
    - Download results when finished.

5. **(Optional) Deploy online:**  
   Fork to [Streamlit Community Cloud](https://streamlit.io/cloud) for point-and-click access—no install needed for users.

---

**Tip:**  
If you use Anaconda, do this first:
```
conda create -n goenrich python=3.9
conda activate goenrich
```
