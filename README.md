## Easy web interface (Streamlit)

1. **Install required packages:**
   ```bash
   pip install streamlit pandas goatools matplotlib scikit-learn numpy
   ```
   *(Or use `pip install -r requirements.txt` for all dependencies)*

2. **Make sure your files are in your working directory:**
    - `streamlit_app.py` (the app script)
    - Your annotation file (eggNOG-mapper CSV file)
    - Your gene list (`.csv` or `.txt`)
    - The `go-basic.obo` ontology file ([download from Gene Ontology](http://purl.obolibrary.org/obo/go/go-basic.obo))

3. **Launch the app:**
   ```bash
   streamlit run streamlit_app.py
   ```

4. **Follow the interface in your browser:**
    - Upload your files and run the analysis.
    - Download output results.

5. **(Optional) Deploy online:**  
   Fork to [Streamlit Community Cloud](https://streamlit.io/cloud) for instant, web-based access—no install needed for users.

---

**Tip:**  
Using Anaconda? First run:
```bash
conda create -n goenrich python=3.9
conda activate goenrich
```

---

This enables anyone—biologists, non-coders, or students—to set up and use your GO enrichment tool quickly!
