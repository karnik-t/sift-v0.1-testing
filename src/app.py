import streamlit as st
import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sift_target import ensure_model, smiles_to_fp, safe, DATA_DIR
from name_to_smiles import name_to_smiles
from pipeline import dock_smiles

st.set_page_config(page_title="SIFT", page_icon="🧬", layout="centered")

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Instrument+Serif:ital@0;1&family=Manrope:wght@400;500;600&display=swap');

.stApp { background-color: #F2F0EF; font-family: 'Manrope', sans-serif; }
h2, h3, p, label, .stMarkdown { color: #606C38 !important; font-family: 'Manrope', sans-serif; }

.stTextInput input {
    background-color: #ffffff;
    color: #606C38;
    border: 1px solid #dcdad7;
    border-radius: 8px;
    font-family: 'Manrope', sans-serif;
}
.stButton>button {
    background-color: #BC6C25;
    color: #FFFFFF !important;
    font-weight: 600;
    border-radius: 8px;
    border: none;
    padding: 0.5em 1.5em;
    font-family: 'Manrope', sans-serif;
}
.stButton>button:hover { background-color: #a35a1f; color: #FFFFFF !important; }
.stButton>button p { color: #FFFFFF !important; }

.stRadio label p { color: #606C38 !important; font-family: 'Manrope', sans-serif; }
.stTabs [data-baseweb="tab"] { color: #a3a29e; font-family: 'Manrope', sans-serif; }
.stTabs [aria-selected="true"] {
    color: #606C38 !important;
    border-bottom-color: #BC6C25 !important;
}
[data-testid="stMetricValue"] { color: #BC6C25; font-family: 'Manrope', sans-serif; }
[data-testid="stMetricLabel"] { color: #606C38; }
[data-testid="stVerticalBlockBorderWrapper"] > div {
    background-color: #ffffff;
    border-radius: 10px;
    border: none;
}
pre, code {
    background-color: #ffffff !important;
    color: #BC6C25 !important;
    border: 1px solid #dcdad7 !important;
}

/* Slider track + thumb */
[data-testid="stSlider"] [role="slider"] {
    background-color: #BC6C25 !important;
    border-color: #BC6C25 !important;
}
[data-testid="stSlider"] div[data-baseweb="slider"] > div > div {
    background-color: #BC6C25 !important;
}
</style>
""", unsafe_allow_html=True)

st.markdown(
    "<div style='display:flex; align-items:baseline; gap:12px; margin-bottom:8px;'>"
    "<span style=\"font-family:'Instrument Serif',serif; font-size:44px; letter-spacing:1px; color:#606C38;\">SIFT</span>"
    "<span style='font-size:13px; color:#8a9264;'>An AI-assisted drug discovery pipeline</span>"
    "</div>",
    unsafe_allow_html=True
)

tab_app, tab_translate, tab_commands = st.tabs(["Sift", "SMILES translator", "Commands"])

with tab_app:
    with st.container(border=True):
        st.markdown("**Target protein**")
        target_name = st.text_input("Target protein", label_visibility="collapsed",
                                     placeholder="EGFR, BRAF, JAK2, mTOR...")

        st.markdown("**What do you want to do?**")
        mode = st.radio("mode", ["Show top candidates", "Screen a specific molecule"],
                         label_visibility="collapsed")

        if mode == "Show top candidates":
            top_n = st.slider("How many candidates?", 5, 20, 10)
            run_screen = st.button("Run screening")
        else:
            molecule_name = st.text_input("Molecule name", placeholder="common or IUPAC name")
            run_screen = st.button("Screen this molecule")

    if mode == "Show top candidates" and run_screen and target_name:
        with st.spinner(f"Screening candidates against {target_name}..."):
            model = ensure_model(target_name)
            df = pd.read_csv(os.path.join(DATA_DIR, "candidates.csv"))
            df["fingerprint"] = df["smiles"].apply(smiles_to_fp)
            df = df.dropna(subset=["fingerprint"])
            X = np.stack(df["fingerprint"].values)
            df["predicted_prob_active"] = model.predict_proba(X)[:, 1]
            top = df.sort_values("predicted_prob_active", ascending=False).head(top_n)

        with st.container(border=True):
            st.markdown(f"**Top {top_n} candidates for {target_name}**")
            st.caption("Score = predicted probability of activity, 0 to 1 (model confidence, not a physical unit)")
            for _, row in top.iterrows():
                st.markdown(
                    f"<div style='display:flex; justify-content:space-between; padding:8px 0; "
                    f"border-bottom:1px solid #f0efed; font-size:14px;'>"
                    f"<span style='color:#606C38;'>{row['name'].title()}</span>"
                    f"<span style='font-weight:600; color:#BC6C25;'>{row['predicted_prob_active']:.3f}</span>"
                    f"</div>", unsafe_allow_html=True
                )

    elif mode == "Screen a specific molecule" and run_screen and target_name and molecule_name:
        with st.spinner("Resolving molecule..."):
            smiles, source = name_to_smiles(molecule_name)
        if smiles is None:
            st.error(f"Could not resolve '{molecule_name}' to a structure.")
        else:
            st.info(f"Resolved via {source}: `{smiles}`")
            with st.spinner(f"Screening against {target_name}..."):
                model = ensure_model(target_name)
                fp = smiles_to_fp(smiles)
                prob = model.predict_proba(fp.reshape(1, -1))[0][1]
            st.metric("Predicted probability of activity (0-1 scale)", f"{prob:.3f}")

            receptor_path = os.path.join(DATA_DIR, "structures", f"{safe(target_name)}_receptor.pdbqt")
            pocket_path = os.path.join(DATA_DIR, "structures", f"{safe(target_name)}_pocket.txt")
            if os.path.exists(receptor_path) and os.path.exists(pocket_path):
                with st.spinner("Docking..."):
                    score = dock_smiles(smiles, safe(molecule_name)[:30])
                st.metric("Predicted binding affinity", f"{score} kcal/mol")
                st.caption("More negative = stronger predicted binding")
            else:
                st.warning(f"No prepared docking structure for {target_name} yet.")

with tab_translate:
    st.markdown("**Translate a molecule name to its chemical structure (SMILES)**")
    st.caption("Works with common drug names or IUPAC systematic names")
    lookup_name = st.text_input("Molecule name", placeholder="e.g. aspirin, or 4-methylpent-1-ene")
    if st.button("Translate"):
        with st.spinner("Looking up..."):
            smiles, source = name_to_smiles(lookup_name)
        if smiles is None:
            st.error(f"Could not resolve '{lookup_name}' via PubChem or OPSIN.")
        else:
            st.success(f"Resolved via {source}")
            st.code(smiles, language=None)

with tab_commands:
    st.markdown("**Screen a target for top candidates:**")
    st.code('python src/sift_target.py "TARGET_NAME"', language="bash")

    st.markdown("**Screen (and dock, if prepared) a specific molecule against a target:**")
    st.code('python src/sift_target.py "TARGET_NAME" "molecule name"', language="bash")

    st.markdown("**Prepare a target for docking (needs a PDB ID + ligand code from rcsb.org):**")
    st.code('python src/prepare_target_structure.py "TARGET_NAME" "PDB_ID" "LIGAND_CODE"', language="bash")

    st.markdown("**Translate a molecule name to its chemical structure (SMILES):**")
    st.code('python src/name_to_smiles.py "molecule name"', language="bash")

    st.markdown("**Train a screening model for a target directly:**")
    st.code('python src/train_target.py "TARGET_NAME"', language="bash")

    st.markdown("**Launch this app:**")
    st.code('streamlit run src/app.py', language="bash")
