# Sift - A drug discovery pipeline (v0.3.1)

## What this is

Sift is a personal project I'm building over summer 2026 to learn, hands-on, how AI and bioinformatics tools are actually used in early-stage drug discovery. I'm a biomedicine student, and going in I knew the theory around drug discovery only from lectures. This project is my way of testing the waters with the real, practical tools researchers use, using EGFR (a well-studied cancer drug target) as a concrete starting case study, entirely with free and public data and tools.

This isn't meant to be a novel scientific contribution. It's more a documented exploration of a real workflow, and a way to figure out where my own interests and skill gaps are before committing to a research direction.

## What it does

Given a target protein name, Sift chains together three real drug-discovery steps:

1. **Bioactivity screening** - trains a random forest classifier on real ChEMBL bioactivity data for that target and ranks candidate molecules by predicted probability of activity
2. **Structural docking** - docks top-ranked candidates against the target's real binding pocket using AutoDock Vina, if a structure has been prepared for it
3. **Literature mining** - pulls recent PubMed abstracts for that target and uses an LLM (Gemini) to extract known drugs and inhibitors, kept separate from the target's natural endogenous ligands, plus a plain-language disease context summary

Two supporting pieces tie it together:
- A **SMILES translator** (`name_to_smiles.py`) that resolves a plain drug or molecule name into its chemical structure via PubChem or OPSIN
- A **Streamlit UI** (`src/app.py`) that brings all of the above into one interface, including a bridge between the Sift and Literature tabs. A drug the literature step finds can be screened and docked inline with a single click.

It started as an EGFR-only proof of concept but has since generalized to arbitrary targets, validated so far on EGFR, BRAF, JAK2, mTOR, and GLUT4, as long as ChEMBL has bioactivity data for the target name given.

## Try it

```bash
streamlit run src/app.py
```

or from the command line:

```bash
python src/sift_target.py "TARGET_NAME"
python src/sift_target.py "TARGET_NAME" "molecule name"
python src/literature_mining.py "TARGET_NAME"
python src/name_to_smiles.py "molecule name"
```

See the Commands tab in the app for the full, current list of CLI entry points.

## Roadmap

- [x] v0.1 - Bioactivity screening model, validated against approved EGFR drugs
- [x] v0.2 - Structural docking validation (AutoDock Vina against EGFR's binding pocket)
- [x] v0.2.1 - Generalized screening and docking to arbitrary targets, plus a SMILES translator and repo reorganization
- [x] v0.2.2 - Streamlit UI (Sift, SMILES translator, Commands tabs)
- [x] v0.3 - LLM-based literature mining (PubMed context and known-ligand extraction, via Gemini)
- [x] v0.3.1 - Literature mining wired into the UI (Literature tab, Sift/Literature bridge, disk caching)
- [ ] v1.0 - Full write-up, polish, and further UI work

## Reflections

**v0.1 - Bioactivity screening**
Watching the model rediscover real drugs (lapatinib, erlotinib) just from fingerprint patterns made "structure-activity relationship" click for the first time. Still an open question how much of that success is real chemistry versus the model picking up on which molecules ChEMBL just happens to have the most data for.

**v0.2 - Docking, and generalizing beyond EGFR**
Docking made the screening scores feel less abstract: from "this looks similar to known drugs" to "here's how it would actually sit in the pocket." Real data is messy, too. Salt forms like lapatinib ditosylate needed cleanup before docking would even run.

Generalizing past EGFR (`train_target.py`, `prepare_target_structure.py`, `sift_target.py`) was the point this stopped feeling like a one-off script. I validated it on BRAF (screening ROC-AUC 0.827, docking redocked vemurafenib at -8.705 kcal/mol), then stress-tested it on JAK2 and GLUT4, including a tiny 129-molecule dataset, just to make sure it wasn't only working because EGFR and BRAF happen to have tons of cancer data.

Reorganizing the repo into `src/`, `data/`, `notebooks/` taught me a real lesson about relative paths: they only worked when scripts ran from inside `src/`, so I switched to `__file__`-based absolute paths.

Then came the Streamlit UI. Mocking up the design first was much faster than iterating through Streamlit's reload cycle, and reviewing it caught a real bug: the candidates view was never actually showing docking scores, even when a structure existed for that target.

Open question: how little ChEMBL data can a target have before I should stop trusting the model's confidence.

**v0.3 - Literature mining**
Picked Gemini's free API to keep the whole pipeline free. Hit friction fast: the model name from the docs, `gemini-2.5-flash`, was already deprecated, and I had to switch to `gemini-3-flash-preview`. Then had a genuinely annoying stretch getting a Python file with nested quotes through shell heredocs without it silently corrupting mid-paste, more than once. Lesson: past a certain quoting depth, it's faster to just rewrite the whole file than to patch it.

Also hit the free tier's real limits: a 503 on my first real run, then a full daily quota cutoff mid-session. Free seems to mean best-effort and rate-limited, not unreliable, but you do have to plan around it.

The best design call was splitting ligands into drug versus endogenous categories. Without it, EGF and epiregulin were mixed in with actual drugs like cetuximab, which would have been misleading.

The Sift-Literature bridge, where you can screen and dock a literature-found drug in one click, ended up being the most satisfying feature so far. It's the same validation logic I've used the whole project, except now the "should be good" list comes from real papers instead of my own picks.

Open question: would a paid-tier model meaningfully improve extraction quality, and is that worth paying for.
