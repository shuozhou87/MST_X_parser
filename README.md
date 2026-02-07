# Monolith X Raw Data App

A local Streamlit app for NanoTemper Monolith X raw data analysis.

## Features

- Drag and drop only the CSV/XLSX files you want to analyze.
- Auto-classify channels by filename suffix: `-650`, `-670`, `-R`.
- Read sample metadata from XLSX (capillary to ligand/analyte mapping).
- Manually merge experiments into custom merged sets.
- Manually select sample groups for comparison.
- Manually set MST laser-on window and compute single-point values.
- Show both raw time-signal series and window single-point comparison.
- Bar chart comparison with ANOVA and pairwise Welch t-tests.

## Quick start

1. Install dependencies:

```bash
python3 -m pip install -r requirements.txt
```

2. Run app:

```bash
streamlit run app.py
```

3. In sidebar:

- Drag files into upload box (`csv` + optional `Sample_Info.xlsx`).
- Choose `Sample info file` from uploaded xlsx files.
- Choose channels / experiments / groups.
- Optionally filter by ligand/analyte.
- Optionally enable manual experiment merge.
- In `Window Single-Point Comparison`, set laser-on window.

## Expected file naming

- Raw data: `Experiment-650.csv`, `Experiment-670.csv`, `Experiment-R.csv`
- Sample info: `.xlsx` containing `Capillary`, `Ligand`, `Analyte` columns (case-insensitive matching)
