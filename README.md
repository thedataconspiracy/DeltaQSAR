# Δ-QSAR: A Knowledge Discovery Suite for Automatic SAR/QSAR Rules Induction

![Python](https://img.shields.io/badge/Python-100%25-green)

Δ-QSAR (Delta-QSAR) is a knowledge discovery suite designed to automate the induction of Structure-Activity Relationship (SAR) and Quantitative Structure-Activity Relationship (QSAR) rules. This tool empowers researchers, chemists, and data scientists to derive meaningful insights from chemical datasets, enabling advancements in drug discovery, material science, and other related fields.

Version: 0.2
- SARpy v2.0 (SAR)
- QSARpy v2.1 (QSAR)

Copyright (C) 2021-2025 Thomas Ferrari

---

## Overview

This repository provides the code, tutorials and an environment specification to reproduce the development runtime. 

## Installation

1. Clone the repository:
```bash
git clone https://github.com/thedataconspiracy/DeltaQSAR.git
cd DeltaQSAR
```

2. Recreate the Conda environment (environment.yml is included in the repo):
```bash
# using conda
conda env create -f environment.yml

# then activate the environment (environment name is defined inside environment.yml; default: DeltaQSAR)
conda activate DeltaQSAR
```

## Usage

Start the graphical user interface:
```bash
# from the repository root, with the conda environment activated
python GUI.py
```

When running on Windows, prebuilt executable binaries are available (see the "Windows binaries" link below)

## Documentation (local files)

Two HTML tutorials are included in this repository. Open them in your browser for step-by-step guidance:

- [SARpy Quickstart](https://thedataconspiracy.github.io/DeltaQSAR/SARpy2_Quickstart.html)
- [QSARpy Manual](https://thedataconspiracy.github.io/DeltaQSAR/QSARpy2_Manual.html)

## Windows executable binaries

Prebuilt Windows executable binaries are available for download here:
https://drive.google.com/file/d/1P3iWHvTH-2kqHzA0jTmFnQXMWHzGSj3A/view

