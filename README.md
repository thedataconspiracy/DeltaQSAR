# Δ-QSAR: A Knowledge Discovery Suite for Automatic SAR/QSAR Rules Induction

![Python](https://img.shields.io/badge/Python-100%25-green)

Δ-QSAR (Delta-QSAR) is a knowledge discovery suite designed to automate the induction of Structure-Activity Relationship (SAR) and Quantitative Structure-Activity Relationship (QSAR) rules. This tool empowers researchers, chemists, and data scientists to derive meaningful insights from chemical datasets, enabling advancements in drug discovery, material science, and other related fields.

Version: 0.2
- SARpy v2.0
- QSARpy v2.1

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

When running on Windows, prebuilt executable binaries are available (see the "Windows binaries" link below) — download and run the appropriate installer or executable.

## Documentation (local files)

Two HTML tutorials are included in this repository. Open them in your browser for step-by-step guidance:

- [SARpy quickstart](https://thedataconspiracy.github.io/DeltaQSAR/SARpy2_Quickstart.html)
- [QSAR manual](https://thedataconspiracy.github.io/DeltaQSAR/QSARpy2_Manual.html)

(You can also view these files directly on GitHub in the repository root.)

## Windows executable binaries

Prebuilt Windows executable binaries are available for download here:
https://drive.google.com/file/d/1P3iWHvTH-2kqHzA0jTmFnQXMWHzGSj3A/view

