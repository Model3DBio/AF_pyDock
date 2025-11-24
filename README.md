# ⚠️ This repository is currently undergoing refactoring. Please wait a few days before using it.
# AF_pyDock
# Analysis and Calculation of Model Confidence and pyDock Energy

This repository contains a Jupyter Notebook (`analysis_protocol.ipynb`) that provides a detailed protocol for the analysis and calculation of Model Confidence and pyDock Energy functions. The notebook is intended for use in structural biology and computational docking studies, where these functions are critical for evaluating the reliability of molecular models.

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Notebook Overview](#notebook-overview)
- [References](#references)
- [License](#license)

## Introduction

Model Confidence and pyDock Energy are important metrics in the field of molecular docking and structural biology. Model Confidence helps in assessing the reliability of a predicted structure, while pyDock Energy evaluates the docking energy, providing insights into the stability and feasibility of the docking complex.

This notebook provides a step-by-step protocol for calculating these metrics from the outputs of pyDock, AlphaFold 2, and AlphaFold 3 using Python, with a focus on reproducibility and ease of use.

## Requirements

Before running the notebook, ensure you have the following dependencies installed:

- Python 3.x
- Jupyter Notebook or JupyterLab
- NumPy
- SciPy
- pyDock
- Matplotlib

You can install the necessary Python packages using pip:

```bash
pip install numpy scipy matplotlib jupyterlab
