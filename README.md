# Validating Genuine Changes in Heartbeat Evoked Potentials using Pseudotrials and Surrogate Procedures

[![DOI](https://img.shields.io/badge/DOI-10.1162%2FIMAG.a.30-blue)](https://direct.mit.edu/imag/article/doi/10.1162/IMAG.a.30/130942/Validating-genuine-changes-in-Heartbeat-Evoked)
[![License](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a+-orange.svg)](https://www.mathworks.com/products/matlab.html)

This repository contains he code that accompanies the publication by Steinfath et al. (2025) in *Imaging Neuroscience*.

## Overview

Each heartbeat generates time-locked brain activity known as heartbeat evoked potentials (HEP). This project addresses a critical methodological challenge: distinguishing genuine heartbeat-related neural activity from heartbeat-independent confounding signals.

## Key Methods

### 1. Pseudotrial Correction

The pseudotrial correction method removes heartbeat-independent confounding activity by:

- Creating control trials with simulated heartbeat timings
- Subtracting pseudotrial activity from genuine HEP responses
- Preserving genuine heartbeat-related neural responses

### 2. Surrogate Heartbeat Analysis

Validates HEP time-locking by:

- Shuffling inter-beat intervals while preserving temporal structure
- Testing whether observed effects persist with disrupted heart-brain timing
- Confirming genuine cardiac-neural coupling

### 3. Simulation Framework


- Phase-randomized EEG data with controlled HEP-P300 relationships
- Multiple relationship types (direct, inverse, none)
- Systematic testing of analysis robustness

## Requirements

### Software Dependencies

- **MATLAB** (R2020a or later)
- **EEGLAB** (v2024.1 or later) with plugins:
  - `bva_io` 
  - `erplab` 
- **HEPLAB** 
- **FieldTrip** 

## Citation

If you use this code in your research, please cite:

```bibtex
@article{steinfath2025validating,
  title={Validating genuine changes in Heartbeat Evoked Potentials using Pseudotrials and Surrogate Procedures},
  author={Steinfath, Paul and Herzog, Nadine and Fourcade, Antonin and Sander, Christian and Nikulin, Vadim and Villringer, Arno},
  journal={Imaging Neuroscience},
  year={2025},
  doi={10.1162/IMAG.a.30},
  url={https://direct.mit.edu/imag/article/doi/10.1162/IMAG.a.30/130942/Validating-genuine-changes-in-Heartbeat-Evoked}
}
```

---

**Keywords:** Heartbeat Evoked Potential, Oddball Task, Pseudotrial Correction, Surrogate Heartbeat Analysis, Interoception, EEG, Heart-Brain Coupling
