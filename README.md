# Validating Genuine Changes in Heartbeat Evoked Potentials using Pseudotrials and Surrogate Procedures


[![License](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a+-orange.svg)](https://www.mathworks.com/products/matlab.html)

This repository contains he code that accompanies this publication: [![DOI](https://img.shields.io/badge/DOI-10.1162%2FIMAG.a.30-blue)](https://direct.mit.edu/imag/article/doi/10.1162/IMAG.a.30/130942/Validating-genuine-changes-in-Heartbeat-Evoked)

## Overview

Each heartbeat generates time-locked brain activity known as heartbeat evoked potentials (HEP). This project addresses the methodological challenge of distinguishing genuine heartbeat-related neural activity from heartbeat-independent confounding signals.

## Key Methods

### 1. Pseudotrial Correction

The pseudotrial correction method removes heartbeat-independent confounding activity by:

- Creating pseudotrials with random heartbeat timings
- Subtracting pseudotrial activity from HEP responses
- Ideally preserving/uncovering genuine heartbeat-related neural responses

### 2. Surrogate Heartbeat Analysis

Validates HEP time-locking by:

- Shuffling R-peak triggers while preserving temporal structure
- Repeating statistics to test whether effects of similar strength persist with disrupted heart-brain timing
- Ideally confirming genuine cardiac-neural associations

## Requirements

### Software Dependencies

- **MATLAB** 
- **EEGLAB** with plugins:
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
