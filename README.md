# DockInspect

**DockInspect** is a Python tool for analyzing and visualizing protein–ligand docking results.  
It provides command-line utilities and visualization features built on top of PyMOL, RDKit and FreeSASA.

---

## Features

- Extract and summarize docking poses, compute hydrogen bonds
- Compute pocket properties like SASA, GRAVY, and net charge
- Compute ligand properties like SASA, logP, and net charge
- Visualize docking results with PyMOL
- Command-line interface for quick access

---

## Installation

Requires Python 3.8+ and PyMOL

```bash
git clone https://github.com/yourusername/dockinspect.git
cd dockinspect
pip install . 
``` 

## Tests

Tests can be run using
```bash
pytest tests
``` 