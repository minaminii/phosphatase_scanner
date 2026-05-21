# Phosphatase Scanner

**HMM-based phosphatase detection with motif-aware report generation**

## Overview

**Phosphatase Scanner** is a HMM analysis pipeline designed to improve functional interpretation of phosphatase-like proteins. While HMM-based tools such as hmmscan are powerful for detecting phosphatase-related folds, they cannot determine whether key catalytic motifs are intact, relaxed, or degenerated.

This project combines:

1. HMMER hmmscan–based fold detection
2. Automated hit summarization
3. Motif-level pattern matching for functional validation
4. User-friendly reporting and exploration via Jupyter Notebook


#### Requirements

- HMMER3 (hmmscan) with attached hmm profile (in `data` directory)
- Python 3
    - tqdm
    - pandas
    - jinja2
    - biopython
- hmmscan_summary.py

## Installation

```bash
# Download the script
git clone https://github.com/minaminii/phosphatase_scanner.git

# Install Python 3 and Python package manager (pip)
sudo apt update
sudo apt install python3 python3-pip

# Install hmmer for hmmscan (from http://hmmer.org/documentation.html)
apt install hmmer

# Install Python dependencies
pip3 install tqdm pandas jinja2 biopython
```

## Usage

```bash
hmmscan_pipeline.sh -i proteins.faa -p sample1 -o results/
```

## Citation & Contact

**Citation:**

Pongpom, M., Wattanasombat, S., Aphiwongcharoen, C., Limsamutpet, K., & Wangsanut, T. (2026). Bridging manual and computational approaches: The Pseudophosphatase Scanner for Genome-Wide fungal pseudophosphatome analysis. ACS Omega. https://doi.org/10.1021/acsomega.6c00759

**Author:**

Sara Wattabasombat
Department of Microbiology, Faculty of Medicine, Chiang Mai University
