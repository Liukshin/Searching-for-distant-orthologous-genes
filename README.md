# Searching for Distant Orthologous Genes

## Description
This project performs iterative search for distant orthologs using profile Hidden Markov Models (HMMs). After each iteration, the algorithm constructs an improved HMM from the detected sequences. The final results are visualized as a phylogenetic tree, with optional clustering of sequences into orthologous groups.

Key features:
- **Iterative refinement**: Updates HMMs dynamically to enhance sensitivity.
- **Phylogenetic visualization**: Uses `phytreeviz` to display evolutionary relationships.
- **Integration with NCBI**: Fetches protein sequences via API requests.
- **Modular design**: Separates alignment, HMM generation, and visualization logic.

## Installation

### Prerequisites
- **Linux** (required for `pyhmmer`)
- **Python 3.9+**

### Dependencies
Install core packages via pip:
```bash
pip install biopython numpy requests pandas scipy matplotlib
