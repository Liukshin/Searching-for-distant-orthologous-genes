# Searching for Distant Orthologous Genes

## Description
This project performs iterative search for distant orthologs using profile Hidden Markov Models (HMMs). After each iteration, the algorithm constructs an improved HMM from the detected sequences. The final results are visualized as a phylogenetic tree, with optional clustering of sequences into orthologous groups.

Key features:
- **Iterative refinement**: Updates HMMs dynamically to enhance sensitivity.
- **Phylogenetic visualization**: Uses `phytreeviz` to display evolutionary relationships.
- **Integration with NCBI**: Fetches protein sequences via API requests.
- **Modular design**: Separates alignment, HMM generation, and visualization logic.

---

##  pihmmi Library

`pihmmi` is a standalone Python library used in this project. It provides programmatic access to protein sequence analysis, multiple alignment using HMMs, and distant ortholog search.  
It is built with modern scientific packages and is installable via [Poetry](https://python-poetry.org/).

---

## Installation

### Prerequisites
- **Linux** (required for `pyhmmer`)
- **Python 3.11+**
###  Recommended (using Poetry)

Make sure you have [Poetry](https://python-poetry.org/docs/#installation) installed.

Then, in the project root (where `pyproject.toml` and `poetry.lock` are located):

```bash
poetry install
```

### Dependencies
You can also install all dependencies manually using pip:
```bash
pip install pyhmmer biopython numpy requests pandas scipy matplotlib phytreeviz scikit-learn
```
### Special Cases
#### OrthoDB_py
```bash
git clone https://gitlab.com/ezlab/orthodb_py.git
cd orthodb_py && pip install .
```
---
# Usage

## Step 1: Fetch Protein Sequences

```python
from pihmmi import ProteinDatabaseHandlerNCBI
handler = ProteinDatabaseHandlerNCBI(df)  # df = DataFrame with organism names
handler.protein_search(email="your_email@example.com")
handler.download_protein(output_dir="data/", file_name="file.fasta")
```
## Step 2: Alignment

```python
from pihmmi import ClustalWAlignment
alignment = ClustalWAlignment(file_name="file.fasta")
aligned_result = alignment.align()
alignment.save_alignment_to_fasta(aligned_result, output_file=os.path.join(prot_dir, "alignedfile.fasta"))
```

## Step 3: Run Iterative HMM Search

```python
from pihmmi import Model

model = Model(filename="file.fasta", dataset="your_dataset.fasta", output_folder="data/")
found_seqs = model.sequential_search(
    iterations=3,
    max_sequences=10,
    initial_threshold=700,
    combine_output=True
)
```
## Step 4: Visualize Results

```python
from pihmmi import create_tree
create_tree(output_file=os.path.join(prot_dir, "file.fasta"))
```
## License

MIT License. See LICENSE for details.


