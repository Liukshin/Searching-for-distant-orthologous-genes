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
```
# Special Cases
## OrthoDB_py
```bash
git clone https://gitlab.com/ezlab/orthodb_py.git
cd orthodb_py && pip install .
```
# pyhmmer (Linux-only)

```bash
pip install pyhmmer
```
# Usage

## Step 1: Fetch Protein Sequences

```python
from DatabaseTool import ProteinDatabaseHandlerNCBI

handler = ProteinDatabaseHandlerNCBI(df)  # df = DataFrame with organism names
handler.protein_search(email="your_email@example.com")
handler.download_protein(output_dir="data/", file_name="file.fasta")
```
## Step 2: Alignment

```python
alignment = ClustalWAlignment(file_name="file.fasta")
aligned_result = alignment.align()
alignment.save_alignment_to_fasta(aligned_result, output_file=os.path.join(prot_dir, "alignedfile.fasta"))
```

## Step 3: Run Iterative HMM Search

```python
from HMM import Model

model = Model(filename="file.fasta", dataset="your_dataset.fasta", output_folder="data/")
found_seqs = model.sequential_search(
    iterations=3,
    max_sequences=10,
    initial_threshold=700,
    combine_output=True
)
```
##Step 3: Visualize Results
```python
create_tree(output_file=os.path.join(prot_dir, "file.fasta"))
```
## License

MIT License. See LICENSE for details.


