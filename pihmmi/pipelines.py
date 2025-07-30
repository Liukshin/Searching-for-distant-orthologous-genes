"""
Main pipelines for protein sequence analysis
"""
import os
import pandas as pd
from typing import Optional, Dict, Any
from pandas import DataFrame
from .hmm import Model

from .database_tool import (
    ProteinDatabaseHandlerNCBI,
    download_orhodb_dataset,
    merge_unique_fasta,
    download_dataset_url,
    filter_dataset
)
from .alignment import GlobalAlignment, ClustalWAlignment
from .visualization import (
    create_table,
    create_tree,
    update_fasta_from_df,
    find_custom_motif,
    clusters_tree,
    find_optimal_threshold_newick
)

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
from Bio import Phylo
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import pdist


def create_dataset(output_dir: str = None,
                   dataset_name: str = None,
                   gen_name_list: list = None,
                   file_name_list: list = None,
                   uniport_url: str = None,
                   orthodb_url: str = None,
                   additional_uniport_file: str = None,
                   additional_orthodb_file: str = None):
    """
    Create dataset from various sources (OrthoDB, UniProt)

    Parameters
    ----------
    output_dir : str, optional
        Directory to save files (default uses prot_dir)
    """


    os.makedirs(output_dir, exist_ok=True)


    print("Downloading data from OrthoDB...")
    for i in range(len(gen_name_list)):
        download_orhodb_dataset(gen_name_list[i], os.path.join(output_dir, file_name_list[i]))

    print("Downloading data from UniProt...")
    download_dataset_url(uniport_url, os.path.join(output_dir, additional_uniport_file))

    print("Downloading data from OrthoDB (combined set)...")
    download_dataset_url(orthodb_url, os.path.join(output_dir, additional_orthodb_file))

    print("Merging files...")
    input_fasta = os.path.join(output_dir, dataset_name)
    merge_unique_fasta(file_name_list, output_dir, input_fasta)

    print("Filtering dataset...")
    input_fasta = os.path.join(output_dir, dataset_name)
    output_fasta = os.path.join(output_dir, f'filtered_{dataset_name}')
    filter_dataset(input_fasta, output_fasta)


    print(f"Dataset created in directory: {output_dir}")


def search_protein(df: DataFrame, email: str = "247034@vut.cz", output_dir: str = None, file_name: str = 'phac_cupr'):
    """
    Search and download protein sequences from NCBI

    Parameters
    ----------
    df : DataFrame
        DataFrame with 'Source' and 'Gene' columns
    email : str
        Email for NCBI Entrez
    output_dir : str, optional
        Directory to save files (default prot_dir)
    file_name : str
        Output file name

    Returns
    -------
    str
        Path to saved FASTA file
    """


    handler = ProteinDatabaseHandlerNCBI(df)
    handler.protein_search(email=email)
    handler.download_protein(output_dir=output_dir, file_name=file_name)



def global_alignment(name_file:str,name_dir:str):
    fasta_file_phac = os.path.join(name_dir, name_file)
    seqs = list(SeqIO.parse(fasta_file_phac, 'fasta'))
    alignment = GlobalAlignment(seqs[0],seqs[1],gap= -2)
    sek1,sek2 = alignment.align()
    record1 = SeqRecord(Seq(sek1), id="seq1")
    record2 = SeqRecord(Seq(sek2), id="seq2")
    output_path = os.path.join(name_dir,f"aligned{name_file}.fasta")

    SeqIO.write([record1, record2], output_path, "fasta")

def multiple_alignment(name_file:str,name_dir:str):
    fasta_file_phac = os.path.join(name_dir, name_file)
    alignment = ClustalWAlignment(file_name=fasta_file_phac)
    aligned_result = alignment.align()
    alignment.save_alignment_to_fasta(aligned_result, output_file=os.path.join(name_dir,f"aligned{name_file}.fasta"))


def motif_search(input_file: str, motif: str = "[GS]-X-C-X-[GA]-G") -> Dict[str, Any]:
    """
    Motif search pipeline in protein sequences

    Parameters
    ----------
    input_file : str
        FASTA file to search
    motif : str
        Motif pattern

    Returns
    -------
    Dict[str, Any]
        Statistics of found motifs
    """
    stats = find_custom_motif(fasta_file=input_file, pattern=motif)

    print(f"Pattern used: {stats['pattern_used']}")
    print(f"Total sequences: {stats['total_sequences']}")
    print(f"Sequences with motif: {stats['sequences_with_motif']}")
    print(f"Total motifs found: {stats['total_motifs']}")
    print(f"Average motifs per sequence: {stats['avg_motifs_per_seq']:.2f}")

    print("\nMotif positions by sequence:")
    for seq_id, motifs in stats['positions'].items():
        organism = stats['organisms'].get(seq_id, "Unknown")
        print(f"{seq_id} ({organism}): {motifs}")

    return stats





def run_pihmmi_pipeline(
        df_class: DataFrame,
        number_pipeline: int,
        email: str = "247034@vut.cz",
        model_name: str = 'phac_cupr',
        iterations: int = 5,
        max_seq: int = 7,
        multi_alignment: bool = True,
        threshold_search: int = 500,
        cluster_tree: bool = True,
        threshold_viz: float = 0.35,
        viz_hit: bool = False,
        silhouette_analysis: bool = False,
        data_set: str = None,
        output_dir: str = None
) -> Dict[str, Any]:
    """
    Main PIHMMI analysis pipeline

    Parameters
    ----------
    df_class : DataFrame
        DataFrame with protein data
    number_pipeline : int
        Pipeline number
    email : str
        Email for NCBI
    iterations : int
        Number of HMM search iterations
    max_seq : int
        Maximum number of sequences
    multi_alignment : bool
        Use multiple alignment
    threshold_search : int
        Search threshold
    threshold_viz : float
        Visualization threshold
    viz_hit : bool
        Show hit visualization
    silhouette_analysis : bool
        silhouette analysis for find optimal clustering threshold
    data_set : str, optional
        Path to dataset
    output_dir : str, optional
        Output directory

    Returns
    -------
    Dict[str, Any]
        Analysis results
    """

    results = {}

    # Step 1: Protein search
    print(f"=== Pipeline {number_pipeline}: Protein Search ===")
    protein_file = search_protein(df_class, email, output_dir)
    results['protein_file'] = protein_file

    # Step 2: Alignment
    print("=== Sequence Alignment ===")

    if multi_alignment:
        #aligned_file = multiple_alignment_pipeline(model_name, output_dir)
        multiple_alignment(name_dir=output_dir,name_file=model_name)
    else:
        #aligned_file = global_alignment_pipeline(model_name, output_dir)
        global_alignment(name_dir=output_dir,name_file=model_name)


    # Step 3: HMM search
    print("=== HMM Search ===")
    model = Model(filename=model_name, dataset=data_set, output_folder=output_dir)

    found_seqs = model.sequential_search(
        iterations=iterations,
        max_sequences=max_seq,
        initial_threshold=threshold_search,
        combine_output=True,
        visualization=viz_hit,
        combine_output_file=f'aligned_final{number_pipeline}.fasta'
    )
    results['found_sequences'] = found_seqs

    # Step 4: Process results
    print("=== Processing Results ===")
    test_fasta = os.path.join(output_dir, f'aligned_final{number_pipeline}.fasta')
    df_t1 = create_table(
        file_result=test_fasta,
        db_file=data_set,
        table_name=f"table_results_phac{number_pipeline}.csv",
        output_path=output_dir
    )

    updated_fasta = os.path.join(output_dir, f'update_phac{number_pipeline}.fasta')
    update_fasta_from_df(test_fasta, df_t1, updated_fasta)
    results['updated_fasta'] = updated_fasta

    # Step 5: Motif search
    print("=== Motif Search ===")
    motif_stats = motif_search(updated_fasta)

    results['motif_stats'] = motif_stats

    # Step 6: Phylogenetic analysis
    print("=== Phylogenetic Analysis ===")
    if not cluster_tree:
        create_tree(os.path.join(output_dir, f'update_phac{number_pipeline}.fasta'))
    else:
        alignment = ClustalWAlignment(updated_fasta)
        newick_str = alignment.neighbor_joining(fyl_tree_viz=True)

        if silhouette_analysis:
            optimal_threshold = find_optimal_threshold_newick(newick_str)
            print(f"Optimal clustering threshold: {optimal_threshold:.3f}")
            threshold_viz = optimal_threshold

        # Create tree and clustering
        clust_dict = clusters_tree(
            newick_str,
            updated_fasta,
            output_image=f"colored_clusters{number_pipeline}.png",
            path=output_dir,
            cluster_threshold=threshold_viz
        )
        results['clusters'] = clust_dict

    print(f"=== Pipeline {number_pipeline} completed ===")
    return results
