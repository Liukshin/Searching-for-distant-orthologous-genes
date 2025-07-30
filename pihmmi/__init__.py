"""
PIHMMI - Bioinformatics library for protein sequence analysis
Based on distant orthologous gene search project
"""

__version__ = "1.6.0"
__author__ = "Maksim Liukshin"
__email__ = "247034@vut.cz"


from .alignment import GlobalAlignment, ClustalWAlignment, MultipleAlignment
from .database_tool import (
    ProteinDatabaseHandlerNCBI,
    download_orhodb_dataset,
    merge_unique_fasta,
    download_dataset_url,
    filter_dataset
)
from .visualization import (
    create_table,
    update_fasta_from_df,
    create_tree,
    clusters_tree,
    plot_hits,
    find_custom_motif,
    find_optimal_threshold_newick
)

from .pipelines import (
    create_dataset,
    search_protein,
    global_alignment,
    multiple_alignment,
    motif_search,
    run_pihmmi_pipeline
)

from .hmm import Model


__all__ = [
    # Main pipelines
    "create_dataset",
    "search_protein",
    "global_alignment",
    "multiple_alignment",
    "motif_search",
    "run_pihmmi_pipeline",

    # Classes for advanced usage
    "GlobalAlignment",
    "ClustalWAlignment",
    "ProteinDatabaseHandlerNCBI",
    "Model",

    # Utilities
    "create_tree",
    "clusters_tree",
    "find_custom_motif",
    "plot_hits",

    # Information
    "__version__"
]

