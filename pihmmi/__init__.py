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


# Import main pipelines
from .pipelines import (
    create_dataset,
    search_protein,
    global_alignment,
    multiple_alignment,
    motif_search,
    run_pihmmi_pipeline
)

# Conditional HMM import (requires pyhmmer)
try:
    from .hmm import Model

    HMM_AVAILABLE = True
except ImportError:
    Model = None
    HMM_AVAILABLE = False

# Check OrthoDB availability
try:
    import orthodb

    ORTHODB_AVAILABLE = True
except ImportError:
    ORTHODB_AVAILABLE = False


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
    "__version__",
    "HMM_AVAILABLE",
    "ORTHODB_AVAILABLE"
]


def check_dependencies():
    """Check availability of optional dependencies"""
    deps = {
        "pyhmmer (for HMM)": HMM_AVAILABLE,
        "orthodb-py (for OrthoDB)": ORTHODB_AVAILABLE,
        "core dependencies": True
    }

    print("Dependency Status:")
    for dep, available in deps.items():
        status = "✓" if available else "✗"
        print(f"  {status} {dep}")

    if not HMM_AVAILABLE:
        print("\nFor full HMM functionality install: conda install -c bioconda pyhmmer")
    if not ORTHODB_AVAILABLE:
        print("For OrthoDB support install: pip install git+https://gitlab.com/ezlab/orthodb_py.git")

    return deps