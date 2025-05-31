from pandas import DataFrame
from DatabaseTool import *
from Alignment import *
import pandas as pd
import re
from collections import defaultdict
from Bio import SeqIO
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from io import StringIO
from phytreeviz import TreeViz
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from tempfile import NamedTemporaryFile
from Bio import Phylo
import numpy as np




def create_table(file_result, db_file:str,table_name: str,output_path: str)->DataFrame:
    """
    Save results to csv
    return DataFrame

    Parameters
    ----------
    file_result : str
        Name of fasta file.
    db_file : str
        Name of fasta file.
    table_name : str
        Name of csv file.
    output_path : str
    """

    target_ids = set()
    with open(file_result) as f1:
        for record in SeqIO.parse(f1, "fasta"):
            target_ids.add(record.id)

    records = []
    with open(db_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            record_id = record.id
            if record_id in target_ids:
                try:
                    json_str = record.description.split(" ", 1)[1]
                    data = json.loads(json_str)
                    data["sequence"] = str(record.seq)
                    data["id"] = record_id
                    records.append(data)
                except Exception as e:
                    print(f"Error {record_id}: {e}")


    df = pd.DataFrame(records)
    df.to_csv(os.path.join(output_path,table_name), index=False)
    return df



def update_fasta_from_df(fasta_file:str, df:DataFrame,output_file:str):
    """
    Overwrites fasta file with organism names

    Parameters
    ----------
    fasta_file : str
        Name of fasta file.
    df : DataFrame
        DataFrame containing organism names.
    output_file : str
        Name of fasta file.

    """

    updated_records = []
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            record_id = record.id.split(":")[0]
            #record_id = record.description.split()[0].split(":")[0]
            organism_name = df.loc[df['organism_taxid'] == record_id, 'organism_name'].values

            if len(organism_name) > 0:

                record.description = f"{organism_name[0]}"
                record.id = record_id
            else:
                print(f"Warning: ID {record_id} not found")

            updated_records.append(record)


    with open(output_file, "w") as f:
        SeqIO.write(updated_records, f, "fasta")

    print(f"'{fasta_file}' upload")





def create_tree(fasta_file: str, output_newick="tree.nwk", output_image="tree.png", output_dir='data'):
    """
    Simple visualization of a phylogenetic tree

    Parameters
    ----------
    fasta_file : str
        Name of fasta file.
    output_newick : str
        Newick file
    output_image : str
        Name of .png file.
    output_dir : str
        Name of directory file
    """

    alignment = ClustalWAlignment(fasta_file)
    newick_str = alignment.neighbor_joining(fyl_tree_viz=True)

    output_newick = os.path.join(output_dir, output_newick)
    output_image = os.path.join(output_dir, output_image)


    with open(output_newick, "w") as f:
        f.write(newick_str)


    tree = Phylo.read(StringIO(newick_str), "newick")


    descriptions = {}
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            record_id = record.id.split(":")[0]
            descriptions[record_id] = record.description


    for clade in tree.get_terminals():
        seq_id = clade.name.split(":")[0]
        if seq_id in descriptions:
            clade.name = descriptions[seq_id]


    plt.figure(figsize=(14, 12), dpi=300)
    ax = plt.gca()
    Phylo.draw(tree, axes=ax, do_show=False,
               branch_labels=lambda c: f"{c.branch_length:.3f}" if c.branch_length else "")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.savefig(output_image, bbox_inches='tight', pad_inches=0)
    plt.show()
    plt.close()

    print(f"Tree saved in: {output_newick}")
    print(f"Image saved in: {output_image}")


def clusters_tree(newick_str, fasta_file, output_image="colored_clusters.png", path='data',
                  cluster_threshold=0.3):
    """
    Colored clusters visualization of a phylogenetic tree

    Parameters
    ----------
    newick_str : str
        Newick file
    fasta_file : str
        Name of fasta file
    output_image : str
        Name of .png file.
    path : str
        Name of directory file
    cluster_threshold: float
        Affects the number of clusters
    Returns
        -------
        clusters:dict
            Dictionary of clusters where:
            - Keys are cluster identifiers (integers),
            - Values are lists of sequence labels grouped in that cluster.
    """

    tree = Phylo.read(StringIO(newick_str), "newick")

    id_to_description = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split(":")[0]
        id_to_description[seq_id] = record.description

    for clade in tree.get_terminals():
        seq_id = clade.name.split(":")[0]
        if seq_id in id_to_description:
            clade.name = id_to_description[seq_id]

    with NamedTemporaryFile("w+", delete=False, suffix=".nwk", dir=path) as tmp_file:
        Phylo.write(tree, tmp_file.name, "newick")
        tree_path = tmp_file.name

    terminals = tree.get_terminals()
    labels = [clade.name for clade in terminals]
    distances = np.zeros((len(terminals), len(terminals)))
    for i, clade1 in enumerate(terminals):
        for j, clade2 in enumerate(terminals):
            if i < j:
                d = tree.distance(clade1, clade2)
                distances[i, j] = d
                distances[j, i] = d

    condensed_dist = pdist(distances)
    linkage_matrix = linkage(condensed_dist, method="average")
    cluster_ids = fcluster(linkage_matrix, cluster_threshold, criterion='distance')

    clusters = {}
    for label, cluster_id in zip(labels, cluster_ids):
        clusters.setdefault(cluster_id, []).append(label)

    cmap = plt.colormaps["tab20"]
    color_list = [cmap(i) for i in range(len(clusters))]

    tv = TreeViz(tree_path, leaf_label_size=10, height=0.3)
    tv.show_scale_bar()
    tv.show_branch_length(color="black")

    legend_handles = []
    for i, (cluster_id, names) in enumerate(clusters.items()):
        color = color_list[i]
        tv.highlight(names, color=color)
        legend_handles.append(Patch(color=color, label=f"Shluk {cluster_id}"))

    fig = tv.plotfig()

    fig.legend(handles=legend_handles, frameon=False, bbox_to_anchor=(0.4, 0.05), loc="lower center", ncols=3)
    output_name = os.path.join(path, output_image)
    fig.savefig(output_name, dpi=600)
    plt.show()
    print(f"Saved in : {output_name}")
    return clusters


def plot_hits(all_hits, output_dir='data', name_graph = "all_iterations.png"):
    """
    Plot bit scores of HMMER hits across multiple iterations.

    Parameters
    ----------
    all_hits : dict
        Dictionary with iteration numbers as keys and lists of Hit objects as values.
    output_dir : str, optional
        Directory to save the output image (default is 'data').
    name_graph : str, optional
        Name of the output .png file (default is 'all_iterations.png').

    """
    plt.style.use('seaborn-v0_8')
    plt.figure(figsize=(10, 6), dpi=300)

    colors = plt.cm.viridis(np.linspace(0, 1, len(all_hits)))

    for iteration, hits in all_hits.items():
        scores = [hit.score for hit in hits]
        x_values = range(1, len(scores) + 1)

        plt.plot(
            x_values,
            scores,
            'o-',  
            linewidth=1.5,
            markersize=6,
            markeredgewidth=0.5,
            label=f'Iterace {iteration+1}',
            color=colors[iteration]
        )

    plt.legend(frameon=True, framealpha=0.8, edgecolor='gray')
    plt.xlabel("Pozice zásahu", fontsize=12)
    plt.ylabel("Bitové skóre", fontsize=12)
    plt.title("Skóre zásahů HMMER při různých iteracích", fontsize=14)

    for spine in plt.gca().spines.values():
        spine.set_linewidth(0.5)
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)


    plt.savefig(os.path.join(output_dir, name_graph), bbox_inches='tight', dpi=300)
    plt.show()



def find_custom_motif(fasta_file, pattern):
    """
    Search for a custom motif in a FASTA file.

    Parameters
    ----------
    fasta_file : str
        Name of fasta file.
    pattern : str
        Motif pattern (e.g., "[GS]-X-C-X-[GA]-G").

    Returns
    -------
    stats : dict
        Dictionary with the following keys:
            - total_sequences : int
                Number of sequences in the input file.
            - sequences_with_motif : int
                Number of sequences containing the motif.
            - total_motifs : int
                Total number of motif occurrences found.
            - motifs_per_seq : dict
                Number of motifs found per sequence (seq_id → count).
            - positions : dict
                Motif positions per sequence (seq_id → list of positions).
            - organisms : dict
                Organism name per sequence (seq_id → name).
            - pattern_used : str
                Motif pattern used for the search.
    """
    regex_pattern = pattern.replace("X", ".").replace("-", "")
    regex_pattern = f"(?=({regex_pattern}))"

    records = list(SeqIO.parse(fasta_file, "fasta"))

    stats = {
        'total_sequences': len(records),
        'sequences_with_motif': 0,
        'total_motifs': 0,
        'motifs_per_seq': defaultdict(int),
        'positions': defaultdict(list),
        'organisms': {},
        'pattern_used': regex_pattern
    }

    for record in records:
        seq = str(record.seq)
        seq_id = record.id
        # Try to extract organism name from description
        description_parts = record.description.split(maxsplit=1)
        organism = description_parts[1] if len(description_parts) > 1 else "Unknown"
        stats['organisms'][seq_id] = organism

        matches = re.finditer(regex_pattern, seq)
        count = 0
        for match in matches:
            start = match.start() + 1
            motif = match.group(1)
            stats['positions'][seq_id].append((start, motif))
            count += 1
        if count > 0:
            stats['sequences_with_motif'] += 1
            stats['motifs_per_seq'][seq_id] = count
            stats['total_motifs'] += count

    stats['avg_motifs_per_seq'] = (
        stats['total_motifs'] / stats['total_sequences']
        if stats['total_sequences'] > 0 else 0
    )

    return stats