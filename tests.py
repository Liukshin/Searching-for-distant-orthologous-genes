import os.path

from DatabaseTool import *
from Alignment import *
from HMM import *
import pandas as pd
from vizualization import *


prot_dir = r'data'
dataset = os.path.join(prot_dir, 'orthodb_phac.fasta')

#handler = ProteinDatabaseHandlerNCBI(None)
#create_tree(os.path.join(prot_dir, "all_iters_combined.fasta"), db_handler=handler)

#create_tree(os.path.join(prot_dir, 'aligned_final.fasta'), handler)

# dataset = os.path.join(prot_dir, 'orthodb_phac.fasta')
# model = Model(filename='phac_cupr', dataset= dataset, output_folder=prot_dir)
#
# hits_v = model.iterative_search(iterations=2, max_sequences=5)
#
#
#
# tree = create_colored_tree(hits_v)
#
#
# Phylo.draw(tree)
# plt.show()
#
#
#
# print(hits_v)
#
# for iteration, hits in hits_v.items():
#     # Extract scores from Hit objects
#     scores = [hit.score for hit in hits]
#     print([hit for hit in hits])
#     print(scores)




# test_fasta = os.path.join(prot_dir, 'all_iters_combined.fasta')
# df = create_table(test_fasta,dataset)
#
# #update_fasta_from_df(test_fasta,df)
# handler = ProteinDatabaseHandlerNCBI(None)
# create_tree(test_fasta)

import re
from collections import defaultdict

from Bio import SeqIO
import re
from collections import defaultdict

# dataset123 = os.path.join(prot_dir, 'orthodb123.fasta')
# file_name_list = ["phac1.fasta", "phac2.fasta", "phac3.fasta","phac4.fasta"]
# merge_unique_fasta(file_name_list,prot_dir, dataset123)

# orthodb_url = "https://data.orthodb.org/current/fasta?id=122391at1224&seqtype=protein&species="
# download_dataset_url(orthodb_url, os.path.join(prot_dir, "phac1234.fasta"))
#
from Bio import SeqIO
import os

input_fasta = os.path.join(prot_dir, 'orthodb123.fasta')
output_fasta = os.path.join(prot_dir, 'filtered_orthodb123.fasta')

unique_organisms = {}
input_count = 0

with open(input_fasta, "r") as in_handle:
    for record in SeqIO.parse(in_handle, "fasta"):
        input_count += 1
        header = record.description
        organism = None
        if '"organism_name":"' in header:
            try:
                organism = header.split('"organism_name":"')[1].split('"')[0]
            except IndexError:
                continue

        if organism and organism not in unique_organisms:
            unique_organisms[organism] = record


with open(output_fasta, "w") as out_handle:
    SeqIO.write(unique_organisms.values(), out_handle, "fasta")


original_size = os.path.getsize(input_fasta) / (1024 * 1024)
filtered_size = os.path.getsize(output_fasta) / (1024 * 1024)
filtered_count = len(unique_organisms)

print(f"seq: {input_count}")
print(f"after filt: {filtered_count}")
print(f"size: {original_size:.2f} MB")
print(f"after filt: {filtered_size:.2f} MB")


# def motiv_search(file='update_phac.fasta',motiv = "[GS]-X-C-X-[GA]-G"):
#
#     stats = find_custom_motif(
#         fasta_file=file,
#         pattern=motiv
#     )
#
#     print(f"Pattern used: {stats['pattern_used']}")
#     print(f"Total sequences: {stats['total_sequences']}")
#     print(f"Sequences with motif: {stats['sequences_with_motif']}")
#     print(f"Total motif occurrences: {stats['total_motifs']}")
#     print(f"Average motifs per sequence: {stats['avg_motifs_per_seq']:.2f}")
#
#     print("\nMotif positions by sequence:")
#     for seq_id, motifs in stats['positions'].items():
#         organism = stats['organisms'].get(seq_id, "Unknown")
#         print(f"{seq_id} ({organism}): {motifs}")
#
#
#
#
# dataset4 = os.path.join(prot_dir, 'filtered_orthodb123.fasta')
# # result processing
# test_fasta = os.path.join(prot_dir, f'aligned_final1.fasta')
# df_t1 = create_table(test_fasta, dataset4, output_path=prot_dir)
# updated_fasta = os.path.join(prot_dir, f'update_phac1.fasta')
# update_fasta_from_df(test_fasta, df_t1, updated_fasta)
# create_tree(os.path.join(prot_dir, f'update_phac1.fasta'))
#
# motiv_search(file=os.path.join(prot_dir, f'update_phac1.fasta'))
#
# update_fasta_files = [os.path.join(prot_dir, f'update_phac1.fasta'),os.path.join(prot_dir, f'update_phac2.fasta'),os.path.join(prot_dir, f'update_phac3.fasta'),os.path.join(prot_dir, f'update_phac4.fasta')]
#
# combined_output = os.path.join(prot_dir, "final_combined_update_phac.fasta")
#
# with open(combined_output, "w") as outfile:
#     for fasta_file in update_fasta_files:
#         for record in SeqIO.parse(fasta_file, "fasta"):
#             SeqIO.write(record, outfile, "fasta")


from Bio import SeqIO
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from io import StringIO
from phytreeviz import TreeViz
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from Bio import Phylo
import numpy as np

# rabotajet
def visualize_clusters_with_phytreeviz(newick_str, fasta_file, output_image="colored_clusters.png",
                                       cluster_threshold=0.3):

    tree = Phylo.read(StringIO(newick_str), "newick")


    id_to_description = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split(":")[0]
        id_to_description[seq_id] = record.description


    for clade in tree.get_terminals():
        seq_id = clade.name.split(":")[0]
        if seq_id in id_to_description:
            clade.name = id_to_description[seq_id]


    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile("w+", delete=False, suffix=".nwk") as tmp_file:
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
    fig.savefig(output_image, dpi=600)
    plt.show()
    print(f"Saved in : {output_image}")
    return clusters


#u vsech bylo 0.35 u 3 0.6
fasta_file = f"data/update_phac4.fasta"
alignment = ClustalWAlignment(fasta_file)
newick_str = alignment.neighbor_joining(fyl_tree_viz=True)

clust_dict = visualize_clusters_with_phytreeviz(newick_str, fasta_file, output_image=f"data/colored_clusters4.png", cluster_threshold=0.35)
print(clust_dict)


