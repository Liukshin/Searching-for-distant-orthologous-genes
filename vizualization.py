import matplotlib.pyplot as plt
import os
from DatabaseTool import *
from Alignment import *
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.path import Path

def create_table(file_result, db_file="orthodb_phac.txt",table_name="table_results_phac.csv",output_path="./"):

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



def update_fasta_from_df(fasta_file, df,output_file):

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


    # fig = plt.figure(figsize=(10, 10), dpi=300)
    # Phylo.draw(tree, do_show=False)
    # plt.savefig(output_image, bbox_inches='tight')
    # plt.close()


    plt.figure(figsize=(18, 14), dpi=400)
    ax = plt.gca()
    Phylo.draw(tree, axes=ax, do_show=False, branch_labels=None)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.savefig(output_image, bbox_inches='tight', pad_inches=0)
    plt.show()
    plt.close()

    print(f"Tree saved in: {output_newick}")
    print(f"Image saved in: {output_image}")


# def create_tree(
#         fasta_file: str,
#         output_newick="tree.nwk",
#         output_image="tree.png",
#         output_dir='data',
#         groups=None,
#         colors=None,
#         link_groups=None
#     ):
#
#     alignment = ClustalWAlignment(fasta_file)
#     newick_str = alignment.neighbor_joining()
#
#     output_newick = os.path.join(output_dir, output_newick)
#     output_image = os.path.join(output_dir, output_image)
#
#     with open(output_newick, "w") as f:
#         f.write(newick_str)
#
#     tree = Phylo.read(StringIO(newick_str), "newick")
#
#
#     descriptions = {}
#     with open(fasta_file, "r") as f:
#         for record in SeqIO.parse(f, "fasta"):
#             record_id = record.id.split(":")[0]
#             descriptions[record_id] = record.description
#
#     for clade in tree.get_terminals():
#         seq_id = clade.name.split(":")[0]
#         if seq_id in descriptions:
#             clade.name = descriptions[seq_id]
#
#     if groups and colors:
#         for clade in tree.get_terminals():
#             full_name = clade.name
#
#
#             for group_name, members in groups.items():
#                 if any(member in full_name for member in members):
#                     color = colors.get(group_name)
#
#                     for node in tree.get_path(clade):
#                         if not hasattr(node, 'color'):
#                             node.color = color
#                     break
#
#
#     fig = plt.figure(figsize=(12, 10), dpi=300)
#     ax = fig.add_subplot(111)
#
#
#     Phylo.draw(tree, axes=ax, do_show=False, branch_labels=None)
#
#     if groups and colors:
#         for group_name, members in groups.items():
#             color = colors.get(group_name, 'gray')
#             for clade in tree.get_terminals():
#                 if clade.name.split()[0] in members:
#
#                     for node in tree.get_path(clade):
#                         if hasattr(node, 'color'):
#                             node.color = color
#                         else:
#                             setattr(node, 'color', color)
#
#
#     Phylo.draw(tree, axes=ax, do_show=False, branch_labels=None)
#
#     if link_groups and groups:
#         for group1, group2 in link_groups:
#
#             pass
#
#
#     if groups and colors:
#         legend_patches = [
#             Patch(label=name, color=color)
#             for name, color in colors.items()
#         ]
#         ax.legend(
#             handles=legend_patches,
#             frameon=False,
#             bbox_to_anchor=(0.3, 0.3),
#             loc="center",
#             ncols=2
#         )
#
#     plt.savefig(output_image, bbox_inches='tight')
#     plt.show()
#     plt.close()
#
#     print(f"Tree saved in: {output_newick}")
#     print(f"Image saved in: {output_image}")




def plot_hits(all_hits, output_dir='data', name_graph = "all_iterations.png"):
    """
    Plots the hits from all iterations with different colors.
    :param all_hits: Dictionary where keys are iteration numbers and values are lists of Hit objects.
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
            label=f'Iteration {iteration}',
            color=colors[iteration]
        )

    plt.legend(frameon=True, framealpha=0.8, edgecolor='gray')
    plt.xlabel("Hit Position", fontsize=12)
    plt.ylabel("Bit Score", fontsize=12)
    plt.title("HMMER Hit Scores Across Iterations", fontsize=14)

    for spine in plt.gca().spines.values():
        spine.set_linewidth(0.5)
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)


    plt.savefig(os.path.join(output_dir, name_graph), bbox_inches='tight', dpi=300)
    plt.show()


