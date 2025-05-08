import matplotlib.pyplot as plt
import os
from DatabaseTool import *
from Alignment import *


def create_tree(fasta_file: str, db_handler: ProteinDatabaseHandlerNCBI,
                           output_newick="tree.nwk", output_image="tree.png", output_dir ='data'):

    alignment = ClustalWAlignment(fasta_file)
    newick_str = alignment.neighbor_joining()


    with open(output_newick, "w") as f:
        f.write(newick_str)


    tree = Phylo.read(StringIO(newick_str), "newick")
    fig = plt.figure(figsize=(10, 10), dpi=300)
    Phylo.draw(tree, do_show=False)
    plt.savefig(output_image, bbox_inches='tight')
    plt.close()

    output_newick_path = os.path.join(output_dir, output_newick)
    output_image_path = os.path.join(output_dir, output_image)

    print(f"Tree save in: {output_newick_path}")
    print(f"Image save in: {output_image_path}")