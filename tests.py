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




test_fasta = os.path.join(prot_dir, 'all_iters_combined.fasta')
df = create_table(test_fasta,dataset)

#update_fasta_from_df(test_fasta,df)
handler = ProteinDatabaseHandlerNCBI(None)
create_tree(test_fasta)
