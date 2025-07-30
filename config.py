import os
import pandas as pd

prot_dir = r'data'
os.makedirs(prot_dir, exist_ok=True)

data_0iter = {
    'Source': [
        'Caulobacter crescentus',
        'Cupriavidus necator',
        'Paracoccus denitrificans',
        'Vibrio cholerae',
        'Chromobacterium violaceum',
        'Hyphomonas sp. CACIAM 19H1',
        'Pseudomonas aeruginosa',
        'Pseudomonas fluorescens',
        'Pseudomonas putida',
        'Haloferax mediterranei',
        'Haloquadratum walsbyi',
        'Synechocystis sp. PCC 6803',
        'Bacillus cereus',
        'Bacillus megaterium'
    ],
    'Gene': [
        'PhaC',
        'PhaC',
        'PhaC',
        'PhaC',
        'PhaC',
        'PhaC',
        'PhaC',
        'PhaC',
        'PhaC',
        'PhaC',
        'PhaC',
        'PhaC',
        'PhaC',
        'PhaC'
    ]
}

df = pd.DataFrame(data_0iter)

df_class_1 = df.iloc[0:4].reset_index(drop=True)
df_class_2 = df.iloc[4:9].reset_index(drop=True)
df_class_2v2 = df.iloc[4:8].reset_index(drop=True)
df_class_3 = df.iloc[10:12].reset_index(drop=True)
df_class_4 = df.iloc[12:14].reset_index(drop=True)


dataset123 = 'orthodb123.fasta'
# dataset123 = os.path.join(prot_dir, 'phac1234.fasta')  #
# dataset123 = os.path.join(prot_dir, 'orthodb_phac.fasta')  #

dataset1 = os.path.join(prot_dir, 'filtered_phac1.fasta')
dataset2 = os.path.join(prot_dir, 'filtered_phac2.fasta')
dataset3 = os.path.join(prot_dir, 'filtered_phac3.fasta')
dataset4 = os.path.join(prot_dir, 'filtered_orthodb123.fasta')


gen_name_phac = ["PhaC class I", "PhaC class II", "PhaC class III"]
file_name_phac = ["phac1.fasta", "phac2.fasta", "phac3.fasta", "phac4.fasta"]

uniport_url_phac4 = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28PhaC+class+IV%29"
orthodb_url_phac4 ="https://data.orthodb.org/current/fasta?id=122391at1224&seqtype=protein&species="

uniport_url_file_name = "phac4.fasta"
orthodb_url_file_name = "phac1234.fasta"

search_file = 'phac_cupr.fasta'

