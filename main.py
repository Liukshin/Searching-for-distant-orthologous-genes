import os.path
from DatabaseTool import *
from Alignment import *
from HMM import *
import pandas as pd
from vizualization import *


def create_dataset():
    dataset123 = os.path.join(prot_dir, 'orthodb123.fasta')
    gen_name_list = ["PhaC class I", "PhaC class II", "PhaC class III"]
    file_name_list = ["phac1.fasta", "phac2.fasta", "phac3.fasta","phac4.fasta"]
    for i in range(len(gen_name_list)):
        download_orhodb_dataset(gen_name_list[i],os.path.join(prot_dir,file_name_list[i]))
    #uniport
    uniport_url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28PhaC+class+IV%29"
    download_dataset_url(uniport_url,os.path.join(prot_dir, "phac4.fasta"))
    #orthodb
    orthodb_url ="https://data.orthodb.org/current/fasta?id=122391at1224&seqtype=protein&species="
    download_dataset_url(orthodb_url,os.path.join(prot_dir, "phac1234.fasta"))
    #one file
    merge_unique_fasta(file_name_list,prot_dir, dataset123)

def glob_alig():
    alig = Protein_alignment(output_dir=prot_dir)
    alig.GlobalAlignment(file_name='phac_cupr', list_match=[2, 1, -1])

def multi_alig():
    fasta_file_phac = os.path.join(prot_dir, 'phac_cupr.fasta')
    alignment = ClustalWAlignment(file_name=fasta_file_phac)
    aligned_result = alignment.align()
    alignment.save_alignment_to_fasta(aligned_result, output_file=os.path.join(prot_dir,"alignedphac_cupr.fasta"))

def search_prot(df):
    handler = ProteinDatabaseHandlerNCBI(df)
    handler.protein_search(email= "247034@vut.cz")
    handler.download_protein(output_dir=prot_dir,file_name='phac_cupr')

def motiv_search(file='update_phac.fasta',motiv = "[GS]-X-C-X-[GA]-G"):

    stats = find_custom_motif(
        fasta_file=file,
        pattern=motiv
    )

    print(f"Pattern used: {stats['pattern_used']}")
    print(f"Total sequences: {stats['total_sequences']}")
    print(f"Sequences with motif: {stats['sequences_with_motif']}")
    print(f"Total motif occurrences: {stats['total_motifs']}")
    print(f"Average motifs per sequence: {stats['avg_motifs_per_seq']:.2f}")

    print("\nMotif positions by sequence:")
    for seq_id, motifs in stats['positions'].items():
        organism = stats['organisms'].get(seq_id, "Unknown")
        print(f"{seq_id} ({organism}): {motifs}")

if __name__ == '__main__':

    prot_dir = r'data'
    os.makedirs(prot_dir, exist_ok=True)

    df_cupr = pd.DataFrame({
        'Source': ['Cupriavidus necator', 'Caulobacter vibrioides'],
        'Gene': ['phaC', 'phaC']})

    # df = pd.DataFrame({
    #     'Source': ['Cupriavidus necator', 'Pseudomonas aeruginosa'],
    #     'Gene': ['phaC', 'phaC']})
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
    df_class_3 = df.iloc[9:12].reset_index(drop=True)
    df_class_4 = df.iloc[12:14].reset_index(drop=True)

    #create_dataset()
    #dataset123 = os.path.join(prot_dir, 'orthodb123.fasta') #net 3 iter
    #dataset123 = os.path.join(prot_dir, 'phac1.fasta') #net 3  iter
    dataset123 = os.path.join(prot_dir, 'phac1234.fasta') #rabotaet 3 iter
    #vse rabotajet
    #dataset123 = os.path.join(prot_dir,'orthodb_phac.fasta')

    dataset1 = os.path.join(prot_dir, 'filtered_phac1.fasta')
    dataset2 = os.path.join(prot_dir, 'filtered_phac2.fasta')
    dataset3 = os.path.join(prot_dir, 'filtered_phac3.fasta')
    dataset4 = os.path.join(prot_dir, 'filtered_orthodb123.fasta')

    df_classes = [df_class_1, df_class_2, df_class_3, df_class_4]
    update_fasta_files = []

    search_prot(df_class_4)
    #multi_alig()
    glob_alig()
    model = Model(filename='phac_cupr', dataset=dataset4, output_folder=prot_dir)
    # model.create_hmm_profil()
    found_seqs = model.sequential_search(iterations=3, max_sequences=10, initial_threshold=100, combine_output=True,
                                         combine_output_file=f'aligned_final4.fasta')

    # result processing
    test_fasta = os.path.join(prot_dir, f'aligned_final4.fasta')
    df_t1 = create_table(test_fasta, dataset4, output_path=prot_dir)
    updated_fasta = os.path.join(prot_dir, f'update_phac4.fasta')
    update_fasta_from_df(test_fasta, df_t1, updated_fasta)
    create_tree(os.path.join(prot_dir, f'update_phac4.fasta'))

    motiv_search(file=os.path.join(prot_dir, f'update_phac4.fasta'))

    update_fasta_files.append(updated_fasta)

    # combined_output = os.path.join(prot_dir, "final_combined_update_phac.fasta")
    #
    # with open(combined_output, "w") as outfile:
    #     for fasta_file in update_fasta_files:
    #         for record in SeqIO.parse(fasta_file, "fasta"):
    #             SeqIO.write(record, outfile, "fasta")


