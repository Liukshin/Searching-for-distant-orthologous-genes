import os.path
from DatabaseTool import *
from Alignment import *
from HMM import *
import pandas as pd
from vizualization import *
from config import (
    prot_dir,
    df, df_class_1, df_class_2, df_class_3, df_class_4,
    dataset123, dataset1, dataset2, dataset3, dataset4
)




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
    #filtr datset
    input_fasta = os.path.join(prot_dir, 'orthodb123.fasta')
    output_fasta = os.path.join(prot_dir, 'filtered_orthodb123.fasta')
    filter_dataset(input_fasta, output_fasta)

def glob_alig():
    fasta_file_phac = os.path.join(prot_dir, 'phac_cupr.fasta')
    seqs = list(SeqIO.parse(fasta_file_phac, 'fasta'))
    alignment = GlobalAlignment(seqs[0],seqs[1],match=2,mismatch=-1,gap= -2)
    sek1,sek2 = alignment.align()
    record1 = SeqRecord(Seq(sek1), id="seq1")
    record2 = SeqRecord(Seq(sek2), id="seq2")
    output_path = os.path.join(prot_dir,"alignedphac_cupr.fasta")

    SeqIO.write([record1, record2], output_path, "fasta")


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



def main(df_class:DataFrame ,number_pipline:int,multi_alignment= True,threshold = 0.35):

    search_prot(df_class)
    if multi_alignment:
        multi_alig()
    else:
        glob_alig()
    model = Model(filename='phac_cupr', dataset=dataset4, output_folder=prot_dir)
    # model.create_hmm_profil()
    found_seqs = model.sequential_search(iterations=3, max_sequences=10, initial_threshold=100, combine_output=True,
                                         combine_output_file=f'aligned_final{number_pipline}.fasta')

    # result processing
    test_fasta = os.path.join(prot_dir, f'aligned_final{number_pipline}.fasta')
    df_t1 = create_table(file_result=test_fasta, db_file=dataset4, table_name=f"table_results_phac{number_pipline}.csv",
                         output_path=prot_dir)
    updated_fasta = os.path.join(prot_dir, f'update_phac{number_pipline}.fasta')
    update_fasta_from_df(test_fasta, df_t1, updated_fasta)
    # simple tree
    #create_tree(os.path.join(prot_dir, f'update_phac{number_pipline}.fasta'))
    motiv_search(file=os.path.join(prot_dir, f'update_phac{number_pipline}.fasta'))

    # colored tree

    fasta_file = os.path.join(prot_dir,f"update_phac{number_pipline}.fasta")
    alignment = ClustalWAlignment(fasta_file)
    newick_str = alignment.neighbor_joining(fyl_tree_viz=True)

    clust_dict = clusters_tree(newick_str, fasta_file, output_image=f"colored_clusters{number_pipline}.png", path=prot_dir,
                               cluster_threshold=threshold)
    print(clust_dict)



if __name__ == '__main__':
    df_classes = [df_class_1, df_class_2, df_class_3, df_class_4]


    #PhaC I
    #main(df_classes[0],number_pipline=1)
    # #PhaC II
    # main(df_classes[1],number_pipline=2)
    # #PhaC III
    #main(df_classes[2],number_pipline=3,multi_alignment=False,threshold = 0.6)
    # #PhaC IV
    main(df_classes[3],number_pipline=4,multi_alignment=False)



