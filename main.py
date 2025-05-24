import os.path
from DatabaseTool import *
from Alignment import *
from HMM import *
import pandas as pd
from vizualization import *


if __name__ == '__main__':

    prot_dir = r'data'
    os.makedirs(prot_dir, exist_ok=True)

    # df = pd.DataFrame({
    #     'Source': ['Cupriavidus necator', 'Caulobacter vibrioides'],
    #     'Gene': ['phaC', 'phaC']})

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

    #create dataset

    #orthodb
    dataset123 = os.path.join(prot_dir, 'orthodb123.fasta')
    gen_name_list = ["PhaC class I", "PhaC class II", "PhaC class III"]
    file_name_list = ["phac1.fasta", "phac2.fasta", "phac3.fasta","phac4.fasta"]
    for i in range(len(gen_name_list)):
        download_orhodb_dataset(gen_name_list[i],os.path.join(prot_dir,file_name_list[i]))


    #uniport
    uniport_url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28PhaC+class+IV%29"
    download_dataset_url(uniport_url,os.path.join(prot_dir, "phac4.fasta"))

    #one file
    merge_unique_fasta(file_name_list,prot_dir, dataset123)


    #pipline 1
    handler = ProteinDatabaseHandlerNCBI(df)
    handler.protein_search(email= "247034@vut.cz")
    handler.download_protein(output_dir=prot_dir,file_name='phac_cupr')

    #glob alig
    alig = Protein_alignment(output_dir=prot_dir)
    alig.GlobalAlignment(file_name='phac_cupr', list_match=[2, 1, -1])

    #multi alig
    fasta_file_phac = os.path.join(prot_dir, 'phac_cupr.fasta')
    alignment = ClustalWAlignment(file_name=fasta_file_phac)
    aligned_result = alignment.align()
    alignment.save_alignment_to_fasta(aligned_result, output_file=os.path.join(prot_dir,"alignedphac_cupr.fasta"))

    #model
    dataset = os.path.join(prot_dir, 'orthodb3.fasta')
    model = Model(filename='phac_cupr', dataset= dataset123, output_folder=prot_dir)
    #model.create_hmm_profil()

    #model.iterative_search(iterations=2, max_sequences=5)

    #seq seq
    found_seqs = model.sequential_search(iterations=3,max_sequences=10, initial_threshold=100,combine_output=True)

    #result processing
    test_fasta = os.path.join(prot_dir, 'aligned_final.fasta')
    # #tyt byl drugoi dataset potom sdelaj s dataset123
    df_t1 = create_table(test_fasta, dataset123,output_path=prot_dir)
    update_fasta_from_df(test_fasta,df_t1,os.path.join(prot_dir, 'update_phac.fasta'))
    create_tree(os.path.join(prot_dir, 'update_phac.fasta'))



