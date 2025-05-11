import os.path
from DatabaseTool import *
from Alignment import *
from HMM import *
import pandas as pd
from vizualization import *


if __name__ == '__main__':

    prot_dir = r'data'
    os.makedirs(prot_dir, exist_ok=True)

    df = pd.DataFrame({
        'Source': ['Cupriavidus necator', 'Caulobacter vibrioides'],
        'Gene': ['phaC', 'phaC']})

    #pipline 1
    handler = ProteinDatabaseHandlerNCBI(df)
    handler.protein_search(email= "247034@vut.cz")
    handler.download_protein(output_dir=prot_dir,file_name='phac_cupr')

    #first alig
    alig = Protein_alignment(output_dir=prot_dir)
    alig.GlobalAlignment(file_name='phac_cupr', list_match=[2, 1, -1])

    #model
    dataset = os.path.join(prot_dir, 'orthodb_phac.fasta')
    model = Model(filename='phac_cupr', dataset= dataset, output_folder=prot_dir)
    #model.create_hmm_profil()


    #model.iterative_search(iterations=2, max_sequences=5)
    #model.upload_search(iterations=3,max_sequences=10,min_sequences_for_update=2)

    #seq seqar
    found_seqs = model.sequential_search(iterations=5,max_sequences=20, initial_threshold=500,combine_output=True)

    test_fasta = os.path.join(prot_dir, 'aligned_final.fasta')
    df = create_table(test_fasta, dataset)

    update_fasta_from_df(test_fasta,df)
    create_tree(test_fasta)



