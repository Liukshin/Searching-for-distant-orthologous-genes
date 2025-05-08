from DatabaseTool import *
from Alignment import *
from HMM import *
import pandas as pd
from vizualization import *


if __name__ == '__main__':

    prot_dir = r'data'

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
    model = Model(filename='phac_cupr', dataset='orthodb_phac.fasta')
    #model.create_hmm_profil()

    # Итеративный поиск
    #model.iterative_search(iterations=2, max_sequences=5)
    #model.upload_search(iterations=3,max_sequences=10,min_sequences_for_update=2)
    found_seqs = model.sequential_search(iterations=3, initial_threshold=500)

    #handler = ProteinDatabaseHandlerNCBI(None)
    #create_tree("all_iters_combined.fasta", db_handler=handler)



