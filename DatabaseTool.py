from Bio import Entrez
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import requests
import orthodb





class ProteinDatabaseHandlerNCBI:
    def __init__(self, protein_name):
        self.protein_name = protein_name
        self.prot_dict = {}


    def protein_search(self,email:str,errors = True):
        Entrez.email = email
        info_dict = self.prot_dict
        protein = self.protein_name
        for i in range(np.shape(protein)[0]):
            handle = Entrez.esearch(db="protein", term=protein['Source'][i] + ' ' + protein['Gene'][i], retmax=10)
            record = Entrez.read(handle)
            if errors:
                if len(record['IdList']) > 0:
                    info_dict[str(protein['Source'][i])] = record['IdList'][0]
                else:
                    raise ValueError(f"{protein['Source'][i]} is not found in NCBI GenBank")

            else:
                info_dict[str(protein['Source'][i])] = record['IdList'][0] if len(record['IdList']) > 0 else 'not Found'


            print(protein['Source'][i] + ' ' + protein['Gene'][i] + f'{record}')
        return info_dict

    def get_names(self,email:str, ids: list):
        Entrez.email = email
        id_to_name = {}

        for pid in ids:
            try:
                handle = Entrez.efetch(db="protein", id=pid, rettype="gb", retmode="text")
                record = SeqIO.read(handle, "gb")
                id_to_name[pid] = record.description
            except Exception as e:
                print(f"Error getting name for {pid}: {e}")
                id_to_name[pid] = pid  # fallback

        return id_to_name



    def download_protein(self,output_dir: str,file_name: str):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        seq_arr = []
        info_dict = self.prot_dict
        for name_org in info_dict:
            handle = Entrez.efetch(db="protein", id=info_dict[name_org], rettype="fasta", retmode="text")
            record = SeqIO.read(handle, 'fasta')
            record2 = SeqRecord(record.seq, id=info_dict[name_org])
            seq_arr.append(record2)

        output_path = os.path.join(output_dir, f'{file_name}.fasta')

        SeqIO.write(seq_arr, output_path, 'fasta')


        print("download.....\033[1;32;40m success  \n")






def download_orthodb_fasta(gene_name: str, output_file: str = "orthologs.fasta", output_dir: str = "./fasta_results"):
    """
    Uses the official OrthoDB Python package to fetch all orthologs of a gene and save them as a FASTA file.

    Args:
        gene_name (str): Name of the gene (e.g., 'PhaC').
        output_file (str): Filename for the resulting FASTA file.
        output_dir (str): Directory to store the output file.
    """
    os.makedirs(output_dir, exist_ok=True)
    full_path = os.path.join(output_dir, output_file)

    # Connect to OrthoDB API
    api = orthodb.OdbAPI()

    search_results = api.search(gene_name)
    if not search_results or len(search_results) == 0:
        raise Exception(f"No ortholog groups found for gene: {gene_name}")

    fasta_data = ""

    print(f"Found {len(search_results)} ortholog group(s) for gene '{gene_name}'")

    for group in search_results:
        group_id = group.get("id")
        if not group_id:
            continue

        try:
            # Get protein FASTA for this group
            fasta = api.fasta(group_id)
            fasta_data += fasta + "\n"
            print(f"Group {group_id} added.")
        except Exception as e:
            print(f"Error with group {group_id}: {e}")

    if not fasta_data:
        raise Exception("No FASTA data downloaded.")

    with open(full_path, "w") as f:
        f.write(fasta_data)

    print(f"All FASTA sequences saved to {full_path}")

