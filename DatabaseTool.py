from Bio import Entrez
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import requests
import orthodb
from Bio.Seq import Seq
import json
from orthodb import OdbAPI
import re
from typing import List
#from pandas.core.interchange.dataframe_protocol import DataFrame


class ProteinDatabaseHandlerNCBI:
    """
        Class created for API requeries in NCBI protein (or others) and downloading the required sequences

        Parameters
        ----------
        protein_name : DataFrame
            Dataframe consisting of columns source and gene

        """
    def __init__(self, protein_name):
        self.protein_name = protein_name
        self.prot_dict = {}


    def protein_search(self,email:str,errors = True, database = "protein" )->dict:
        """
            Search for a protein sequence in the NCBI database.

            Parameters
            ----------
            email : str
                User email required by NCBI Entrez
            errors : bool
                If True, raises an exception on failure. Default is True
            database : str
                NCBI database to search in. Default is "protein"

            Returns
            -------
            dict
                Dictionary containing search results
        """

        Entrez.email = email
        info_dict = self.prot_dict
        protein = self.protein_name
        for i in range(np.shape(protein)[0]):
            handle = Entrez.esearch(db= database, term=protein['Source'][i] + ' ' + protein['Gene'][i], retmax=10)
            record = Entrez.read(handle)
            if errors:
                if len(record['IdList']) > 0:
                    info_dict[str(protein['Source'][i])] = record['IdList'][0]
                else:
                    raise ValueError(f"{protein['Source'][i]} is not found in NCBI {database}")

            else:
                info_dict[str(protein['Source'][i])] = record['IdList'][0] if len(record['IdList']) > 0 else 'not Found'


            print(protein['Source'][i] + ' ' + protein['Gene'][i] + f'{record}')
        return info_dict



    def download_protein(self,output_dir: str,file_name: str):
        """
            Download a protein sequence in FASTA format from NCBI.

            Parameters
            ----------
            output_dir : str
                Directory where the file will be saved.
            file_name : str
                Name of the output FASTA file.
        """

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



def download_orhodb_dataset(gen_name= "PhaC",output_file="phac.fasta"):
    """
        Download protein sequences from OrthoDB by gene name.

        Parameters
        ----------
        gen_name : str
            Target gene name to search for
        output_file : str
            Name of the output FASTA file
    """

    odb = OdbAPI()


    search_result = odb.search(gen_name)

    with open(output_file, "w") as f:
        for record in search_result.entries:
            cluster_id = record.cluster_id.id
            try:

                fasta_url = f"https://data.orthodb.org/current/fasta?id={cluster_id}"
                response = requests.get(fasta_url)

                if response.status_code == 200:
                    f.write(response.text)
                    f.write("\n")
                else:
                    print(f"Warning: Failed to download FASTA for cluster {cluster_id} (HTTP {response.status_code})")

            except Exception as e:
                print(f"Error processing cluster {cluster_id}: {str(e)}")



def merge_unique_fasta(input_files: List[str], file_path: str, output_file: str) -> None:
    """
        Merge multiple FASTA files into one, keeping only unique sequences.

        Parameters
        ----------
        input_files : List[str]
            List of input FASTA file paths.
        file_path : str
            Directory where the output file will be saved.
        output_file : str
            Name of the resulting merged FASTA file.
    """

    unique = set()
    records_to_write = []
    pattern = re.compile(r'{"pub_og_id":"(\d+at\d+)"')

    for file in input_files:
        full_path = os.path.join(file_path, file)
        for record in SeqIO.parse(full_path, "fasta"):
            header = record.description
            sequence = str(record.seq)

            match = pattern.search(header)
            seq_id = match.group(1) if match else header

            key = (seq_id, sequence)
            if key not in unique:
                unique.add(key)

                new_record = SeqRecord(
                    record.seq,
                    id=record.id,
                    description=record.description
                )
                records_to_write.append(new_record)


    with open(output_file, "w") as out_f:
        for record in records_to_write:
            SeqIO.write(record, out_f, "fasta")
            out_f.write("\n")

    print(f"Merged {len(records_to_write)} unique sequences into: {output_file}")


def download_dataset_url(url, output_file):
    """
        Download a dataset from a given URL and save it to a file.

        Parameters
        ----------
        url : str
            URL of the dataset to download.
        output_file : str
            Path to save the downloaded file.
    """

    try:
        response = requests.get(url)
        response.raise_for_status()

        with open(output_file, "w") as file:
            file.write(response.text)
        print(f"File saved: {output_file}")
    except requests.RequestException as e:
        print(f"Error: {e}")


def filter_dataset(input_fasta = 'orthodb123.fasta', output_fasta = 'filtered_orthodb123.fasta'):
    """
        Filter protein sequences from a FASTA file by specific criteria and save the result.

        Parameters
        ----------
        input_fasta : str, optional
            Path to the input FASTA file (default is 'orthodb123.fasta').
        output_fasta : str, optional
            Path to the output filtered FASTA file (default is 'filtered_orthodb123.fasta').
    """

    unique_organisms = {}
    input_count = 0

    with open(input_fasta, "r") as in_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            input_count += 1
            header = record.description
            organism = None
            if '"organism_name":"' in header:
                try:
                    organism = header.split('"organism_name":"')[1].split('"')[0]
                except IndexError:
                    continue

            if organism and organism not in unique_organisms:
                unique_organisms[organism] = record

    with open(output_fasta, "w") as out_handle:
        SeqIO.write(unique_organisms.values(), out_handle, "fasta")

    original_size = os.path.getsize(input_fasta) / (1024 * 1024)
    filtered_size = os.path.getsize(output_fasta) / (1024 * 1024)
    filtered_count = len(unique_organisms)

    print(f"seq: {input_count}")
    print(f"after filt: {filtered_count}")
    print(f"size: {original_size:.2f} MB")
    print(f"after filt: {filtered_size:.2f} MB")



# def get_organism_sequences(hits, threshold=0.01, max_sequences=100, output_file="hmm_search_results.fasta"):
#     hit_sequences = []
#     unique_names = set()
#     organism_names = []
#     sequence_count = 0
#
#     for hit in hits:
#         if sequence_count >= max_sequences or hit.evalue > threshold:
#             break
#
#         hit_description = json.loads(hit.description.decode())
#         organism_name = hit_description.get("organism_name", "Unknown organism")
#         hit_name = hit.name.decode()
#
#         if hit_name not in unique_names:
#             unique_names.add(hit_name)
#             organism_names.append(organism_name)
#
#             for domain in hit.domains:
#                 seq_record = SeqRecord(
#                     Seq(str(domain.alignment.target_sequence)),
#                     id=hit_name,
#                     description=f"{organism_name} | evalue={hit.evalue:.2g}"
#                 )
#                 hit_sequences.append(seq_record)
#                 sequence_count += 1
#                 break
#
#     if output_file:
#         SeqIO.write(hit_sequences, output_file, "fasta")
#
#     return hit_sequences




