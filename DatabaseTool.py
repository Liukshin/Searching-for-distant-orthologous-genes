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




class ProteinDatabaseHandlerNCBI:
    def __init__(self, protein_name):
        self.protein_name = protein_name
        self.prot_dict = {}


    def protein_search(self,email:str,errors = True, database = "protein" ):
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



def download_orhodb_dataset(gen_name= "PhaC",output_file="phac.fasta"):
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

# def merge_unique_fasta(input_files: List[str],file_path, output_file: str) -> None:
#
#     unique_sequences = set()
#     pattern = re.compile(r'{"pub_og_id":"(\d+at\d+)"')
#
#     with open(output_file, 'w') as out_f:
#         for file in input_files:
#             with open(os.path.join(file_path,file), 'r') as in_f:
#                 current_header = None
#                 current_sequence = []
#
#                 for line in in_f:
#                     line = line.strip()
#                     if line.startswith('>'):
#
#                         if current_header and current_sequence:
#                             seq_id_match = pattern.search(current_header)
#                             seq_id = seq_id_match.group(1) if seq_id_match else current_header
#                             sequence = ''.join(current_sequence)
#
#                             if (seq_id, sequence) not in unique_sequences:
#                                 unique_sequences.add((seq_id, sequence))
#                                 out_f.write(f"{current_header}\n{sequence}\n")
#
#
#                         current_header = line
#                         current_sequence = []
#                     else:
#                         current_sequence.append(line)
#
#
#                 if current_header and current_sequence:
#                     seq_id_match = pattern.search(current_header)
#                     seq_id = seq_id_match.group(1) if seq_id_match else current_header
#                     sequence = ''.join(current_sequence)
#
#                     if (seq_id, sequence) not in unique_sequences:
#                         unique_sequences.add((seq_id, sequence))
#                         out_f.write(f"{current_header}\n{sequence}\n")
#
#     print(f"merge {len(unique_sequences)} unique sequences in file {output_file}")





def merge_unique_fasta(input_files: List[str], file_path: str, output_file: str) -> None:
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
    try:
        response = requests.get(url)
        response.raise_for_status()

        with open(output_file, "w") as file:
            file.write(response.text)
        print(f"File saved: {output_file}")
    except requests.RequestException as e:
        print(f"Error: {e}")


def get_organism_sequences(hits, threshold=0.01, max_sequences=100, output_file="hmm_search_results.fasta"):
    hit_sequences = []
    unique_names = set()
    organism_names = []
    sequence_count = 0

    for hit in hits:
        if sequence_count >= max_sequences or hit.evalue > threshold:
            break

        hit_description = json.loads(hit.description.decode())
        organism_name = hit_description.get("organism_name", "Unknown organism")
        hit_name = hit.name.decode()

        if hit_name not in unique_names:
            unique_names.add(hit_name)
            organism_names.append(organism_name)

            for domain in hit.domains:
                seq_record = SeqRecord(
                    Seq(str(domain.alignment.target_sequence)),
                    id=hit_name,
                    description=f"{organism_name} | evalue={hit.evalue:.2g}"
                )
                hit_sequences.append(seq_record)
                sequence_count += 1
                break

    if output_file:
        SeqIO.write(hit_sequences, output_file, "fasta")

    return hit_sequences




