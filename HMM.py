import pyhmmer
import os
import shutil
import itertools
from Alignment import *
from vizualization import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


class Model:
    def __init__(self, filename: str,dataset:str, output_folder:str):

        self.output_folder = output_folder
        self.filename = os.path.join(self.output_folder, f'aligned{filename}.fasta')
        self.hmm_profil = os.path.join(self.output_folder, f'aligned{filename}.fasta.hmm')
        self.dataset = dataset
        self.hmm = None


    def create_hmm_profil(self):
        alphabet = pyhmmer.easel.Alphabet.amino()
        with pyhmmer.easel.MSAFile(self.filename, digital=True, alphabet=alphabet) as msa_file:
            msa = msa_file.read()
            msa.name = b"My_HMM"
            builder = pyhmmer.plan7.Builder(alphabet)
            background = pyhmmer.plan7.Background(alphabet)
            self.hmm, _, _ = builder.build_msa(msa, background)
            with open(f"{self.filename}.hmm", "wb") as hmm_file:
                self.hmm.write(hmm_file)


    def hmm_search(self)->pyhmmer.plan7.TopHits:
        with pyhmmer.plan7.HMMFile(self.hmm_profil) as hmm_file:
            hmm = hmm_file.read()

        with pyhmmer.easel.SequenceFile(self.dataset, digital=True) as seq_file:
            sequences = seq_file.read_block()

        pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
        hits = pipeline.search_hmm(hmm, sequences)
        return hits





    def viterbi(self, sequence):
        """
         Viterbi algorithm to find the most probable path of hidden states in a profile HMM.

         Parameters
         ----------
         sequence : str
             Input amino acid sequence.

         Returns
         -------
         best_path : list[str]
             Most probable sequence of HMM states (by index).
         best_path_prob : float
             Log-probability score of the best path.

         """

        alphabet = pyhmmer.easel.Alphabet.amino()

        if self.hmm is None:
            raise ValueError("The HMM model is not loaded. First call create_hmm_profil().")

        num_states = self.hmm.M  # numbers of state Match
        seq_len = len(sequence)

        # Table initialization
        viterbi = np.zeros((num_states, seq_len))
        #viterbi = np.full((num_states, seq_len), -np.inf)
        backpointer = np.zeros((num_states, seq_len), dtype=int)

        symbol = sequence[0]
        valid_states = np.arange(num_states)

        if symbol == "-":
            viterbi[:, 0] = ([self.hmm.transition_probabilities[state, 5] for state in range(num_states)])  # d->m
        else:
            symbol_index = int(alphabet.encode(symbol)[0])
            if symbol_index >= 20:
                print(f"Incorrect symbol {symbol}, index {symbol_index} out of range. Skip symbol {symbol}.")
                return  # skip the current sequence
            trans_prob =np.array([self.hmm.transition_probabilities[state, 0] for state in range(num_states)])
            match_prob = np.array([self.hmm.match_emissions[state, symbol_index] for state in range(num_states)])# m->m
            viterbi[:, 0] = (trans_prob + match_prob)


        for t in range(1, seq_len):
            symbol = sequence[t]

            if symbol != "-":
                symbol_index = alphabet.encode(symbol)[0]
                if symbol_index >= 20:
                    print(f"Incorrect symbol {symbol}, index {symbol_index} out of range. Skip symbol {symbol}.")
                    continue  # skip symbol and continue
                emission_probs = np.array([self.hmm.match_emissions[state, symbol_index] for state in range(num_states)])
            else:
                emission_probs = np.zeros(num_states)  # gap (-)

            # transitions
            transition_probs = np.array([self.hmm.transition_probabilities[state, 0]for state in range(num_states)])  # m->m по умолчанию
            transition_probs[symbol == "-"] = np.array([self.hmm.transition_probabilities[state, 6]for state in range(num_states)])  # d->d для '-'

            # vectorized update
            viterbi[:, t] = np.max(viterbi[:, t - 1][:, np.newaxis] + (transition_probs) + (emission_probs), axis=0)
            backpointer[:, t] = np.argmax(viterbi[:, t - 1][:, np.newaxis] + np.array(transition_probs) + np.array(emission_probs), axis=0)
        best_path_prob = np.nan  # start with none value
        best_last_state = 0

        # seqrch max prob in the last colmn  viterbi
        for state in range(num_states):
            if np.isnan(best_path_prob) or viterbi[state, seq_len - 1] > best_path_prob:
                best_path_prob = viterbi[state, seq_len - 1]
                best_last_state = state

        # if the best prob= -∞ , its mistake
        if np.isnan(best_path_prob):
            print("Error: All stats have -∞ probability in the last column Viterbi.")
            return [], -1  #  return emty arr and -1 for prob

        # path restoration
        best_path = [best_last_state]
        for t in range(seq_len - 1, 0, -1):
            best_last_state = backpointer[best_last_state, t]
            best_path.insert(0, best_last_state)

        return best_path, best_path_prob




    def sequential_search(self, iterations=5, initial_threshold=800, max_sequences=10, combine_output=True,
                          visualization=True,combine_output_file = "aligned_final.fasta"):
        """
        Sequential search, where each successive iteration uses a model,
        created from the results of the previous iteration.

        Logic:
        1. iteration 1: Search with the initial model (M0) → create M1
        2. iteration 2: Search with M1 → create M2
        3. iteration 3: Search with M2 → create M3.
        ...
        Parameters
        ----------
        iterations : int
            Number of search iterations to perform.
        initial_threshold : int
            Score threshold for filtering hits in the first iteration.
        max_sequences : int
            Maximum number of sequences to retain per iteration.
        combine_output : bool
            Whether to merge results from all iterations into a single file.
        visualization : bool
            Whether to generate visualization for the final output.
        combine_output_file : str
            Output file name for combined aligned sequences in FASTA format.

        Returns
        -------
        found_sequences:list[str]
            List of sequences found during the search process.

        """
        found_sequences = []
        combined_sequences = []
        current_model = self.filename  # start with model (М0)
        all_hits = {i: [] for i in range(iterations)}

        for i in range(iterations):
            print(f"\n=== Iteration {i + 1}/{iterations} (model: {current_model}) ===")

            # load model
            self.filename = current_model
            self.create_hmm_profil()

            # search
            hits = list(itertools.islice(self.hmm_search(), max_sequences * 10))
            new_seqs = []
            iteration_hits = []

            for hit in hits:
                try:
                    seq_id = hit.name.decode()
                except:
                    continue

                if seq_id in found_sequences:
                    continue

                seq = str(hit.domains[0].alignment.target_sequence)
                _, prob = self.viterbi(seq)

                if prob is None:
                    continue

                current_threshold = initial_threshold * (0.7 ** i)
                if prob >= current_threshold:
                    new_seqs.append((seq_id, seq))
                    iteration_hits.append(hit)
                    found_sequences.append(seq_id)
                    print(f"Find: {seq_id} (p={prob:.2f}, threshold={current_threshold:.2f})")
                    if len(new_seqs) >= max_sequences:
                        break

            all_hits[i] = iteration_hits[:max_sequences]

            if not new_seqs:
                print("No new sequences found. End the search.")
                break

            iter_filename = os.path.join(self.output_folder, f"iter_{i + 1}.fasta")
            iter_records = [SeqRecord(Seq(seq), id=seq_id, description="") for seq_id, seq in new_seqs]
            with open(iter_filename, "w") as f:
                SeqIO.write(iter_records, f, "fasta")

            # alignment and create new model for next iteration
            clustalw = ClustalWAlignment(iter_filename)
            aligned_sequences = clustalw.align()

            aligned_filename = os.path.join(self.output_folder, f"aligned_iter_{i + 1}.fasta")
            aligned_records = [SeqRecord(Seq(seq), id=seq_id, description="") for seq_id, seq in
                               aligned_sequences.items()]
            with open(aligned_filename, "w") as f:
                SeqIO.write(aligned_records, f, "fasta")

            current_model = aligned_filename

            if combine_output:
                combined_sequences.extend(new_seqs)

        # final
        if combine_output and combined_sequences:
            final_filename = os.path.join(self.output_folder, "final_combined.fasta")
            final_records = [SeqRecord(Seq(seq), id=seq_id, description="") for seq_id, seq in combined_sequences]
            with open(final_filename, "w") as out:
                SeqIO.write(final_records, out, "fasta")


            clustalw = ClustalWAlignment(final_filename)
            aligned_sequences = clustalw.align()

            final_aligned = os.path.join(self.output_folder, combine_output_file)
            aligned_final_records = [SeqRecord(Seq(seq), id=seq_id, description="") for seq_id, seq in
                                     aligned_sequences.items()]
            with open(final_aligned, "w") as f:
                SeqIO.write(aligned_final_records, f, "fasta")


            self.filename = final_aligned
            self.create_hmm_profil()
            print(f"\n final model: {final_aligned}")

        print(f"Total unique sequences found: {len(found_sequences)}")

        if visualization:
            plot_hits(all_hits, output_dir=self.output_folder)

        return found_sequences


