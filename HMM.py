import pyhmmer
import os
import shutil
import itertools
from Alignment import *
from vizualization import *
#from clustalo import clustalo


class Model:
    def __init__(self, filename: str,dataset:str, output_folder:str):

        self.output_folder = output_folder
        #self.filename = f'aligned{filename}.fasta'
        #self.hmm_profil = f'aligned{filename}.fasta.hmm'
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
        algorithm viterbi for HMM.

        :param sequence: Input sequence.
        :return: Most probable sequence of states.

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




    # def sequential_search(self, iterations=5, initial_threshold=800, max_sequences=10, combine_output=True):
    #     """
    #     Sequential search, where each successive iteration uses a model,
    #     created from the results of the previous iteration.
    #
    #     Logic:
    #     1. iteration 1: Search with the initial model (M0) → create M1
    #     2. iteration 2: Search with M1 → create M2
    #     3. iteration 3: Search with M2 → create M3.
    #     ...
    #     """
    #     found_sequences = []
    #     combined_sequences = []
    #     current_model = self.filename  # start with model (М0)
    #
    #
    #     for i in range(iterations):
    #         print(f"\n=== Iteration {i + 1}/{iterations} (model: {current_model}) ===")
    #
    #         # load model
    #         self.filename = current_model
    #         self.create_hmm_profil()
    #
    #         # search
    #         hits = list(itertools.islice(self.hmm_search(), max_sequences * 10))
    #         new_seqs = []
    #
    #         for hit in hits:
    #             try:
    #                 seq_id = hit.name.decode()
    #             except:
    #                 continue
    #
    #             if seq_id in found_sequences:
    #                 continue
    #
    #             seq = str(hit.domains[0].alignment.target_sequence)
    #             _, prob = self.viterbi(seq)
    #
    #             if prob is None:
    #                 continue
    #
    #             current_threshold = initial_threshold * (0.7 ** i)
    #             if prob >= current_threshold:
    #                 new_seqs.append((seq_id, seq))
    #                 found_sequences.append(seq_id)
    #                 print(f"Find: {seq_id} (p={prob:.2f}, threshold={current_threshold:.2f})")
    #                 if len(new_seqs) >= max_sequences:
    #                     break
    #
    #         if not new_seqs:
    #             print("No new sequences found. End the search.")
    #             break
    #
    #
    #         iter_filename = os.path.join(self.output_folder,f"iter_{i + 1}.fasta")
    #         with open(iter_filename, "w") as f:
    #             for seq_id, seq in new_seqs:
    #                 f.write(f">{seq_id}\n{seq}\n")
    #
    #         # alignmetn and create new model for next iteration
    #         clustalw = ClustalWAlignment(iter_filename)
    #         aligned_sequences = clustalw.align()
    #
    #         aligned_filename = os.path.join(self.output_folder,f"aligned_iter_{i + 1}.fasta")
    #         with open(aligned_filename, "w") as f:
    #             for seq_id, seq in aligned_sequences.items():
    #                 f.write(f">{seq_id}\n{seq}\n")
    #
    #         current_model = aligned_filename  # upload model for next iteration
    #
    #         if combine_output:
    #             combined_sequences.extend(new_seqs)
    #
    #     # final
    #     if combine_output and combined_sequences:
    #         final_filename = os.path.join(self.output_folder,"final_combined.fasta")
    #         with open(final_filename, "w") as out:
    #             for seq_id, seq in combined_sequences:
    #                 out.write(f">{seq_id}\n{seq}\n")
    #
    #         # create final model
    #         clustalw = ClustalWAlignment(final_filename)
    #         aligned_sequences = clustalw.align()
    #
    #         final_aligned = os.path.join(self.output_folder,"aligned_final.fasta")
    #         with open(final_aligned, "w") as f:
    #             for seq_id, seq in aligned_sequences.items():
    #                 f.write(f">{seq_id}\n{seq}\n")
    #
    #         self.filename = final_aligned
    #         self.create_hmm_profil()
    #         print(f"\n final model: {final_aligned}")
    #
    #     print(f"Total unique sequences found: {len(found_sequences)}")
    #     return found_sequences




    # def sequential_search(self, iterations=5, initial_threshold=800, max_sequences=10, combine_output=True,
    #                       visualization=True,switch_mode = False):
    #     """
    #     Sequential search, where each successive iteration uses a model,
    #     created from the results of the previous iteration.
    #
    #     Logic:
    #     1. iteration 1: Search with the initial model (M0) → create M1
    #     2. iteration 2: Search with M1 → create M2
    #     3. iteration 3: Search with M2 → create M3.
    #     ...
    #     """
    #     found_sequences = []
    #     combined_sequences = []
    #     current_model = self.filename  # start with model (М0)
    #     all_hits = {i: [] for i in range(iterations)}
    #     #new
    #     remaining_iterations = iterations
    #
    #     for i in range(iterations):
    #         print(f"\n=== Iteration {i + 1}/{iterations} (model: {current_model}) ===")
    #
    #         # load model
    #         self.filename = current_model
    #         self.create_hmm_profil()
    #
    #         # search
    #         hits = list(itertools.islice(self.hmm_search(), max_sequences * 10))
    #         new_seqs = []
    #         iteration_hits = []
    #
    #         for hit in hits:
    #             try:
    #                 seq_id = hit.name.decode()
    #             except:
    #                 continue
    #
    #             if seq_id in found_sequences:
    #                 continue
    #
    #             seq = str(hit.domains[0].alignment.target_sequence)
    #             _, prob = self.viterbi(seq)
    #
    #             if prob is None:
    #                 continue
    #
    #             current_threshold = initial_threshold * (0.7 ** i)
    #             if prob >= current_threshold:
    #                 new_seqs.append((seq_id, seq))
    #                 iteration_hits.append(hit)
    #                 found_sequences.append(seq_id)
    #                 print(f"Find: {seq_id} (p={prob:.2f}, threshold={current_threshold:.2f})")
    #                 if len(new_seqs) >= max_sequences:
    #                     break
    #
    #         all_hits[i] = iteration_hits[:max_sequences]
    #
    #         if not new_seqs and not switch_mode:
    #             print("No new sequences found. End the search.")
    #             break
    #
    #         #new
    #         if not new_seqs and switch_mode:
    #             remaining = iterations - (i + 1)
    #             if remaining > 0:
    #                 print(
    #                     f"No new sequences found. Switching to iterative search for {remaining_iterations - (i + 1)} iterations.")
    #
    #
    #                 iterative_hits, iterative_result = self.iterative_search(
    #                     iterations=remaining,
    #                     initial_threshold=current_threshold,
    #                     max_sequences=max_sequences,
    #                     combine_output=False,
    #                     initial_model=current_model,
    #                     vizualization=False,
    #                 )
    #
    #                 #found_sequences.extend(iterative_result)
    #                 combined_sequences.extend(iterative_result)
    #                 found_sequences.extend([seq_id for seq_id, _ in iterative_result])
    #                 for j, hits in iterative_hits.items():
    #                     all_hits[i + 1 + j] = hits
    #             else:
    #                 print("No iterations left for iterative search.")
    #                 break
    #         #
    #
    #         iter_filename = os.path.join(self.output_folder, f"iter_{i + 1}.fasta")
    #         with open(iter_filename, "w") as f:
    #             for seq_id, seq in new_seqs:
    #                 f.write(f">{seq_id}\n{seq}\n")
    #
    #         # alignment and create new model for next iteration
    #         clustalw = ClustalWAlignment(iter_filename)
    #         aligned_sequences = clustalw.align()
    #
    #         aligned_filename = os.path.join(self.output_folder, f"aligned_iter_{i + 1}.fasta")
    #         with open(aligned_filename, "w") as f:
    #             for seq_id, seq in aligned_sequences.items():
    #                 f.write(f">{seq_id}\n{seq}\n")
    #
    #         current_model = aligned_filename  # upload model for next iteration
    #
    #         if combine_output:
    #             combined_sequences.extend(new_seqs)
    #
    #     # final
    #     if combine_output and combined_sequences:
    #         final_filename = os.path.join(self.output_folder, "final_combined.fasta")
    #         with open(final_filename, "w") as out:
    #             for seq_id, seq in combined_sequences:
    #                 out.write(f">{seq_id}\n{seq}\n")
    #
    #         # create final model
    #         if len(combined_sequences) >= 2:
    #             clustalw = ClustalWAlignment(final_filename)
    #             aligned_sequences = clustalw.align()
    #
    #             final_aligned = os.path.join(self.output_folder, "aligned_final.fasta")
    #             with open(final_aligned, "w") as f:
    #                 for seq_id, seq in aligned_sequences.items():
    #                     f.write(f">{seq_id}\n{seq}\n")
    #
    #             self.filename = final_aligned
    #             self.create_hmm_profil()
    #             print(f"\n final model: {final_aligned}")
    #
    #     print(f"Total unique sequences found: {len(found_sequences)}")
    #
    #
    #     # if visualization:
    #     #     plot_hits(all_hits, output_dir=self.output_folder)
    #
    #     #new
    #     if visualization:
    #
    #         combined_hits = {}
    #         for key in sorted(all_hits.keys()):
    #             combined_hits[key] = all_hits[key]
    #         plot_hits(combined_hits, output_dir=self.output_folder)
    #     #
    #
    #     return found_sequences

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
            with open(iter_filename, "w") as f:
                for seq_id, seq in new_seqs:
                    f.write(f">{seq_id}\n{seq}\n")

            # alignment and create new model for next iteration
            clustalw = ClustalWAlignment(iter_filename)
            aligned_sequences = clustalw.align()

            aligned_filename = os.path.join(self.output_folder, f"aligned_iter_{i + 1}.fasta")
            with open(aligned_filename, "w") as f:
                for seq_id, seq in aligned_sequences.items():
                    f.write(f">{seq_id}\n{seq}\n")

            current_model = aligned_filename  # upload model for next iteration

            if combine_output:
                combined_sequences.extend(new_seqs)

        # final
        if combine_output and combined_sequences:
            final_filename = os.path.join(self.output_folder, "final_combined.fasta")
            with open(final_filename, "w") as out:
                for seq_id, seq in combined_sequences:
                    out.write(f">{seq_id}\n{seq}\n")

            # create final model
            clustalw = ClustalWAlignment(final_filename)
            aligned_sequences = clustalw.align()

            final_aligned = os.path.join(self.output_folder, combine_output_file)
            with open(final_aligned, "w") as f:
                for seq_id, seq in aligned_sequences.items():
                    f.write(f">{seq_id}\n{seq}\n")

            self.filename = final_aligned
            self.create_hmm_profil()
            print(f"\n final model: {final_aligned}")

        print(f"Total unique sequences found: {len(found_sequences)}")

        if visualization:
            plot_hits(all_hits, output_dir=self.output_folder)

        return found_sequences

    # def iterative_search(self, iterations=5, initial_threshold=800, max_sequences=10, combine_output=True,
    #                      vizualization=True,initial_model=None):
    #     """
    #     Simplified iterative search.
    #
    #     """
    #     #original_hmm = self.filename
    #     #new
    #     original_hmm = initial_model if initial_model else self.filename
    #     #
    #     found_sequences = []
    #     combined_sequences = []
    #     all_hits = {i: [] for i in range(iterations)}
    #
    #     for i in range(iterations):
    #         print(f"\n=== iteration {i + 1}/{iterations} ===")
    #
    #         self.filename = original_hmm
    #         self.create_hmm_profil()
    #
    #         hits = list(itertools.islice(self.hmm_search(), max_sequences * 10))
    #         new_seqs = []
    #         iteration_hits = []
    #
    #         for hit in hits:
    #             seq_id = hit.name.decode()
    #             if seq_id in found_sequences:
    #                 continue
    #
    #             seq = str(hit.domains[0].alignment.target_sequence)
    #             _, prob = self.viterbi(seq)
    #
    #             #if prob >= initial_threshold * (0.5 ** i):
    #             #new
    #             current_iter_threshold = initial_threshold * (
    #                         0.7 ** i)
    #             if prob >= current_iter_threshold:
    #             #
    #                 new_seqs.append((seq_id, seq))
    #                 iteration_hits.append(hit)
    #                 found_sequences.append(seq_id)
    #                 print(f"Find: {seq_id} (prob: {prob:.2f})")
    #                 if len(new_seqs) >= max_sequences:
    #                     break
    #
    #
    #         all_hits[i] = iteration_hits[:max_sequences]
    #
    #         if not new_seqs:
    #             print("No sequences were found. Skip the iteration.")
    #             continue
    #
    #         if combine_output:
    #             combined_sequences.extend(new_seqs)
    #
    #         with open(os.path.join(self.output_folder, f"iter_{i + 1}.fasta"), "w") as f:
    #             for seq_id, seq in new_seqs:
    #                 f.write(f">{seq_id}\n{seq}\n")
    #
    #     if combine_output and combined_sequences:
    #         fasta_filename = os.path.join(self.output_folder, "all_iters_combined.fasta")
    #         with open(fasta_filename, "w") as out:
    #             for seq_id, seq in combined_sequences:
    #                 out.write(f">{seq_id}\n{seq}\n")
    #
    #         clustalw = ClustalWAlignment(fasta_filename)
    #         aligned_sequences = clustalw.align()
    #         with open(fasta_filename, "w") as f:
    #             for seq_id, seq in aligned_sequences.items():
    #                 f.write(f">{seq_id}\n{seq}\n")
    #         print("\nAll iterations are combined in all_iters_combined.fasta")
    #
    #     print(f"\nFound unique orthologists: {len(found_sequences)}")
    #
    #     if vizualization:
    #         plot_hits(all_hits, output_dir=self.output_folder)
    #
    #     #return all_hits
    #     return all_hits, [(seq_id, seq) for seq_id, seq in combined_sequences]

    def iterative_search(self, iterations=5, initial_threshold=800, max_sequences=10, combine_output=True,
                         vizualization=True, initial_model=None):
        """
        Simplified iterative search with controlled output
        """
        original_hmm = initial_model if initial_model else self.filename
        found_sequences = []
        combined_sequences = []
        all_hits = {i: [] for i in range(iterations)}

        for i in range(iterations):
            print(f"\n=== iteration {i + 1}/{iterations} ===")

            self.filename = original_hmm
            self.create_hmm_profil()

            hits = list(itertools.islice(self.hmm_search(), max_sequences * 10))
            new_seqs = []
            iteration_hits = []

            for hit in hits:
                seq_id = hit.name.decode()
                if seq_id in found_sequences:
                    continue

                seq = str(hit.domains[0].alignment.target_sequence)
                _, prob = self.viterbi(seq)
                current_iter_threshold = initial_threshold * (0.7 ** i)

                if prob >= current_iter_threshold:
                    new_seqs.append((seq_id, seq))
                    iteration_hits.append(hit)
                    found_sequences.append(seq_id)
                    print(f"Find: {seq_id} (prob: {prob:.2f})")
                    if len(new_seqs) >= max_sequences:
                        break

            all_hits[i] = iteration_hits[:max_sequences]

            # Сохраняем последовательности ТОЛЬКО если явно указано
            if combine_output:  # <-- Ключевое изменение!
                if new_seqs:
                    with open(os.path.join(self.output_folder, f"iter_{i + 1}.fasta"), "w") as f:
                        for seq_id, seq in new_seqs:
                            f.write(f">{seq_id}\n{seq}\n")

                    combined_sequences.extend(new_seqs)

            # Если не сохраняем в файлы - добавляем напрямую
            else:
                combined_sequences.extend(new_seqs)

            if not new_seqs:
                print("No sequences were found. Skip the iteration.")
                continue

        # Финализация ТОЛЬКО при combine_output=True
        if combine_output and combined_sequences:
            fasta_filename = os.path.join(self.output_folder, "all_iters_combined.fasta")
            with open(fasta_filename, "w") as out:
                for seq_id, seq in combined_sequences:
                    out.write(f">{seq_id}\n{seq}\n")

            clustalw = ClustalWAlignment(fasta_filename)
            aligned_sequences = clustalw.align()
            with open(fasta_filename, "w") as f:
                for seq_id, seq in aligned_sequences.items():
                    f.write(f">{seq_id}\n{seq}\n")
            print("\nAll iterations are combined in all_iters_combined.fasta")

        print(f"\nFound unique orthologists: {len(found_sequences)}")

        if vizualization:
            plot_hits(all_hits, output_dir=self.output_folder)

        return all_hits, combined_sequences  # Возвращаем сырые данные вместо путей