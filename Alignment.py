from io import StringIO
from Bio import Phylo
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import os



class GlobalAlignment:
    """Needleman-Wunsch algorthm

    Parameters
    ----------
    seq1 : str
        First input sequence
    seq2 : str
        Second input sequence
    match : int
        Score for a match
    mismatch : int
        Penalty for a mismatch
    gap : int
        Penalty for a gap

    """

    def __init__(self, seq1: str, seq2: str, match: int, mismatch: int, gap: int):
        self.seq1 = seq1
        self.seq2 = seq2
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

    def align(self):
        """Global alignment
            Aligns only two sequences
        """
        seq1, seq2 = self.seq1, self.seq2
        len1, len2 = len(seq1), len(seq2)
        arr_s = np.zeros((len1 + 1, len2 + 1))

        for i in range(1, len2 + 1):
            arr_s[0][i] = arr_s[0][i - 1] - abs(self.gap)
        for i in range(1, len1 + 1):
            arr_s[i][0] = arr_s[i - 1][0] - abs(self.gap)

        for i in range(1, len1 + 1):
            for j in range(1, len2 + 1):
                arr_s[i][j] = max(
                    arr_s[i - 1][j] - abs(self.gap),
                    arr_s[i][j - 1] - abs(self.gap),
                    arr_s[i - 1][j - 1] + (self.match if seq1[i - 1] == seq2[j - 1] else self.mismatch)
                )

        aligned1, aligned2 = "", ""
        i, j = len1, len2

        while i > 0 or j > 0:
            if i > 0 and arr_s[i][j] == arr_s[i - 1][j] - abs(self.gap):
                aligned1 = seq1[i - 1] + aligned1
                aligned2 = "-" + aligned2
                i -= 1
            elif j > 0 and arr_s[i][j] == arr_s[i][j - 1] - abs(self.gap):
                aligned1 = "-" + aligned1
                aligned2 = seq2[j - 1] + aligned2
                j -= 1
            else:
                aligned1 = seq1[i - 1] + aligned1
                aligned2 = seq2[j - 1] + aligned2
                i -= 1
                j -= 1

        return aligned1, aligned2



class MultipleAlignment:
    """Basic class for MSA"""

    def __init__(self, file_name: str):
        self.file_name = file_name
        self.sequences = self.read_fasta()

    def read_fasta(self):

        sequences = list(SeqIO.parse(self.file_name, 'fasta'))
        return sequences

    def align(self):

        raise NotImplementedError("This method must be overridden in the subclass")



class ClustalWAlignment(MultipleAlignment):

    def __init__(self, file_name: str, match: int = 2, mismatch: int = -1, gap: int = -2):
        super().__init__(file_name)
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        #kal dorabotka
        self.original_ids = [seq.id for seq in self.sequences]

    def compute_distance_matrix(self):

        num_seqs = len(self.sequences)
        distance_matrix = np.zeros((num_seqs, num_seqs))

        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                seq1, seq2 = str(self.sequences[i].seq), str(self.sequences[j].seq)
                aligned_seq1, aligned_seq2 = GlobalAlignment(seq1, seq2,match=2,mismatch=-1,gap=-2).align()

                matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
                max_length = max(len(aligned_seq1), len(aligned_seq2))
                identity = matches / max_length if max_length > 0 else 0

                distance = 1 - identity
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance

        return distance_matrix

    def _create_id_mapping(self):

        self.id_to_safe = {}
        self.safe_to_id = {}

        for i, seq in enumerate(self.sequences):
            original_id = seq.id
            safe_id = f"seq_{i}"
            self.id_to_safe[original_id] = safe_id
            self.safe_to_id[safe_id] = original_id


    def neighbor_joining(self,fyl_tree_viz = False):
        """Builds a phylogenetic tree using the neighbor-joining algorithm

            Parameter
            ----------
            fyl_tree_viz : bool
            responsible for adding original names of sequences (needed when creating a phylogenetic tree in the visualization module)

        """

        distance_matrix = self.compute_distance_matrix()
        if fyl_tree_viz:
            labels = [seq.id for seq in self.sequences]
        else:

            labels = [f"seq_{i}" for i in range(len(self.sequences))]
        #

        D = np.array(distance_matrix, dtype=float)
        clusters = {label: label for label in labels}

        while len(labels) > 2:
            n = len(D)

            # 1. Compute Q matrix
            total_distances = np.sum(D, axis=1)
            Q = (n - 2) * D - total_distances[:, np.newaxis] - total_distances[np.newaxis, :]
            np.fill_diagonal(Q, np.inf)

            #2. Find pair with minimal Q[i][j]
            i, j = np.unravel_index(np.argmin(Q), Q.shape)
            if i == j:  # Safety check
                continue


            if i > j:
                i, j = j, i

            #3 Calculate branch lengths for the new node
            delta = (total_distances[i] - total_distances[j]) / (n - 2)
            limb_i = 0.5 * (D[i, j] + delta)
            limb_j = 0.5 * (D[i, j] - delta)

            #4. Create new cluster in Newick format
            new_cluster = f"({clusters[labels[i]]}:{limb_i:.3f}, {clusters[labels[j]]}:{limb_j:.3f})"
            new_label = f"Cluster_{len(labels)}"
            clusters[new_label] = new_cluster

            #5 Update the distance matrix
            new_distances = 0.5 * (D[i, :] + D[j, :] - D[i, j])
            new_distances = np.delete(new_distances, [i, j])


            D = np.delete(D, [j, i], axis=0)
            D = np.delete(D, [j, i], axis=1)


            D = np.vstack([D, new_distances])
            new_distances = np.append(new_distances, 0.0)
            D = np.column_stack([D, new_distances])

            #6. Update labels list
            labels.pop(j)
            labels.pop(i)
            labels.append(new_label)

        # 7. Connect the last two nodes
        if len(labels) == 2:
            tree = f"({clusters[labels[0]]}:{D[0, 1] / 2:.3f}, {clusters[labels[1]]}:{D[0, 1] / 2:.3f})"
        else:
            tree = clusters[labels[0]]

        return tree


    def read_newick(self):
        newick_str = self.neighbor_joining()

        
        tree = Phylo.read(StringIO(newick_str), "newick")

        terminals = list(tree.get_terminals())
        if len(terminals) != len(self.original_ids):
            raise ValueError("Number of terminals doesn't match number of sequences")

        for terminal, original_id in zip(terminals, self.original_ids):
            terminal.name = original_id

        return tree


    def progressive_alignment(self):
        """Progressive alignment with guaranteed equal lengths"""

        tree = self.read_newick()
        sequences_dict = {seq.id: str(seq.seq) for seq in self.sequences}

        # 1. Initialization - store both original and aligned sequences
        for clade in tree.find_clades():
            if clade.is_terminal():
                clade.original_seq = sequences_dict[clade.name]
                clade.aligned_seq = sequences_dict[clade.name]

        #2 post-order tree traversal
        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal():
                if len(clade.clades) != 2:
                    raise ValueError("Tree must be strictly binary")

                left, right = clade.clades

                #3 Align representatives (without gaps)
                aligner = GlobalAlignment(
                    left.aligned_seq.replace('-', ''),
                    right.aligned_seq.replace('-', ''),
                    match=self.match,
                    mismatch=self.mismatch,
                    gap=self.gap
                )
                aligned_left, aligned_right = aligner.align()

                # 4. Find all gap positions in the alignment
                gap_positions = set()
                for i, (a, b) in enumerate(zip(aligned_left, aligned_right)):
                    if a == '-' or b == '-':
                        gap_positions.add(i)

                #5 Apply gaps to all sequences in subtrees
                def apply_global_gaps(node, gaps):
                    if node.is_terminal():

                        seq_list = list(node.original_seq)
                        for pos in sorted(gaps, reverse=True):
                            seq_list.insert(pos, '-')
                        node.aligned_seq = ''.join(seq_list)
                    else:
                        for child in node.clades:
                            apply_global_gaps(child, gaps)

                apply_global_gaps(left, [i for i, char in enumerate(aligned_left) if char == '-'])
                apply_global_gaps(right, [i for i, char in enumerate(aligned_right) if char == '-'])

                # 6. Store the merged alignment
                clade.aligned_seq = aligned_left

        # 7 compile final alignment
        msa = {}
        for clade in tree.get_terminals():
            msa[clade.name] = clade.aligned_seq

        # 8. Length validation
        lengths = {k: len(v) for k, v in msa.items()}
        if len(set(lengths.values())) != 1:


            max_len = max(lengths.values())
            for name in msa:
                msa[name] = msa[name].ljust(max_len, '-')

        return msa

    def align(self):
        """Main ClustalW algorithm"""
        print("Step 1: Computing distance matrix")
        distance_matrix = self.compute_distance_matrix()

        print("Step 2: Building guide tree")
        tree = self.neighbor_joining()
        print(f"Guide tree: {tree}")

        print("Step 3: Progressive alignment")
        aligned_result = self.progressive_alignment()
        return aligned_result


    def save_alignment_to_fasta(self,alignment_seq, headers=False, output_file="alignment_result.fasta"):

        records = []
        for seq_id, aligned_seq in alignment_seq.items():
            if headers:
                record = SeqRecord(Seq(aligned_seq),
                                   id=seq_id,
                                   description=f"{headers[seq_id]} Aligned length: {len(aligned_seq)}")
            else:
                record = SeqRecord(Seq(aligned_seq), id=seq_id)

            records.append(record)

        with open(output_file, "w") as handle:
            SeqIO.write(records, handle, "fasta")

        print(f"Alignment saved to {os.path.abspath(output_file)}")
        return output_file