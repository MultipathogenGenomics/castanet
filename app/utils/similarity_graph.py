import warnings
import numpy as np
import matplotlib.pyplot as plt
import biotite.sequence.align as align
import biotite.sequence.io.fasta as fasta

from app.utils.utility_fns import read_fa
from app.utils.shell_cmds import stoperr


class SimilarityGraph:
    '''Adapted from BioTite Documentation:
    doi: 10.1186/s12859-018-2367-z
    https://www.biotite-python.org/examples/gallery/sequence/pi3k_alignment.html#sphx-glr-examples-gallery-sequence-pi3k-alignment-py'''

    def __init__(self, ExpName, RefOrg, in_fname, out_fname, is_eval=True) -> None:
        if is_eval:
            fol = "evaluation"
        else:
            fol = f"consensus_data/{RefOrg}/"
        self.a = {
            "folder_stem": f"experiments/{ExpName}/",
            "ref_org": f"{RefOrg}",
            "in_fname": in_fname,
            "out_fname": f"experiments/{ExpName}/{fol}/{out_fname}"
        }
        self.bins = 200
        self.figsize = (18, 3.0)

    def construct_matrix(self, trace_code):
        '''Accepts codified sequence, generates subst matrix and calls similarity scorer'''
        matrix = align.SubstitutionMatrix.std_nucleotide_matrix()
        similarities = np.zeros(trace_code.shape)
        for i in range(similarities.shape[0]):
            for j in range(similarities.shape[1]):
                similarities[i, j] = self.get_average_normalized_similarity(
                    trace_code, matrix.score_matrix(), i, j)

        with warnings.catch_warnings():
            # Catch warnings about empty slice for gap-only parts
            warnings.simplefilter("ignore")
            similarities = self.calculate_bins(similarities, self.bins)

        return similarities

    def get_average_normalized_similarity(self, trace_code, matrix, seq_i, pos_i):
        '''Calculate avg norm similarity score'''
        code1 = trace_code[seq_i, pos_i]
        if code1 == -1:
            return np.nan
        similarities = np.zeros(trace_code.shape[0])
        for i in range(trace_code.shape[0]):
            code2 = trace_code[i, pos_i]
            if code2 == -1:
                similarities[i] = 0
            else:
                sim = matrix[code1, code2]
                # Normalize (range 0.0 - 1.0)
                min_sim = np.min(matrix[code1])
                max_sim = np.max(matrix[code1])
                sim = (sim - min_sim) / (max_sim - min_sim)
                similarities[i] = sim
        # Delete self-similarity
        similarities = np.delete(similarities, seq_i)
        return np.average(similarities)

    def calculate_bins(self, similarities, bin_count):
        '''Auto-scale bins to data'''
        edges = np.linspace(0, similarities.shape[1], bin_count, dtype=int)
        edges = np.append(edges, similarities.shape[1])
        binned_similarities = np.zeros(similarities.shape)
        for i in range(similarities.shape[0]):
            for j in range(len(edges) - 1):
                binned_similarities[i, edges[j]:edges[j+1]] = \
                    np.nanmean(similarities[i, edges[j]:edges[j+1]])
        return binned_similarities

    def draw_figure(self, seq_dict, similarities):
        '''Draw figure'''
        figure = plt.figure(figsize=self.figsize)
        ax = figure.add_subplot(111)
        heatmap = ax.pcolor(
            similarities, cmap="RdYlGn", vmin=0.0, vmax=1.0
        )
        cbar = figure.colorbar(heatmap)
        cbar.set_label("Average normalized similarity")
        ax.set_xlabel(
            f"Alignment position (bins: {round(similarities.shape[1] / self.bins)} bp)")
        ax.set_yticks(np.arange(0+0.5, len(seq_dict.keys())))
        ax.set_yticklabels(seq_dict.keys())
        figure.tight_layout()
        figure.savefig(f"{self.a['out_fname']}.png")

    def main(self):
        print(
            f"INFO: Building and graphing similarity matrix: {self.a['ref_org']}\n({self.a['out_fname']})")
        '''Load aln file'''
        seq_dict = {i[0]: i[1] for i in read_fa(f"{self.a['in_fname']}")}
        if len(seq_dict) == 0:
            stoperr(
                "There were no sequences in the input alignment file! (Similarity Graph)")
        '''Extract codes, construct sim matrix & fill similarity scores'''
        similarities = self.construct_matrix(
            align.get_codes(fasta.get_alignment(seq_dict)))
        '''Plot'''
        self.draw_figure(seq_dict, similarities)


def call_graph(seq_name, org, aln_file, out_fname, is_eval=True) -> None:
    '''Call average normalised similarity graph'''
    cls = SimilarityGraph(
        seq_name,
        org,
        aln_file,
        out_fname,
        is_eval
    )
    cls.main()


if __name__ == "__main__":
    cls = SimilarityGraph(
        "ERR10812875",
        "Paramyxoviridae_RSV",
        "consensus_alignment.aln",
        "contig_vs_ref_consensus_alignments"
    )
    cls.main()
