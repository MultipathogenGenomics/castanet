import os

from app.utils.shell_cmds import make_dir, shell
from app.utils.utility_fns import read_fa, save_fa, get_reference_org
from app.utils.similarity_graph import SimilarityGraph
from app.utils.get_genbank import DownloadGenBankFile

# For each fol in consensus_data
# If GT
# 1 Load gt flat consensus
# 2 Load gt ref-guided consensus
# 3 Load gt ref
# Run comparisons: 1-3, 2-3
# Save to csv/png - eval/gt
# Load flat consensus
#


class Evaluate:
    def __init__(self, payload) -> None:
        self.a = payload
        self.a["folder_stem"] = f"experiments/{self.a['ExpName']}/"
        self.aln_fname = f"{self.a['folder_stem']}evaluation/{self.a['GtOrg']}_consensuses.aln"
        make_dir(f"mkdir {self.a['folder_stem']}evaluation/")

    def get_consensus_seqs(self):
        names = [f">{self.a['GtOrg']}_GenBank_seq", f">{self.a['GtOrg']}_flattened_consensus",
                 f">{self.a['GtOrg']}_remapped_consensus", f">{self.a['GtOrg']}_gold_standard_consensus"]
        all_seqs = [get_reference_org(
            self.a['GtFile'], self.a['SeqName'], self.a['folder_stem'])]
        # Un-remapped flattened
        all_seqs.append(
            *read_fa(f"experiments/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_consensus_sequence.fasta"))
        all_seqs.append(
            *read_fa(f"experiments/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_remapped_consensus_sequence.fasta"))  # Remapped
        # Ground truth
        all_seqs.append(
            *read_fa(f"experiments/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_ref_adjusted_consensus.fasta"))
        all_seqs[0][0] = names[0]
        all_seqs[1][0] = names[1]
        all_seqs[2][0] = names[2]
        all_seqs[3][0] = names[3]
        with open(f"{self.a['folder_stem']}evaluation/consensus_seqs.fasta", "w") as f:
            [f.write(f"{i[0]}\n{i[1]}\n") for i in all_seqs]

    def call_alignment(self):
        shell(f"mafft --thread {os.cpu_count()} {self.a['folder_stem']}evaluation/consensus_seqs.fasta > {self.aln_fname}",
              "Mafft align consensus seqs (EVAL.PY)")

    def call_graph(self, aln_file, out_fname):
        cls = SimilarityGraph(
            self.a['SeqName'],
            self.a['GtOrg'],
            aln_file,
            out_fname
        )
        cls.main()

    def main(self):
        '''Retrieve all varieties of cons seq, graph'''
        self.get_consensus_seqs()
        self.call_alignment()
        self.call_graph(self.aln_fname, "consensus_vs_true_seq")

        '''Call graph on contig consensuses'''
        self.call_graph(
            f"experiments/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_consensus_alignment.aln", "contig_vs_ref_consensus_alignments")


if __name__ == "__main__":
    cls = Evaluate(
        {
            "ExpName": "ERR10812875",
            "SeqName": "ERR10812875",
            "GtOrg": "Paramyxoviridae_RSV",
            "GtFile": "data/rsv_set_metadata.csv"
        }
    )
    cls.main()
