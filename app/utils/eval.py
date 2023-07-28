import os
import time
import pandas as pd

from app.utils.shell_cmds import make_dir, shell
from app.utils.utility_fns import read_fa, save_fa, get_reference_org
from app.utils.similarity_graph import SimilarityGraph
from app.utils.report_gen import GenerateReport


class Evaluate:
    def __init__(self, payload) -> None:
        self.a = payload
        if not "StartTime" in self.a.keys():
            self.a["StartTime"] = time.time()
        self.a["folder_stem"] = f"experiments/{self.a['ExpName']}/"
        self.aln_fname = f"{self.a['folder_stem']}evaluation/{self.a['GtOrg']}_consensuses.aln"
        self.seq_names = [f">{self.a['GtOrg']}_GenBank_seq", f">{self.a['GtOrg']}_flattened_consensus",
                          f">{self.a['GtOrg']}_remapped_consensus", f">{self.a['GtOrg']}_gold_standard_consensus"]
        self.contig_graph_fname = f"{self.a['folder_stem']}evaluation/contig_vs_ref_consensus_alignments.png"
        self.consen_graph_fname = f"{self.a['folder_stem']}evaluation/consensus_vs_true_seq.png"
        make_dir(
            f"mkdir {self.a['folder_stem']}evaluation/ {self.a['folder_stem']}evaluation/seqs/")

    def get_consensus_seqs(self) -> None:
        '''Read in reference and consensus seqs, rename seqs, save versions to eval folder for stats and graphing'''
        all_seqs = [get_reference_org(
            self.a['GtFile'], self.a['SeqName'], self.a['folder_stem'])]
        all_seqs.append(  # Un-remapped flattened
            *read_fa(f"experiments/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_consensus_sequence.fasta"))
        all_seqs.append(  # Remapped
            *read_fa(f"experiments/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_remapped_consensus_sequence.fasta"))
        all_seqs.append(  # Ground truth
            *read_fa(f"experiments/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_ref_adjusted_consensus.fasta"))
        '''Rename seqs, save individual fa, build mash sketch and combined fa (for aln)'''
        for i in range(len(all_seqs)):
            all_seqs[i][0] = self.seq_names[i]
            save_fa(f"{self.a['folder_stem']}evaluation/seqs/{self.seq_names[i].replace('>','')}.fa",
                    f">{all_seqs[i][0]}\n{all_seqs[i][1]}\n")
            shell(
                f"mash sketch {self.a['folder_stem']}evaluation/seqs/{self.seq_names[i].replace('>','')}.fa")

        with open(f"{self.a['folder_stem']}evaluation/consensus_seqs.fasta", "w") as f:
            [f.write(f"{i[0]}\n{i[1]}\n") for i in all_seqs]

    def call_alignment(self) -> None:
        '''Build MAFFT MSA of all consensuses and "true" sequence'''
        shell(f"mafft --thread {os.cpu_count()} {self.a['folder_stem']}evaluation/consensus_seqs.fasta > {self.aln_fname}",
              "Mafft align consensus seqs (EVAL.PY)")

    def call_graph(self, aln_file, out_fname) -> None:
        '''Call average normalised similarity graph'''
        cls = SimilarityGraph(
            self.a['SeqName'],
            self.a['GtOrg'],
            aln_file,
            out_fname
        )
        cls.main()

    def call_mash_dist(self) -> None:
        '''Call mash distances on (1) all consensuses vs original viral genome, (2) flat/remapped cons vs gold standard cons'''
        shell(f"mash dist {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_GenBank_seq.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_flattened_consensus.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_remapped_consensus.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_gold_standard_consensus.fa.msh > {self.a['folder_stem']}evaluation/mash_results_cons_vs_true_genome.csv")
        shell(f"mash dist {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_gold_standard_consensus.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_flattened_consensus.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_remapped_consensus.fa.msh > {self.a['folder_stem']}evaluation/mash_results_cons_vs_gold_standard.csv")

    def collate_stats(self) -> pd.DataFrame:
        '''Grab MASH output, put in dataframe, save'''
        df = pd.read_csv(f"{self.a['folder_stem']}evaluation/mash_results_cons_vs_true_genome.csv",
                         delimiter="\t", header=None, names=["ref", "sample", "mash_dist", "p", "matching_hashes"])
        df = pd.concat([df, pd.read_csv(f"{self.a['folder_stem']}evaluation/mash_results_cons_vs_gold_standard.csv",
                       delimiter="\t", header=None, names=["ref", "sample", "mash_dist", "p", "matching_hashes"])])
        df.to_csv(f"{self.a['folder_stem']}evaluation/mash_results_full.csv")
        return df

    def create_report(self, df):
        '''Call report generator'''
        cls = GenerateReport(
            self.contig_graph_fname,
            self.consen_graph_fname,
            self.a,
            df
        )
        cls.main()

    def main(self) -> None:
        '''Retrieve all varieties of cons seq, graph'''
        self.get_consensus_seqs()
        self.call_alignment()
        self.call_graph(self.aln_fname, "consensus_vs_true_seq")

        '''Call graph on contig consensuses'''
        self.call_graph(
            f"experiments/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_consensus_alignment.aln", "contig_vs_ref_consensus_alignments")

        '''Do stats'''
        self.call_mash_dist()
        stats_df = self.collate_stats()
        self.create_report(stats_df)


if __name__ == "__main__":
    cls = Evaluate(
        {
            "ExpName": "ERR10812875",
            "SeqName": "ERR10812875",
            "GtOrg": "Paramyxoviridae_RSV",
            "GtFile": "data/rsv_set_metadata.csv",
            "StartTime": time.time()
        }
    )
    cls.main()
