import os
import time
import pandas as pd
import pickle as p
import rapidfuzz as rf

from app.utils.shell_cmds import make_dir, shell, stoperr, loginfo, logerr
from app.utils.system_messages import end_sec_print
from app.utils.utility_fns import read_fa, save_fa, get_reference_org
from app.utils.similarity_graph import call_graph
from app.utils.report_gen import GenerateReport


class Evaluate:
    '''RM < DEPRECATED< USED FOR ORIGINAL EVALUATION'''

    def __init__(self, payload) -> None:
        self.a = payload
        if not "StartTime" in self.a.keys():
            self.a["StartTime"] = time.time()
        self.a["folder_stem"] = f"{self.a['SaveDir']}/{self.a['ExpName']}/"
        self.aln_fname = f"{self.a['folder_stem']}evaluation/{self.a['GtOrg']}_consensuses.aln"
        self.seq_names = [f">{self.a['GtOrg']}_GenBank_seq", f">{self.a['GtOrg']}_flattened_consensus",
                          f">{self.a['GtOrg']}_remapped_consensus", f">{self.a['GtOrg']}_gold_standard_consensus"]
        self.contig_graph_fname = f"{self.a['folder_stem']}evaluation/contig_vs_ref_consensus_alignments.png"
        self.consen_graph_fname = f"{self.a['folder_stem']}evaluation/consensus_vs_true_seq.png"
        self.additional_stats = {}
        make_dir(
            f"mkdir {self.a['folder_stem']}evaluation/ {self.a['folder_stem']}evaluation/seqs/")

    def get_consensus_seqs(self) -> None:
        '''Read in reference and consensus seqs, rename seqs, save versions to eval folder for stats and graphing'''
        ref_seq_present = True
        all_seqs = [get_reference_org(
            self.a['GtFile'], self.a['SeqName'], self.a['folder_stem'])]
        if all_seqs[0][0] == ">NO REFERENCE AVAILABLE":
            '''If no GT viral sequence, remove empty entry from aln'''
            ref_seq_present = False
        all_seqs.append(  # [1] Un-remapped flattened
            *read_fa(f"{self.a['SaveDir']}/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_consensus_sequence.fasta"))
        all_seqs.append(  # [2] Remapped
            *read_fa(f"{self.a['SaveDir']}/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_remapped_consensus_sequence.fasta"))
        all_seqs.append(  # [3] Ground truth
            *read_fa(f"{self.a['SaveDir']}/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_ref_adjusted_consensus.fasta"))

        if not ref_seq_present:
            '''Remove reference and ref-adjusted seqs if GT not present'''
            del (all_seqs[-1])
            del (all_seqs[0])

        '''Generate additional eval stats from raw seqs'''
        self.additional_stats["gt_len"] = len(all_seqs[0][1])
        self.additional_stats["remap_len"] = len(all_seqs[-2][1])
        self.additional_stats["lev_d"] = rf.distance.Levenshtein.distance(
            all_seqs[0][1], all_seqs[-2][1])
        '''Fuzzy ratio = normalized indel distance'''
        self.additional_stats["fuzz_r"] = round(
            rf.fuzz.ratio(all_seqs[0][1].lower(), all_seqs[-2][1].lower()), 1)
        '''JaroWinkler distance = edit distance for transposition only'''
        self.additional_stats["jw_d"] = round(
            rf.distance.JaroWinkler.normalized_similarity(all_seqs[0][1].lower(), all_seqs[-2][1].lower()), 2)

        '''Rename seqs, save individual fa, build mash sketch and combined fa (for aln)'''
        for i in range(len(all_seqs)):
            all_seqs[i][0] = self.seq_names[i]
            save_fa(f"{self.a['folder_stem']}evaluation/seqs/{self.seq_names[i].replace('>','')}.fa",
                    f"{all_seqs[i][0]}\n{all_seqs[i][1]}\n")
            shell(
                f"mash sketch {self.a['folder_stem']}evaluation/seqs/{self.seq_names[i].replace('>','')}.fa")

        with open(f"{self.a['folder_stem']}evaluation/consensus_seqs.fasta", "w") as f:
            [f.write(f"{i[0]}\n{i[1]}\n") for i in all_seqs]

        return ref_seq_present

    def call_alignment(self) -> None:
        '''Build MAFFT MSA of all consensuses and "true" sequence'''
        loginfo(f"Calling Mafft alignment of consensus seqs")
        shell(f"mafft --auto --thread {os.cpu_count()} {self.a['folder_stem']}evaluation/consensus_seqs.fasta > {self.aln_fname}",
              "Mafft align consensus seqs (EVAL.PY)")
        if len(read_fa(self.aln_fname)) == 0:
            stoperr(f"Alignmnet is empty - check Mafft didn't crash.")

    def call_mash_dist(self) -> None:
        '''Call mash distances on (1) all consensuses vs original viral genome, (2) flat/remapped cons vs gold standard cons'''
        shell(f"mash dist {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_GenBank_seq.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_flattened_consensus.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_remapped_consensus.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_gold_standard_consensus.fa.msh > {self.a['folder_stem']}evaluation/mash_results_cons_vs_true_genome.csv")
        shell(f"mash dist {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_gold_standard_consensus.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_flattened_consensus.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_remapped_consensus.fa.msh > {self.a['folder_stem']}evaluation/mash_results_cons_vs_gold_standard.csv")

    def do_unreferenced_eval(self) -> None:
        '''Mash distances in the absence of references'''
        shell(f"mash dist {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_flattened_consensus.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_remapped_consensus.fa.msh > {self.a['folder_stem']}evaluation/mash_results_cons_vs_true_genome.csv")
        shell(f"mash dist {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_flattened_consensus.fa.msh {self.a['folder_stem']}evaluation/seqs/{self.a['GtOrg']}_remapped_consensus.fa.msh > {self.a['folder_stem']}evaluation/mash_results_cons_vs_gold_standard.csv")

    def collate_stats(self) -> pd.DataFrame:
        '''Collect stats sets on: target consensus, remapped consensus, MASH values between consensus types'''
        '''Stats on target consensuses. b = base, m = map, cov = coverage, d = depth'''
        t_df = pd.read_csv(f"{self.a['folder_stem']}consensus_data/{self.a['GtOrg']}/target_consensus_coverage.csv",
                           names=["tar_name", "start_pos", "end_pos", "n_reads", "cov_bs", "cov", "mean_d",  "mean_b_q", "mean_m_q"])
        t_df = t_df[t_df["n_reads"] != 0]  # must have > 1 read
        t_df = t_df.round(2)
        t_stats = [t_df.columns.tolist()] + t_df.values.tolist()

        ''''Stats on remapped consensus'''
        c_df = pd.read_csv(
            f"{self.a['folder_stem']}consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_consensus_pos_counts.tsv")  # remapped cons stats
        gc = round((c_df["G"].sum() + c_df["C"].sum()) /
                   c_df["Total"].sum() * 100, 2)
        missing = c_df[c_df["Total"] == 0].shape[0]
        ambigs = c_df[(c_df["-"] != 0) & (c_df["A"] == 0) & (c_df["C"]
                                                             == 0) & (c_df["T"] == 0) & (c_df["G"] == 0)].shape[0]
        coverage = round(
            (1 - ((missing + ambigs) / c_df["Total"].sum())) * 100, 2)

        '''Get additional stats on consensus remapping'''
        with open(f"{self.a['folder_stem']}consensus_data/{self.a['GtOrg']}/supplementary_stats.p", 'rb') as f:
            self.additional_stats["n_remapped_seqs"] = p.load(
                f)["filtered_collated_read_num"]
        cons_vs_ref_len = round(
            (self.additional_stats['remap_len'] / self.additional_stats['gt_len']) * 100)
        c_stats = [["n_bases", "n_reads", "len_%", "gc_pc", "missing_bs", "ambig_bs", "cov", "lev_d", "fuzz_r", "jw_d"],
                   [c_df["Total"].sum(), self.additional_stats['n_remapped_seqs'], cons_vs_ref_len, gc, missing, ambigs, coverage,
                    self.additional_stats['lev_d'], self.additional_stats['fuzz_r'], self.additional_stats['jw_d']]]

        '''Grab MASH output, put in dataframe, save'''
        m_df = pd.read_csv(f"{self.a['folder_stem']}evaluation/mash_results_cons_vs_true_genome.csv",
                           delimiter="\t", header=None, names=["ref", "sample", "mash_dist", "p", "matching_hashes"])
        m_df = pd.concat([m_df, pd.read_csv(f"{self.a['folder_stem']}evaluation/mash_results_cons_vs_gold_standard.csv",
                                            delimiter="\t", header=None, names=["ref", "sample", "mash_dist", "p", "matching_hashes"])])
        m_df["p"] = m_df["p"].round(3)
        m_df.to_csv(f"{self.a['folder_stem']}evaluation/mash_results_full.csv")

        return t_stats, c_stats, m_df

    def create_report(self, t_stats, c_stats, stats_df):
        '''Call report generator'''
        cls = GenerateReport(
            self.contig_graph_fname,
            self.consen_graph_fname,
            self.a,
            t_stats, c_stats, stats_df
        )
        cls.main()

    def main(self) -> None:
        '''Retrieve all varieties of cons seq, graph'''
        try:
            ref_seq_present = self.get_consensus_seqs()
            self.call_alignment()

            call_graph(self.a['SeqName'], self.a['GtOrg'],
                       self.aln_fname, "consensus_vs_true_seq")

            '''Call graph on contig consensuses'''
            call_graph(self.a['SeqName'], self.a['GtOrg'],
                       f"{self.a['SaveDir']}/{self.a['ExpName']}/consensus_data/{self.a['GtOrg']}/{self.a['GtOrg']}_consensus_alignment.aln", "contig_vs_ref_consensus_alignments")

            '''Do stats'''
            if ref_seq_present:
                self.call_mash_dist()
            else:
                self.do_unreferenced_eval()

            t_stats, c_stats, stats_df = self.collate_stats()
            self.create_report(t_stats, c_stats, stats_df)

        except Exception as e:
            logerr(
                f"WARNING: I didn't do post hoc evaluation on this sample, but this might have been deliberate if you didn't specify a ground truth organism/metadata file. Exception: {e}")
        end_sec_print("Eval complete")


if __name__ == "__main__":
    cls = Evaluate(
        {
            "ExpName": "ERR10812876",
            "SeqName": "ERR10812876",
            "GtOrg": "Paramyxoviridae_RSV",
            "GtFile": "data/rsv_set_metadata.csv",
            "StartTime": time.time()
        }
    )
    cls.main()
