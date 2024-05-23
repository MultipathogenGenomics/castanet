import os
import numpy as np
import pandas as pd
import pickle as p
from Bio import AlignIO
from collections import Counter
import plotly.express as px

from app.utils.timer import timing
from app.utils.shell_cmds import shell, make_dir, loginfo, stoperr
from app.utils.utility_fns import read_fa, save_fa, get_reference_org
from app.utils.fnames import get_consensus_fnames
from app.utils.system_messages import end_sec_print
from app.utils.basic_cli_calls import (
    samtools_index, bwa_index, find_and_delete, rm, samtools_read_num)
from app.utils.error_handlers import error_handler_consensus_ref_corrected
from app.utils.similarity_graph import call_graph


class Consensus:
    '''Take all targets in one probetype/species aggregation, call consensus for each,
    flatten consensuses into single sequence.'''

    def __init__(self, payload) -> None:
        self.a = payload
        self.a["folder_stem"] = f"{self.a['ExpRoot']}/{self.a['ExpName']}/"
        self.target_consensuses = {}
        self.insufficient_coverage_orgs = []
        self.refs = read_fa(self.a["RefStem"])
        self.probe_names = pd.read_csv(
            f"{self.a['ExpRoot']}/{self.a['ExpName']}/probe_aggregation.csv")
        self.fnames = get_consensus_fnames(self.a)
        self.eval_stats = {}
        make_dir(f"mkdir {self.a['folder_stem']}consensus_data/")

    def filter_bam(self, tar_name) -> None:
        '''Filter bam to specific target, call consensus sequence for sam alignment records, grouped by target'''
        loginfo(f"Calling consensuses on all targets for: {tar_name}")
        group_bam_fname = f"{self.a['folder_stem']}grouped_reads/{tar_name}/{tar_name}.bam"
        shell(f"samtools view -b {self.fnames['master_bam']} {tar_name} "
              f"> {group_bam_fname}")
        match_str = f"(^|\s){tar_name}"
        if len(tar_name) < 99:
            match_str = f"{match_str}($|\s)"
        out = shell(f"samtools coverage {group_bam_fname} | grep -E '{match_str}'",
                    "Coverage, consensus filter bam", ret_output=True)

        try:
            out = out.decode().replace("\n", "").split("\t")[6:9]
        except Exception as ex:
            raise loginfo(
                f"Could not generate a consensus for target: {tar_name}\nException: {ex}")

        assert len(
            out) > 0, f"Could not generate a consensus for target: {tar_name}\nNo coverage detected for this target, this usually happens because your probe naming scheme is incompatible with castanet"

        if float(out[0]) < self.a['ConsensusMinD'] or float(out[2]) < self.a["ConsensusMapQ"]:
            '''If coverage/depth don't surpass threshold, delete grouped reads dir'''
            loginfo(
                f"Not adding subconsensus for {tar_name} to consensus for organism, as min D was {out[0]} and Map Q was {out[2]}")
            # TODO DISABLED FOR TEST
            shell(f"rm -r {self.a['folder_stem']}grouped_reads/{tar_name}/")
            return
        else:
            '''Else, call consensus on this target'''
            shell(
                f"samtools consensus --call-fract 0.9 --min-depth {self.a['ConsensusMinD']} -f fasta {group_bam_fname} -o {self.a['folder_stem']}grouped_reads/{tar_name}/consensus_seqs_{tar_name}.fasta")

    def collate_consensus_seqs(self, tar_name) -> None:
        '''Read and collate consensus seqs from per target to per organism'''
        try:
            seqs_and_refs = [i for i in read_fa(
                f"{self.a['folder_stem']}grouped_reads/{tar_name}/consensus_seqs_{tar_name}.fasta") if tar_name in i[0]]
            seqs_and_refs = [[self.aggregate_to_probename(
                i[0]), i[0], i[1]] for i in seqs_and_refs]    # aggn, refn, seq
        except Exception as ex:
            print(f"***{ex}")
            return

        if len(seqs_and_refs) == 0:
            return

        consensus_org = seqs_and_refs[0][0]
        if not consensus_org in self.target_consensuses.keys():
            self.target_consensuses[consensus_org] = []

        for i in seqs_and_refs:
            self.target_consensuses[consensus_org].append({
                "tar_name": i[1],
                "consensus_seq":  ', '.join([i[2] for i in seqs_and_refs])
            })

    def aggregate_to_probename(self, ref) -> str:
        '''Group targets to organism via probe name (compiled in analysis.py)'''
        match = self.probe_names.iloc[np.where(
            np.isin(self.probe_names["orig_target_id"], ref.replace(">", "")))[0]]
        if match.empty:
            print(
                f"WARNING: Couldn't match reads to probe name: {self.probe_names['target_id']}")
            return f"Unmatched"
        else:
            return f"{match['probetype'].item()}"

    def call_flat_consensus(self, org_name) -> None:
        '''Create consensus sequences'''
        '''Make folder and dictionary key for supplementary stats'''
        self.eval_stats[org_name] = {}
        if not os.path.isdir(f"{self.a['folder_stem']}consensus_data/{org_name}/"):
            shell(f"mkdir {self.a['folder_stem']}consensus_data/{org_name}/")

        '''Filter bam to organism-specific targets, further filter by coverage %'''
        coverage_filter = self.filter_bam_to_organism(org_name)

        if len(coverage_filter) == 0:
            loginfo(
                f"No remapped consensus will be generated for {org_name} as coverage was too low on all target consensues")
            self.insufficient_coverage_orgs.append(org_name)
            return

        '''Filter tar consensuses on coverage, re-make target alignment and consensus to filtered list, save'''
        self.filter_tar_consensuses(org_name, coverage_filter)
        self.build_msa_requisites(org_name)
        flat_consensus = self.flatten_consensus(org_name)
        save_fa(f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta",
                f">{org_name}_consensus\n{flat_consensus}")

        '''Remap to re-made flat consensus, to make `re-mapped consensus`'''
        self.remap_flat_consensus(org_name)

        '''Dump any additional stats to pickle'''
        self.dump_stats(org_name)

    def dump_stats(self, org_name) -> None:
        with open(f"{self.a['folder_stem']}consensus_data/{org_name}/supplementary_stats.p", 'wb') as f:
            p.dump(self.eval_stats[org_name], f, protocol=p.HIGHEST_PROTOCOL)

    def build_msa_requisites(self, org_name) -> None:
        '''Create fasta files containing target reference seqs and consensus seqs, for downstream MSA'''
        ref_seq_names = list(set([i["tar_name"].replace(">", "")
                                  for i in self.target_consensuses[org_name]]))
        ref_seqs = [ref for ref in self.refs if ref[0].replace(
            ">", "") in ref_seq_names]
        with open(self.fnames['flat_cons_refs'], "w") as f:
            [f.write(f"{i[0]}\n{i[1]}\n") for i in ref_seqs]
        with open(self.fnames['flat_cons_seqs'], "w") as f:
            [f.write(f"{i['tar_name']}_CONS\n{i['consensus_seq']}\n")
             for i in self.target_consensuses[org_name]]
        shell(
            f"cat {self.fnames['flat_cons_seqs']} {self.fnames['flat_cons_refs']} > {self.a['folder_stem']}/consensus_data/unaligned_consensuses_and_refs.fna")

    def flatten_consensus(self, org_name) -> str:
        '''Make MSA of references, then add fragments from target consensuses'''
        loginfo(f"making consensus alignments for target group: {org_name}")
        ref_aln_fname = f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_ref_alignment.aln"
        shell(f"mafft --thread {os.cpu_count()} --auto {self.fnames['flat_cons_refs']} > {ref_aln_fname}",
              "Mafft align ref seqs (CONSENSUS.PY)")
        refs = read_fa(ref_aln_fname)
        if len(refs) == 0:
            '''If only 1 reference seq, the alignment wouldn't have worked - defer to temp refs file in these cases'''
            ref_aln_fname = self.fnames['flat_cons_refs']
        shell(f"mafft --thread {os.cpu_count()} --auto --addfragments {self.fnames['flat_cons_seqs']} {ref_aln_fname}"
              f"> {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln",
              "Mafft align consensus with ref seqs (CONSENSUS.PY)")
        try:
            call_graph(self.a["ExpName"], org_name, f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln",
                       f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_target_consensus_alignment", self.a["ExpRoot"], is_eval=False)
        except FileNotFoundError:
            raise SystemError(
                "Castanet couldn't construct a consensus alignment graph")

        '''Return flat consensus'''
        return self.dumb_consensus(f"{self.a['folder_stem']}consensus_data/{org_name}/", org_name)

    @timing
    def dumb_consensus_deprecated(self, aln, org_name) -> str:
        '''DEPRECATED. Constrcut flat consensus to no reference'''
        aln = np.array([list(i[1]) for i in read_fa(
            f"{aln}{org_name}_consensus_alignment.aln")])
        cons, len_max = "", aln.shape[1]
        for i in range(len_max):
            hits, cnt = np.unique(aln[:, i], return_counts=True)
            cons += hits[np.argsort(cnt)[-1]]

        return cons

    @timing
    def dumb_consensus(self, alnfpath, org_name) -> list:
        '''Produce an un-referenced/`flat` consensus sequence for file of target and target ref seqs'''
        def base_cons(s):
            ''' Return strict consensus for a set of bases (eg. column in alignment), ignoring gaps. '''
            len_max = len(s)
            s = s.replace('-', '')
            if not s or (len(s) <= 0.1 * len_max):
                return ('', np.nan)
            consbase, consnum = Counter(s.lower()).most_common()[0]
            return consbase, float(consnum)/len(s)

        aln = AlignIO.read(
            f"{alnfpath}{org_name}_consensus_alignment.aln", 'fasta')
        cluster_cons = pd.DataFrame(
            pd.Series(base_cons(aln[:, i])) for i in range(len(aln[0]))).dropna()

        '''Plot identity for QC'''
        cluster_cons.columns = ['cons', 'ident']
        cluster_cons["ident"].rolling(120).mean().plot()
        fig = px.line(x=cluster_cons.index,
                      y=cluster_cons["ident"], title="Flat consensus identity")
        fig.write_image(f"{alnfpath}{org_name}_flat_consensus_identity.png")

        return "".join(cluster_cons["cons"].tolist())

    def filter_bam_to_organism(self, org_name) -> list:
        '''Retrieve directories for all targets in org grouping for consequent BAM merge'''
        target_dirs = [i for i in os.listdir(
            self.fnames['grouped_reads_dir']) if i in
            [i["tar_name"].replace(">", "") for i in self.target_consensuses[org_name] if i["tar_name"].startswith(">")]]

        '''Merge all bam files in grouped reads dir where they correspond to current target group'''
        shell(f"""samtools merge {self.a['folder_stem']}consensus_data/{org_name}/collated_reads_unf.bam {' '.join([f'{self.a["folder_stem"]}grouped_reads/{i}/{i}.bam' for i in target_dirs])}""",
              "Samtools merge, ref-adjusted consensus call (CONSENSUS.PY)")
        '''Output coverage stats for target consensuses'''
        try:
            shell(f"samtools coverage {self.a['folder_stem']}consensus_data/{org_name}/collated_reads_unf.bam | "
                  f"grep -E '{'|'.join(self.probe_names[self.probe_names['probetype'] == org_name]['orig_target_id'].tolist())}' > {self.a['folder_stem']}consensus_data/{org_name}/target_consensus_coverage.csv")
        except OSError as e:
            stoperr(f"We couldn't call a consensus because your list of probes is too long. This is a kernel limitation. Can you reduce the size of your probe panel?")

        '''Get coverage for each consensus, filter collated bam by consensus coverage and map q'''
        coverage_df = pd.read_csv(f"{self.a['folder_stem']}consensus_data/{org_name}/target_consensus_coverage.csv", sep="\t",
                                  names=["tar_name", "start_pos", "end_pos", "n_reads", "cov_bs", "cov", "mean_d",  "mean_b_q", "mean_m_q"])
        coverage_df = coverage_df[(coverage_df["cov"] >= self.a['ConsensusCoverage']) & (
            coverage_df["mean_m_q"] >= self.a['ConsensusMapQ'])]
        coverage_df.to_csv(
            f"{self.a['folder_stem']}consensus_data/{org_name}/target_consensus_coverage.csv", index=False, header=False)
        coverage_filter = coverage_df["tar_name"].tolist()

        if len(coverage_filter) == 0:
            return []
        else:
            '''Index unfiltered bam, then filter for targets with sufficient coverage'''
            samtools_index(
                f"{self.a['folder_stem']}consensus_data/{org_name}/collated_reads_unf.bam")
            shell(f"samtools view -b {self.a['folder_stem']}consensus_data/{org_name}/collated_reads_unf.bam {' '.join([f'{i}' for i in coverage_filter])} "
                  f"> {self.a['folder_stem']}consensus_data/{org_name}/collated_reads.bam")

            '''Estimate number of mapped reads in the final alignment (get just primary mapped reads, div 2 to average F & R strands)'''
            self.eval_stats[org_name]["filtered_collated_read_num"] = round(samtools_read_num(
                f"{self.a['folder_stem']}consensus_data/{org_name}", "collated_reads", '-F 0x904 -q 20') / 2)
            return coverage_filter

    def filter_tar_consensuses(self, org_name, filter) -> None:
        '''Purge target consensus from master list if coverage was lower than threshold (aln is consequently remade)'''
        to_del = [i for i in range(len(self.target_consensuses[org_name]))
                  if not self.target_consensuses[org_name][i]["tar_name"].replace(">", "") in filter]

        to_del.reverse()
        for i in to_del:
            '''Not done in loop and in reverse to not break the iterator'''
            del self.target_consensuses[org_name][i]

    @timing
    def remap_flat_consensus(self, org_name) -> None:
        '''Remap reads to flattened consensus, save, call stats, remove raw fastas'''
        bwa_index(
            f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta")
        shell(f"samtools fastq {self.a['folder_stem']}consensus_data/{org_name}/collated_reads.bam |"
              f"bwa-mem2 mem -t {self.a['NThreads']} {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta - | "
              f"./ViralConsensus/viral_consensus -i - -r {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta -o {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_remapped_consensus_sequence.fasta --min_depth {self.a['ConsensusMinD']} --out_pos_counts {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_pos_counts.tsv")
        try:
            self.fix_terminal_gaps(f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_pos_counts.tsv",
                                f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_remapped_consensus_sequence.fasta")
        except FileNotFoundError as ex:
            stoperr(f"Castanet call to ViralConsensus produced empty output. Check that it's installed using the dependency_check endpoint (see readme)")
        find_and_delete(
            f"{self.a['folder_stem']}consensus_data/{org_name}/", f"*.fasta.*")

    @timing
    def call_ref_corrected_consensus(self, tar_name) -> None:
        '''Construct a `conventional` consensus from grouped reads, with reference to a complete reference genome'''
        '''Load grouped aligned first consensus seqs to retrieve each target name'''
        if error_handler_consensus_ref_corrected(self.a, tar_name):
            return

        '''Set dynamic fnames, make folders'''
        ref_adj_cons_fname = f"{self.a['folder_stem']}consensus_data/{tar_name}/{tar_name}_ref_adjusted_consensus.fasta"
        outcounts_fname = f"{self.a['folder_stem']}consensus_data/{tar_name}/{tar_name}_refadjconsensus_pos_counts.tsv"
        shell(f"mkdir {self.fnames['temp_folder']}")

        '''Retrieve GT seq'''
        ref_seq = get_reference_org(
            self.a['GtFile'], self.a["ExpName"], self.a['folder_stem'])

        '''Save ref seq and index, then run bwa mem of ref seq vs contigs'''
        loginfo(
            f"generating reference-adjusted consensus for target group / reference: {tar_name} / {ref_seq[0]}")
        save_fa(f"{self.fnames['temp_ref_seq']}",
                f">{ref_seq[0]}\n{ref_seq[1]}\n")
        bwa_index(f"{self.fnames['temp_ref_seq']}")

        '''Run alignment and flatten consensus'''
        shell(f"samtools fastq {self.a['folder_stem']}consensus_data/{tar_name}/collated_reads.bam | "
              f"bwa-mem2 mem {self.fnames['temp_ref_seq']} - -t {self.a['NThreads']} | ./ViralConsensus/viral_consensus -i - -r {self.fnames['temp_ref_seq']} -o {ref_adj_cons_fname} --out_pos_counts {outcounts_fname}")
        self.fix_terminal_gaps(outcounts_fname, ref_adj_cons_fname)

    def fix_terminal_gaps(self, in_fname, out_fname) -> None:
        '''Trim terminal gaps'''
        cons = pd.read_csv(in_fname, sep="\t")
        n_pos = cons.shape[0]
        '''If total reads at pos x < threshold AND in leading/trailing 5% of reads, mark for deletion'''
        cons["del"] = cons.apply(lambda x: np.where(x["Total"] < 30 and (
            x["Pos"] < n_pos * 0.05 or x["Pos"] > n_pos * 0.95), 1, 0), axis=1)
        '''Rm terminal gaps, re-index, re-call consensus.'''
        cons = cons[cons["del"] == 0]
        cons = cons.drop(columns=["Pos", "Total", "del"])
        cons["con"] = cons.apply(lambda x: x.idxmax(), axis=1)
        '''Re-do index and totals'''
        cons["Pos"] = np.arange(1, cons.shape[0] + 1)
        # RM < TODO PLOT CONS COVERAGE
        cons["Total"] = cons.apply(
            lambda x: x["A"] + x["T"] + x["C"] + x["G"] + x["-"], axis=1)
        '''Fix artificial adnylation where totals = 0 (pd idxmax annoyingly picks first col in this case)'''
        cons["con"] = cons.apply(
            lambda x: x["con"] if not x["Total"] == 0 else "-", axis=1)
        cons["-"] = cons.apply(lambda x: 1 if x["Total"] == 0 else 0, axis=1)
        cons["con"] = cons["con"].astype(str)
        cons.to_csv(in_fname)
        save_fa(out_fname, f">CONSENSUS\n{''.join(cons['con'].tolist())}")

    def clean_incomplete_consensus(self) -> None:
        '''If we had insufficient coverage for organism x, clean it from self vars and folder tree'''
        for org_name in self.insufficient_coverage_orgs:
            del self.target_consensuses[org_name]
            shell(f"rm -r {self.a['folder_stem']}/consensus_data/{org_name}/")

    def tidy(self) -> None:
        '''Remove intermediate files to save disc space'''
        rm(f"{self.fnames['collated_reads_fastq']}")
        rm(f"{self.fnames['temp_folder']}", "-r")
        rm(f"{self.fnames['flat_cons_seqs']} {self.fnames['flat_cons_refs']}")
        find_and_delete(
            f"{self.a['folder_stem']}grouped_reads/", "*.bam")
        if self.a['ConsensusCleanFiles']:
            find_and_delete(
                f"{self.a['folder_stem']}consensus_data/", "*.bam")

    def generate_summary(self, org) -> None:
        dfpath = f"{self.a['folder_stem']}/consensus_seq_stats.csv"
        cols = ["target", "n_bases", "n_reads",
                "gc_pc", "missing_bs", "ambig_bs", "cov"]

        if not os.path.isfile(dfpath):
            df = pd.DataFrame(columns=cols)
        else:
            df = pd.read_csv(dfpath)

        # remapped cons stats
        c_df = pd.read_csv(
            f"{self.a['folder_stem']}/consensus_data/{org}/{org}_consensus_pos_counts.tsv")
        gc = round((c_df["G"].sum() + c_df["C"].sum()) /
                   c_df["Total"].sum() * 100, 2)
        missing = c_df[c_df["Total"] == 0].shape[0]
        ambigs = c_df[(c_df["-"] != 0) & (c_df["A"] == 0) & (c_df["C"]
                                                             == 0) & (c_df["T"] == 0) & (c_df["G"] == 0)].shape[0]
        coverage = round(
            (1 - ((missing + ambigs) / c_df["Total"].sum())) * 100, 2)

        '''Get additional stats on consensus remapping'''
        additional_stats = {}
        with open(f"{self.a['folder_stem']}/consensus_data/{org}/supplementary_stats.p", 'rb') as f:
            additional_stats["n_remapped_seqs"] = p.load(
                f)["filtered_collated_read_num"]

        c_stats = pd.DataFrame([[org, c_df["Total"].sum(
        ), additional_stats['n_remapped_seqs'], gc, missing, ambigs, coverage]], columns=cols)
        df = pd.concat([df, c_stats], axis=0, ignore_index=True)
        df.to_csv(dfpath)

        '''Plot consensus coverage'''
        px.line(c_df, x="Pos", y="Total", title=f"Consensus coverage, {org} ({self.a['ExpName']})",
                labels={"Pos": "Position", "Total": "Num Reads"}).write_image(
            f"{self.a['folder_stem']}/consensus_data/{org}/{org}_consensus_coverage.png")

    def main(self) -> None:
        '''Entrypoint. Index main bam, filter it, make target consensuses, then create flattened consensus'''
        end_sec_print(
            "INFO: Calling consensus sequences\nThis may take a little while...")
        samtools_index(f"{self.a['folder_stem']}{self.a['ExpName']}.bam")

        for tar_name in os.listdir(f"{self.a['folder_stem']}grouped_reads/"):
            self.filter_bam(tar_name)

        '''Consensus for each thing target group'''
        [self.collate_consensus_seqs(tar_name) for tar_name in os.listdir(
            f"{self.a['folder_stem']}/grouped_reads/")]
        [self.call_flat_consensus(i) for i in self.target_consensuses.keys()]
        self.clean_incomplete_consensus()
        [self.call_ref_corrected_consensus(tar_name)
            for tar_name in self.target_consensuses.keys()]

        '''Tidy up'''
        self.tidy()

        '''Call CSV summary generator'''
        [self.generate_summary(i) for i in os.listdir(
            f"{self.a['folder_stem']}/consensus_data/") if not "GROUND_TRUTH" in i and not ".fna" in i]

        end_sec_print("INFO: Consensus calling complete")


if __name__ == "__main__":
    cls = Consensus()
    cls.main()
