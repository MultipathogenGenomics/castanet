import os
import numpy as np
import pandas as pd

from app.utils.timer import timing
from app.utils.shell_cmds import shell, make_dir
from app.utils.utility_fns import read_fa, save_fa, get_reference_org
from app.utils.fnames import get_consensus_fnames
from app.utils.system_messages import end_sec_print
from app.utils.basic_cli_calls import (samtools_index, bwa_index, bam_to_fastq,
                                       find_and_delete, rm)
from app.utils.error_handlers import error_handler_consensus_ref_corrected


class Consensus:
    '''Take all targets in one probetype/species aggregation, call consensus for each,
    flatten consensuses into single sequence.'''

    def __init__(self, payload) -> None:
        self.a = payload
        self.a["folder_stem"] = f"experiments/{self.a['ExpName']}/"
        self.target_consensuses = {}
        self.refs = read_fa(self.a["RefStem"])
        self.probe_names = pd.read_csv(
            f"experiments/{self.a['ExpName']}/probe_aggregation.csv")
        self.fnames = get_consensus_fnames(self.a)
        make_dir(f"mkdir {self.a['folder_stem']}consensus_data/")

    def filter_bam(self, tar_name) -> None:
        '''Filter bam to specific target, call consensus sequence for sam alignment records, grouped by target'''
        print(f"INFO: "
              f"Calling consensuses on all targets for: {tar_name}")
        shell(f"samtools view -b {self.fnames['master_bam']} {tar_name} "
              f"> {self.a['folder_stem']}grouped_reads/{tar_name}/{tar_name}.bam")
        shell(
            f"samtools consensus --call-fract 0.9 --min-depth {self.a['ConsensusMinD']} -f fasta {self.a['folder_stem']}grouped_reads/{tar_name}/{tar_name}.bam -o {self.a['folder_stem']}grouped_reads/{tar_name}/consensus_seqs_{tar_name}.fasta")

    def collate_consensus_seqs(self, tar_name) -> None:
        '''Read and collate consensus seqs from per target to per organism'''
        seqs_and_refs = [i for i in read_fa(
            f"{self.a['folder_stem']}grouped_reads/{tar_name}/consensus_seqs_{tar_name}.fasta") if tar_name in i[0]]
        seqs_and_refs = [[self.aggregate_to_probename(
            i[0]), i[0], i[1]] for i in seqs_and_refs]    # aggn, refn, seq

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
            np.isin(self.probe_names["target_id"], ref.replace(">", "")))[0]]
        if match.empty:
            return f"Unmatched"
        else:
            return f"{match['probetype'].item()}"

    def call_flat_consensus(self, org_name) -> None:
        '''Create consensus sequences'''
        if not os.path.isdir(f"{self.a['folder_stem']}consensus_data/{org_name}/"):
            shell(f"mkdir {self.a['folder_stem']}consensus_data/{org_name}/")

        if len(self.target_consensuses[org_name]) == 1:
            '''No flattening to be done, just a single target for this organism'''
            print(f"INFO: "
                  f"Only 1 target found for organism: {org_name} (current strategy is to not flatten)")
            flat_consensus = "".join(
                self.target_consensuses[org_name]["consensus_seq"])

        else:
            '''Otherwise, flatten consensus'''
            # RM < TODO pad consensus seqs to ~same length?
            self.build_msa_requisites(org_name)
            flat_consensus = self.flatten_consensus(org_name)

        '''Save'''
        save_fa(f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta",
                f">{org_name}_consensus\n{flat_consensus}")

        '''Filter bam to organism-specific targets, further filter by coverage %'''
        coverage_filter = self.filter_bam_to_organism(org_name)
        self.filter_tar_consensuses(org_name, coverage_filter)

        if len(self.target_consensuses[org_name]) > 0:
            '''Re-make target alignment and consensus to filtered list, save'''
            self.build_msa_requisites(org_name)
            flat_consensus = self.flatten_consensus(org_name)
            save_fa(f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta",
                    f">{org_name}_consensus\n{flat_consensus}")

            '''Remap to re-made flat consensus, to make `re-mapped consensus`'''
            self.remap_flat_consensus(org_name)
        else:
            print(
                f"INFO: No remapped consensus will be generated for {org_name} as coverage was too low on all target consensues")

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

    def flatten_consensus(self, org_name) -> str:
        '''Make MSA of references, then add fragments from target consensuses'''
        print(f"INFO: "
              f"making consensus alignments for target group: {org_name}")
        shell(f"mafft --thread {os.cpu_count()} {self.fnames['flat_cons_refs']} > {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_ref_alignment.aln",
              "Mafft align ref seqs (CONSENSUS.PY)")
        shell(f"mafft --thread {os.cpu_count()} --6merpair --addfragments {self.fnames['flat_cons_seqs']} {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_ref_alignment.aln "
              f"> {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln",
              "Mafft align consensus with ref seqs (CONSENSUS.PY)")

        '''Make flat consensus'''
        return self.dumb_consensus(np.array([list(i[1]) for i in read_fa(
            f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln")]), False)

    def dumb_consensus(self, aln, gap) -> str:
        '''Constrcut flat consensus to no reference'''
        cons, len_max = "", aln.shape[1]
        for i in range(len_max):
            hits, cnt = np.unique(aln[:, i], return_counts=True)
            tophit = hits[np.argsort(cnt)[-1]]
            if gap:
                '''If allowing gaps, take best hit'''
                cons += tophit
            else:
                '''If filling gaps (i.e. from ref)...'''
                if hits.shape[0] == 0:
                    '''... if gap only hit, put gap...'''
                    cons += tophit
                else:
                    '''...else take best non gap hit'''
                    cons += tophit if tophit != "-" else hits[np.argsort(
                        cnt)[-2]]
        return cons

    def filter_bam_to_organism(self, org_name) -> list:
        '''Retrieve directories for all targets in org grouping for consequent BAM merge'''
        target_dirs = [i for i in os.listdir(
            self.fnames['grouped_reads_dir']) if i in
            [i["tar_name"].replace(">", "") for i in self.target_consensuses[org_name] if i["tar_name"].startswith(">")]]

        '''Merge all bam files in grouped reads dir where they correspond to current target group'''
        shell(f"""samtools merge {self.a['folder_stem']}consensus_data/{org_name}/collated_reads_unf.bam {' '.join([f'{self.a["folder_stem"]}grouped_reads/{i}/{i}.bam' for i in target_dirs])}""",
              "Samtools merge, ref-adjusted consensus call (CONSENSUS.PY)")
        '''Output coverage stats for target consensuses'''
        shell(f"samtools coverage {self.a['folder_stem']}consensus_data/{org_name}/collated_reads_unf.bam | "
              f"grep {org_name} > {self.a['folder_stem']}consensus_data/{org_name}/target_consensus_coverage.csv")

        '''Get coverage for each consensus, further filter collated bam to targets > n coverage'''
        coverage_df = pd.read_csv(f"{self.a['folder_stem']}consensus_data/{org_name}/target_consensus_coverage.csv", sep="\t",
                                  names=["tar_name", "start_pos", "end_pos", "n_reads", "cov_bs", "cov", "mean_d",  "mean_b_q", "mean_m_q"])
        coverage_df = coverage_df[coverage_df["cov"]
                                  >= self.a['ConsensusCoverage']]
        coverage_df.to_csv(
            f"{self.a['folder_stem']}consensus_data/{org_name}/target_consensus_coverage.csv", index=False, header=False)
        coverage_filter = coverage_df["tar_name"].tolist()

        '''Index unfiltered bam, then filter for targets with sufficient coverage'''
        samtools_index(
            f"{self.a['folder_stem']}consensus_data/{org_name}/collated_reads_unf.bam")
        shell(f"samtools view -b {self.a['folder_stem']}consensus_data/{org_name}/collated_reads_unf.bam {' '.join([f'{i}' for i in coverage_filter])} "
              f"> {self.a['folder_stem']}consensus_data/{org_name}/collated_reads.bam")

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
              f"./bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta - | "
              f"viral_consensus -i - -r {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta -o {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_remapped_consensus_sequence.fasta --min_depth {self.a['ConsensusMinD']} --out_pos_counts {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_pos_counts.tsv")
        shell(f"samtools stats {self.a['folder_stem']}consensus_data/{org_name}/collated_reads.bam "
              f"> {self.a['folder_stem']}consensus_data/{org_name}/collated_reads_stats.stats")
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
        shell(f"mkdir {self.fnames['temp_folder']}")

        '''Retrieve GT seq'''
        ref_seq = get_reference_org(
            self.a['GtFile'], self.a["SeqName"], self.a['folder_stem'])

        '''Save ref seq and index, then run bwa mem of ref seq vs contigs'''
        print(f"INFO: "
              f"generating reference-adjusted consensus for target group / reference: {tar_name} / {ref_seq[0]}")
        save_fa(f"{self.fnames['temp_ref_seq']}",
                f">{ref_seq[0]}\n{ref_seq[1]}\n")
        bwa_index(f"{self.fnames['temp_ref_seq']}")

        '''Run alignment and flatten consensus'''
        bam_to_fastq(f"{self.a['folder_stem']}consensus_data/{tar_name}/collated_reads.bam",
                     f"{self.fnames['collated_reads_fastq']}")
        shell(f"./bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem {self.fnames['temp_ref_seq']} {self.fnames['collated_reads_fastq']} | "
              f"viral_consensus -i - -r {self.fnames['temp_ref_seq']} -o {ref_adj_cons_fname} --min_depth {self.a['ConsensusMinD']}")

    def tidy(self) -> None:
        '''Remove intermediate files to save disc space'''
        rm(f"{self.fnames['collated_reads_fastq']}")
        rm(f"{self.fnames['temp_folder']}", "-r")
        rm(f"{self.fnames['flat_cons_seqs']} {self.fnames['flat_cons_refs']}")
        find_and_delete(
            f"{self.a['folder_stem']}grouped_reads/ {self.a['folder_stem']}consensus_data/", "*.bam")

    def main(self) -> None:
        '''Entrypoint. Index main bam, filter it, make target consensuses, then create flattened consensus'''
        end_sec_print(
            "INFO: Calling consensus sequences\nThis may take a little while...")
        samtools_index(f"{self.a['folder_stem']}{self.a['SeqName']}.bam")
        for tar_name in os.listdir(f"{self.a['folder_stem']}grouped_reads/"):
            self.filter_bam(tar_name)

        '''Consensus for each thing target group'''
        [self.collate_consensus_seqs(tar_name) for tar_name in os.listdir(
            f"{self.a['folder_stem']}/grouped_reads/")]
        [self.call_flat_consensus(i) for i in self.target_consensuses.keys()]
        [self.call_ref_corrected_consensus(tar_name)
            for tar_name in self.target_consensuses.keys()]

        '''Tidy up'''
        self.tidy()
        end_sec_print("INFO: Consensus calling complete")


if __name__ == "__main__":
    cls = Consensus()
    cls.main()
