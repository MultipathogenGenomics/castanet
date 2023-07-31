import os
import numpy as np
import pandas as pd

from app.utils.timer import timing
from app.utils.shell_cmds import shell, make_dir
from app.utils.utility_fns import read_fa, save_fa, get_reference_org
from app.utils.fnames import get_consensus_fnames
from app.utils.system_messages import end_sec_print


class Consensus:
    '''Take all targets in one probetype/species aggregation, call consensus for each,
    flatten consensuses into single sequence.'''

    def __init__(self, payload) -> None:
        self.a = payload
        self.a["folder_stem"] = f"experiments/{self.a['ExpName']}/"
        self.consensus_seqs, self.consensus_refs = {}, {}
        self.refs = read_fa(self.a["RefStem"])
        self.probe_names = pd.read_csv(
            f"experiments/{self.a['ExpName']}/probe_aggregation.csv")
        self.fnames = get_consensus_fnames(self.a)
        make_dir(f"mkdir {self.a['folder_stem']}consensus_data/")
        if "GtFile" in payload.keys():
            self.do_ref_adjusted_cons = True
        else:
            self.do_ref_adjusted_cons = False

    def filter_bam(self, tar_name) -> None:
        '''Filter bam to specific target, call consensus sequence for sam alignment records, grouped by target'''
        print(f"INFO: "
              f"Calling consensuses on all targets for: {tar_name}")
        shell(
            f"samtools view -b {self.fnames['master_bam']} {tar_name} > {self.a['folder_stem']}grouped_reads/{tar_name}/{tar_name}.bam")
        shell(
            f"samtools consensus --min-depth {self.a['ConsensusMinD']} -f fasta {self.a['folder_stem']}grouped_reads/{tar_name}/{tar_name}.bam -o {self.a['folder_stem']}grouped_reads/{tar_name}/consensus_seqs_{tar_name}.fasta")

    def collate_consensus_seqs(self, tar_name):
        '''Read and collate to self var the consensus seqs from per target to per organism'''
        seqs_and_refs = [i for i in read_fa(
            f"{self.a['folder_stem']}grouped_reads/{tar_name}/consensus_seqs_{tar_name}.fasta") if tar_name in i[0]]
        # aggn, refn, seq
        seqs_and_refs = [[self.aggregate_to_probename(
            i[0]), i[0], i[1]] for i in seqs_and_refs]

        if len(seqs_and_refs) == 0:
            return

        consensus_org = seqs_and_refs[0][0]

        if not consensus_org in self.consensus_seqs.keys():
            self.consensus_seqs[consensus_org] = [i[2] for i in seqs_and_refs]
            self.consensus_refs[consensus_org] = [
                [[i[0], i[1]]] for i in seqs_and_refs]
        else:
            self.consensus_seqs[consensus_org].append(
                ', '.join([i[2] for i in seqs_and_refs]))
            self.consensus_refs[consensus_org].append(
                [[i[0], i[1]] for i in seqs_and_refs])

    def aggregate_to_probename(self, ref):
        match = self.probe_names.iloc[np.where(
            np.isin(self.probe_names["target_id"], ref.replace(">", "")))[0]]
        if match.empty:
            return f"Unmatched"
        else:
            return f"{match['probetype'].item()}"

    def call_flat_consensus(self, org_name):
        '''Create consensus sequences'''
        if not os.path.isdir(f"{self.a['folder_stem']}consensus_data/{org_name}/"):
            shell(f"mkdir {self.a['folder_stem']}consensus_data/{org_name}/")

        if len(self.consensus_seqs[org_name]) == 1:
            '''No flattening to be done, just a single target for this organism'''
            print(f"INFO: "
                  f"Only 1 target found for organism: {org_name} (current strategy is to not flatten)")
            flat_consensus = "".join(self.consensus_seqs[org_name])

        else:
            '''Otherwise, flatten consensus'''
            # RM < TODO pad consensus seqs to ~same length?
            ref_seq_names = list(set([i[0][1].replace(">", "")
                                      for i in self.consensus_refs[org_name]]))
            ref_seqs = [ref for ref in self.refs if ref[0].replace(
                ">", "") in ref_seq_names]
            with open(self.fnames['flat_cons_refs'], "w") as f:
                [f.write(f"{i[0]}\n{i[1]}\n") for i in ref_seqs]
            with open(self.fnames['flat_cons_seqs'], "w") as f:
                [f.write(f">TARGET_CONSENSUS_{i}\n{self.consensus_seqs[org_name][i]}\n") for i in range(
                    len(self.consensus_seqs[org_name]))]
            flat_consensus = self.flatten_consensus(org_name)

        '''Save'''
        save_fa(f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta",
                f">{org_name}_consensus\n{flat_consensus}")

        '''Remap, re-call consensus, save if multiple consensus merged'''
        self.filter_bam_to_organism(org_name)
        self.remap_flat_consensus(org_name)

    def flatten_consensus(self, org_name):
        '''Make MSA of references, then add fragments from target consensuses'''
        print(f"INFO: "
              f"making consensus alignments for target group: {org_name}")
        shell(f"mafft --thread {os.cpu_count()} {self.fnames['flat_cons_refs']} > {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_ref_alignment.aln",
              "Mafft align ref seqs (CONSENSUS.PY)")
        shell(f"mafft --thread {os.cpu_count()} --6merpair --addfragments {self.fnames['flat_cons_seqs']} {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_ref_alignment.aln > {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln",
              "Mafft align consensus with ref seqs (CONSENSUS.PY)")

        '''Make flat consensus'''
        flat_consensus = self.rich_consensus(np.array([list(i[1]) for i in read_fa(
            f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln")]), False)
        shell(
            f"rm {self.fnames['flat_cons_seqs']} {self.fnames['flat_cons_refs']}")
        return flat_consensus

    def rich_consensus(self, aln, gap):
        '''Constrcut `dumb` consensus'''
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

    def filter_bam_to_organism(self, org_name):
        '''Retrieve directories for all targets in org grouping for consequent BAM merge'''
        targets_for_group = [i[0][1].replace(
            ">", "") for i in self.consensus_refs[org_name] if i[0][1].startswith(">") and not "TARGET" in i]
        target_dirs = [i for i in os.listdir(
            self.fnames['grouped_reads_dir']) if i in targets_for_group]
        '''Merge all bam files in grouped reads dir where they correspond to current target group'''
        shell(f"""samtools merge {self.a['folder_stem']}consensus_data/{org_name}/collated_reads.bam {' '.join([f'{self.a["folder_stem"]}grouped_reads/{i}/{i}.bam' for i in target_dirs])}""",
              "Samtools merge, ref-adjusted consensus call (CONSENSUS.PY)")

    @timing
    def remap_flat_consensus(self, org_name):
        '''Remap reads to flattened consensus, save'''
        shell(
            f"./bwa-mem2-2.2.1_x64-linux/bwa-mem2 index {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta")
        shell(f"samtools fastq {self.a['folder_stem']}consensus_data/{org_name}/collated_reads.bam | ./bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta - | viral_consensus -i - -r {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta -o {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_remapped_consensus_sequence.fasta --min_depth {self.a['ConsensusMinD']}")
        shell(
            f"find {self.a['folder_stem']}consensus_data/{org_name}/ -name '*.fasta.*' -delete")

    @timing
    def call_ref_corrected_consensus(self, tar_name):
        '''Construct a `conventional` consensus from grouped reads, with reference to a complete reference genome'''
        '''Load grouped aligned first consensus seqs to retrieve each target name'''
        self.a["GtOrg"] = "Paramyxoviridae_RSV"
        self.a["GtFile"] = "data/rsv_set_metadata.csv"

        if tar_name != self.a['GtOrg']:
            print(f"Target {tar_name} is not the GT organism")
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
        shell(
            f"./bwa-mem2-2.2.1_x64-linux/bwa-mem2 index {self.fnames['temp_ref_seq']}")

        '''Run alignment and flatten consensus'''
        shell(
            f"samtools fastq {self.a['folder_stem']}consensus_data/{tar_name}/collated_reads.bam > {self.fnames['collated_reads_fastq']}")
        shell(
            f"./bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem {self.fnames['temp_ref_seq']} {self.fnames['collated_reads_fastq']} | viral_consensus -i - -r {self.fnames['temp_ref_seq']} -o {ref_adj_cons_fname} --min_depth {self.a['ConsensusMinD']}")

        '''Save'''
        shell(f"rm {self.fnames['collated_reads_fastq']}")
        shell(f"rm -r {self.fnames['temp_folder']}")

    def main(self):
        end_sec_print(
            "INFO: Calling consensus sequences\nThis may take a little while...")
        shell(f"samtools index {self.a['folder_stem']}{self.a['SeqName']}.bam",
              "Samtools Index Call (CONSENSUS.PY)")

        for tar_name in os.listdir(f"{self.a['folder_stem']}grouped_reads/"):
            self.filter_bam(tar_name)

        [self.collate_consensus_seqs(tar_name) for tar_name in os.listdir(
            f"{self.a['folder_stem']}/grouped_reads/")]
        [self.call_flat_consensus(i) for i in self.consensus_seqs.keys()]

        # if self.do_ref_adjusted_cons:
        [self.call_ref_corrected_consensus(tar_name)
            for tar_name in self.consensus_seqs.keys()]
        shell(
            f"find {self.a['folder_stem']}grouped_reads/ -name '*.bam' -delete")
        shell(
            f"find {self.a['folder_stem']}consensus_data/ -name '*.bam' -delete")
        end_sec_print("INFO: Consensus calling complete")


if __name__ == "__main__":
    cls = Consensus()
    cls.main()
