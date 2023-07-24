import os
import re
import pysam
import numpy as np
import pandas as pd

from app.utils.timer import timing
from app.utils.shell_cmds import shell, make_dir
from app.utils.utility_fns import get_gene_orgid, read_fa
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
        make_dir(f"mkdir {self.a['folder_stem']}consensus_data/")

    def filter_bam(self, tar_name) -> None:
        '''Take list of QNAME ids, filter and make new bam specific to target'''
        fq = set([i.replace("\n", "") for i in open(
            f"{self.a['folder_stem']}grouped_reads/{tar_name}/{tar_name}.lst").readlines()])
        infile = pysam.AlignmentFile(
            f"{self.a['folder_stem']}{self.a['SeqName']}.bam")
        outfile = pysam.AlignmentFile(
            f"{self.a['folder_stem']}grouped_reads/{tar_name}/{tar_name}.bam", template=infile, mode='wb')
        # RM < TODO CHECK WE'RE NOT OVER MATCHING, e.g. "ERR10812876.61077" > "ERR10812876.610776". Append space?
        [outfile.write(aln) for aln in infile if aln.query_name in fq]
        infile.close(), outfile.close()

    def call_consensus(self, tar_name) -> None:
        '''Call consensus sequence for sam alignment records, grouped by target'''
        # RM < TODO we probably want FQ so we can filter low quality?
        shell(f"samtools consensus -f fasta {self.a['folder_stem']}grouped_reads/{tar_name}/{tar_name}.bam -o {self.a['folder_stem']}grouped_reads/{tar_name}/consensus_seqs_{tar_name}.fasta",
              "Samtools consensus call (CONSENSUS.PY)")

    def collate_consensus_seqs(self, tar_name):
        '''Read and collate to self var the consensus seqs from per target to per organism'''
        seqs_and_refs = [i for i in read_fa(
            f"{self.a['folder_stem']}grouped_reads/{tar_name}/consensus_seqs_{tar_name}.fasta") if tar_name in i[0]]
        # aggn, refn, seq
        seqs_and_refs = [[self.aggregate_to_probename(
            i[0]), i[0], i[1]] for i in seqs_and_refs]

        if len(seqs_and_refs) == 0:
            # RM < TODO TEST PROPERLY with [i[0] for i in seqs_and_refs]
            return
        consensus_org = seqs_and_refs[0][0]

        if not consensus_org in self.consensus_seqs.keys():
            self.consensus_seqs[consensus_org] = [i[2] for i in seqs_and_refs]
            self.consensus_refs[consensus_org] = [
                [i[0], i[1]] for i in seqs_and_refs]
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
            print(
                f"INFO: Only 1 target found for organism: {org_name} (current strategy is to not flatten)")
            flat_consensus = "".join(self.consensus_seqs[org_name])

        else:
            '''Otherwise, flatten with MAFFT'''
            # RM < TODO pad consensus seqs to ~same length?
            '''Retrieve matching ref seqs and save to persistent files'''
            ref_seq_names = list(set([i[0][1].replace(">", "")
                                      for i in self.consensus_refs[org_name]]))
            ref_seqs = [ref for ref in self.refs if ref[0].replace(
                ">", "") in ref_seq_names]
            with open(f"{self.a['folder_stem']}consensus_data/temp_refs.fasta", "w") as f:
                [f.write(f"{i[0]}\n{i[1]}\n") for i in ref_seqs]
            with open(f"{self.a['folder_stem']}consensus_data/temp_seqs.fasta", "w") as f:
                [f.write(f">TARGET_CONSENSUS_{i}\n{self.consensus_seqs[org_name][i]}\n") for i in range(
                    len(self.consensus_seqs[org_name]))]

            flat_consensus = self.call_flattened_consensus(org_name)

        '''Save'''
        with open(f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta", "w") as f:
            f.write(f">{org_name}_consensus\n{flat_consensus}")

    def call_flattened_consensus(self, org_name):
        '''Make MSA of references, then add fragments from target consensuses'''
        print(
            f"INFO: making reference alignments for target group: {org_name}")
        shell(f"mafft --thread -1 {self.a['folder_stem']}consensus_data/temp_refs.fasta > {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_ref_alignment.aln",
              "Mafft align ref seqs (CONSENSUS.PY)")
        print(
            f"INFO: adding consensuses to alignment for organism: {org_name}")
        shell(f"mafft --thread -1 --6merpair --addfragments {self.a['folder_stem']}consensus_data/temp_seqs.fasta {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_ref_alignment.aln > {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln",
              "Mafft align consensus with ref seqs (CONSENSUS.PY)")

        '''Make flat consensus'''
        flat_consensus = self.rich_consensus(np.array([list(i[1]) for i in read_fa(
            f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln")]), False)
        shell(
            f"rm {self.a['folder_stem']}consensus_data/temp_seqs.fasta {self.a['folder_stem']}consensus_data/temp_refs.fasta")
        return flat_consensus

    @timing
    def rich_consensus(self, aln, gap):
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

    def main(self):
        # RM <TODO SWAP MKDIRS FOR UTIL FN
        end_sec_print(
            "INFO: Calling consensus sequences\nThis may take a little while...")
        shell(f"samtools index {self.a['folder_stem']}{self.a['SeqName']}.bam",
              "Samtools Index Call (CONSENSUS.PY)")
        for tar_name in os.listdir(f"{self.a['folder_stem']}grouped_reads/"):
            self.filter_bam(tar_name)
            self.call_consensus(tar_name)

        [self.collate_consensus_seqs(tar_name) for tar_name in os.listdir(
            f"{self.a['folder_stem']}/grouped_reads/")]
        [self.call_flat_consensus(i) for i in self.consensus_seqs.keys()]
        [self.call_ref_corrected_consensus(tar_name)
         for tar_name in self.consensus_seqs.keys()]

    def call_ref_corrected_consensus(self, tar_name):
        '''TODO DOCSTRING'''
        targets_for_group = [i.replace("\n", "").replace(">", "") for i in open(
            f"{self.a['folder_stem']}consensus_data/{tar_name}/{tar_name}_consensus_alignment.aln").readlines() if i.startswith(">") and not "TARGET" in i]
        target_dirs = [i for i in os.listdir(
            f"{self.a['folder_stem']}grouped_reads/") if i in targets_for_group]
        bamout_fname = f"{self.a['folder_stem']}consensus_data/ref_cons_merge.bam"
        # can't use both types of quote mark inside an fstring
        fol_stem = self.a['folder_stem']

        shell(
            f"""samtools merge {bamout_fname} {' '.join([f'{fol_stem}grouped_reads/{i}/{i}.bam' for i in target_dirs])}""")

        # Maybe some kallisto in the middle to select ref seq?
        ref_seq = "LR699734.1"  # RM < TODO HARD CODED FOR EX
        # << MAKE BLAST DB https://ftp.ncbi.nlm.nih.gov/blast/db/ https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html
        bam_name = f"{self.a['folder_stem']}consensus_data/temp_bam.bam"
        fas_name = f"{self.a['folder_stem']}consensus_data/temp_fasta.fasta"

        # samtools consensus with ref in
        shell(
            f"samtools consensus -a --show-ins no {bamout_fname} -o {self.a['folder_stem']}consensus_data/ref_cons_merge.fasta")

        contigs = [i for i in read_fa(
            f"{self.a['folder_stem']}consensus_data/ref_cons_merge.fasta") if i[0].replace(">", "") in targets_for_group]
        contigs.append(
            read_fa(f"./pathogen_seq_dbs/{tar_name}/{ref_seq}/{ref_seq}.fasta")[0])
        with open(f"{self.a['folder_stem']}consensus_data/ref_cons_merge.fasta", "w") as f:
            [f.write(f"{i[0]}\n{i[1]}\n") for i in contigs]

        print(
            f"INFO: generating reference-adjusted consensus for target group / reference: {tar_name} / {ref_seq}")
        shell(
            f"./bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem ./pathogen_seq_dbs/{tar_name}/{ref_seq}/{ref_seq}.fasta {self.a['folder_stem']}consensus_data/ref_cons_merge.fasta > {bam_name} && samtools fasta {bam_name} > {fas_name}")
        flat_consensus = self.rich_consensus(np.array([[i[1]] for i in read_fa(
            fas_name) if i[0].replace(">", "") in targets_for_group]), False)

        '''Save'''
        with open(f"{self.a['folder_stem']}consensus_data/{tar_name}/{tar_name}_ref_adjusted_consensus_sequence.fasta", "w") as f:
            f.write(f">{tar_name}_consensus\n{flat_consensus}")
        shell(
            f"rm {bam_name} {fas_name} {self.a['folder_stem']}consensus_data/ref_cons_merge.bam {self.a['folder_stem']}consensus_data/ref_cons_merge.fasta")


if __name__ == "__main__":
    cls = Consensus()
    cls.main()
