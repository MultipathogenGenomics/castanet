import os
import re
import pysam
import numpy as np

from app.utils.timer import timing
from app.utils.shell_cmds import shell
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
        if not os.path.isdir(f"{self.a['folder_stem']}consensus_data/"):
            shell(f"mkdir {self.a['folder_stem']}consensus_data/")

    def filter_bam(self, tar_name) -> None:
        '''Take list of QNAME ids, filter and make new bam specific to target'''
        fq = set([i.replace("\n", "") for i in open(
            f"{self.a['folder_stem']}grouped_reads/{tar_name}/{tar_name}.lst").readlines()])
        infile = pysam.AlignmentFile(
            f"{self.a['folder_stem']}{self.a['SeqName']}.bam")
        outfile = pysam.AlignmentFile(
            f"{self.a['folder_stem']}grouped_reads/{tar_name}/{tar_name}.bam", template=infile, mode='wb')
        [outfile.write(aln) for aln in infile if aln.query_name in fq]
        infile.close(), outfile.close()

    def call_consensus(self, tar_name) -> None:
        '''Call consensus sequence for sam alignment records, grouped by target'''
        # RM < TODO we probably want FQ so we can filter low quality?
        shell(f"samtools consensus -f fasta {self.a['folder_stem']}grouped_reads/{tar_name}/{tar_name}.bam -o {self.a['folder_stem']}grouped_reads/{tar_name}/consensus_seqs_{tar_name}.fasta",
              "Samtools consensus call (CONSENSUS.PY)")

    def collate_consensus_seqs(self, tar_name, consensus_org):
        '''Read and collate to self var the consensus seqs from per target to per organism'''
        seqs_and_refs = [i for i in read_fa(
            f"{self.a['folder_stem']}grouped_reads/{tar_name}/consensus_seqs_{tar_name}.fasta") if tar_name in i[0]]
        if not consensus_org in self.consensus_seqs.keys():
            self.consensus_seqs[consensus_org] = [i[1] for i in seqs_and_refs]
            self.consensus_refs[consensus_org] = [i[0] for i in seqs_and_refs]
        else:
            self.consensus_seqs[consensus_org].append(
                ', '.join([i[1] for i in seqs_and_refs]))
            self.consensus_refs[consensus_org].append(
                ', '.join([i[0] for i in seqs_and_refs]))

    def flatten_consensus(self, org_name):
        '''Create flat consensus for an organism'''
        if not os.path.isdir(f"{self.a['folder_stem']}consensus_data/{org_name}/"):
            shell(f"mkdir {self.a['folder_stem']}consensus_data/{org_name}/")

        if len(self.consensus_seqs[org_name]) == 1:
            '''No flattening to be done, just a single target for this organism'''
            print(
                f"INFO: Only 1 target found for organism: {org_name} (current strategy is to not flatten)")
            flat_consensus = "".join(self.consensus_seqs[org_name])

        else:
            '''Otherwise, flatten with MAFFT'''
            # RM < TODO pad?
            '''Retrieve matching ref seqs and save to persistent files'''
            ref_seq_names = list(set([i.replace(">", "")
                                 for i in self.consensus_refs[org_name]]))
            ref_seqs = [ref for ref in self.refs if any(
                ref_seq_name in ref[0] for ref_seq_name in ref_seq_names)]
            with open(f"{self.a['folder_stem']}consensus_data/temp_refs.fasta", "w") as f:
                [f.write(f"{i[0]}\n{i[1]}\n") for i in ref_seqs]
            with open(f"{self.a['folder_stem']}consensus_data/temp_seqs.fasta", "w") as f:
                [f.write(f">TARGET_CONSENSUS_{i}\n{self.consensus_seqs[org_name][i]}\n") for i in range(
                    len(self.consensus_seqs[org_name]))]

            '''Make MSA of references, then add fragments from target consensuses'''
            print(
                f"INFO: making reference alignments for target group: {org_name}")
            shell(f"mafft {self.a['folder_stem']}consensus_data/temp_refs.fasta > {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_ref_alignment.aln",
                  "Mafft align ref seqs (CONSENSUS.PY)")
            print(
                f"INFO: adding consensuses to alignment for organism: {org_name}")
            shell(f"mafft --6merpair --addfragments {self.a['folder_stem']}consensus_data/temp_seqs.fasta {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_ref_alignment.aln > {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln",
                  "Mafft align consensus with ref seqs (CONSENSUS.PY)")

            '''Make flat consensus'''
            flat_consensus = self.rich_consensus(np.array([list(i[1]) for i in read_fa(
                f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln")]), False)
            shell(
                f"rm {self.a['folder_stem']}consensus_data/temp_seqs.fasta {self.a['folder_stem']}consensus_data/temp_refs.fasta")

        '''Save'''
        with open(f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_sequence.fasta", "w") as f:
            f.write(f">{org_name}_consensus\n{flat_consensus}")

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
        end_sec_print("Calling consensus sequences")
        shell(f"samtools index {self.a['folder_stem']}{self.a['SeqName']}.bam",
              "Samtools Index Call (CONSENSUS.PY)")
        for tar_name in os.listdir(f"{self.a['folder_stem']}grouped_reads/"):
            self.filter_bam(tar_name)
            self.call_consensus(tar_name)

        [self.collate_consensus_seqs(tar_name, get_gene_orgid(tar_name)[
                                     0]) for tar_name in os.listdir(f"{self.a['folder_stem']}/grouped_reads/")]
        [self.flatten_consensus(i) for i in self.consensus_seqs.keys()]


class Consensus_maptorefs:
    def __init__(self, SeqName, RefPathogen, GroupedAlignmentsFname) -> None:
        self.Bwa = './bwa-mem2-2.2.1_x64-linux/bwa-mem2'
        self.SeqName = SeqName
        self.RefPathogen = RefPathogen
        self.PathogenSeqDbs = "./pathogen_seq_dbs/"
        self.GroupedAlignmentsFname = GroupedAlignmentsFname
        self.GroupedAlignmentsFpath = f'experiments/{SeqName}/grouped_alignments/'
        self.RemapDir = f'{self.GroupedAlignmentsFpath}remaps/'
        self.ConsensusDir = f"{self.GroupedAlignmentsFpath}consensus_seqs/"
        self.ConsensusCandidates = f"{self.ConsensusDir}consensus_candidates_{self.RefPathogen}.fasta"
        if not os.path.isdir(self.PathogenSeqDbs):
            shell(f"mkdir {self.PathogenSeqDbs}")
        if not os.path.isdir(self.RemapDir):
            shell(f"mkdir {self.RemapDir}")
        if not os.path.isdir(self.ConsensusDir):
            shell(f"mkdir {self.ConsensusDir}")
            shell(f"touch {self.ConsensusCandidates}")
        else:
            shell(f"> {self.ConsensusCandidates}")

    def get_ref_db(self):
        # RM < TODO MAKE DYNAMIC
        return "pathogen_seq_dbs/hbv/all/all.fasta"

    def get_ref_seq(self, refseq):
        # RM < TODO MAKE DYNAMIC
        return f'hbv/{refseq.replace(".fasta", "")}/{refseq}'

    def reparse(self, refdb):
        '''Re-map to bam, stream SAM to new parse functions'''
        shell(f"{self.Bwa} mem {refdb} {self.GroupedAlignmentsFpath}{self.GroupedAlignmentsFname}.fasta | samtools view -F4 -Sb - | samtools sort - | samtools view -F2048 -F4 - > {self.GroupedAlignmentsFpath}{self.GroupedAlignmentsFname}.sam")
        shell(
            f"python3 -m app.src.parse_bam -Mode reparse -SeqName {self.SeqName}",
            "Call to Castanet parse_bam class (CONSENSUS.PY)")

    def make_consensus(self, fname):
        print(f"Generating candidate consensus for: {fname}")
        # Check if we have a db for this pathogen
        refseq = self.get_ref_seq(fname)
        # map vs genome
        shell(f"{self.Bwa} mem {self.PathogenSeqDbs}{refseq} {self.RemapDir}{fname} | samtools view -F4 -Sb - | samtools sort - 1> {self.RemapDir}{fname}.bam",
              "BWA mem on")
        # build consensus
        shell(
            f"viral_consensus -i {self.RemapDir}{fname}.bam -r {self.PathogenSeqDbs}/{refseq} -o - >> {self.ConsensusCandidates}")

    def qa(self):
        # QA Consensus Genomes
        COV_THRESH = 0.8
        good_seqs = {}
        with open(self.ConsensusCandidates, "r") as f:
            for l in f:
                if l[0] == ">":
                    '''Header'''
                    find_acc = re.findall(
                        r"[A-Z]{1,2}[0-9]{5,6}|[A-Z]{4}[0-9]{6,8}|[A-Z]{2}_[0-9]{6}", l)
                    if find_acc:
                        assert len(
                            find_acc) > 1, f"Can't find two Accession IDs in the fasta header for line {l}"
                        assert find_acc[0] == find_acc[
                            1], f"Accession IDs were mismatched between input bam and ref for line {l}"
                        acc = find_acc[0]
                else:
                    '''Seq'''
                    cov = round((len(l) - l.count('N')) / len(l) * 100)
                    if cov * 0.01 > COV_THRESH:
                        good_seqs[acc] = l
                        print(
                            f"Accepting consensus for hit: {acc}, as genome coverage is {cov}%.")
                    else:
                        print(
                            f"Rejecting consensus for hit: {acc}, as genome coverage is {cov}%.")

        for key in good_seqs.keys():
            with open(f"{self.ConsensusDir}/consensus_{key}.fasta", "w") as file:
                file.write(f">{key}_consensus\n{good_seqs[key]}\n")

    def main(self):
        self.reparse(self.get_ref_db())
        if len(os.listdir(self.RemapDir)) == 0:
            print(
                f"Failed to remap grouped alignments for grouping: {self.RefPathogen} to reference genomes")
            return
        [self.make_consensus(i) for i in os.listdir(
            self.RemapDir) if not "bam" in i]
        self.qa()


if __name__ == "__main__":
    # cls = Consensus_maptorefs(SeqName="HBV13_2",
    #                 RefPathogen="HBV",
    #                 GroupedAlignmentsFname="pHBV")
    # cls.main()
    cls = Consensus()
    cls.main()
