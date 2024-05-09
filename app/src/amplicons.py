import os
import re
import numpy as np
import pandas as pd
from app.utils.shell_cmds import shell, loginfo, logerr
from app.utils.utility_fns import read_fa

'''
REQS: BAM and pair of FASTQs in dir, activate Castanet env
Get TSV of a bam, pull out read numbers, extract same reads from original fqs.
Purpose is to filter original FASTQ by only reads that map to something in the wide panel.
These can then be re-mapped to a more specific panel.
'''


class Amplicons:
    def __init__(self, payload) -> None:
        self.a = payload
        self.min_amp_len = 40  # RM < TODO Parameterise
        self.do_aln_graphs = True  # RM < TODO Parameterise
        self.results = {}
        self.bam_fname = f"{self.a['ExpRoot']}/{self.a['ExpName']}/{self.a['ExpName']}.bam"
        self.amp_folder = f"{self.a['ExpRoot']}/{self.a['ExpName']}/amplicon_data/"
        if not os.path.exists(self.amp_folder):
            os.mkdir(self.amp_folder)

    def get_tsvs(self):
        '''Convert bam to readable TSV format and input FASTQ to FASTA'''
        tsv_fname = self.bam_fname.replace(".bam", ".tsv")
        # RM TODO < TEST OUTPUT
        shell(f"samtools view -F 0x40 {self.bam_fname} > {tsv_fname}")

        '''Read in BAM view TSV and FASTAs'''
        try:
            df = pd.read_csv(tsv_fname, sep="\t", header=None,
                             error_bad_lines=False)
        except:
            print(f"No file to read for {self.bam_fname}")
            return

        # [seqid, refname, cigar, seq]
        seqs = df[[0, 2, 5, 9]].values.tolist()
        read_stats = df[2].value_counts().to_dict()

        return seqs, read_stats

    def crunch(self, seqs):
        for seq in seqs:
            self.filt(seq)

    def remove_char(self, input_string, index):
        first_part = input_string[:index]
        second_part = input_string[index+1:]
        return first_part + second_part

    def filt(self, row):
        cigar = row[2]
        seq = str(row[3])
        cig_parts = [i for i in re.split(
            '([0-9]+[A-Z])', cigar) if not i == ""]
        cnt = 0
        kill_indices = []
        for part in cig_parts:
            part_len = int(re.split(r'[A-Z]', part)[0])
            if "S" in part:
                kill_indices.append(np.arange(cnt, cnt + part_len))
            elif "D" in part:
                continue
            cnt = cnt + part_len

        ar = np.zeros((len(seq),), dtype=str)
        for idx, val in enumerate(seq):
            ar[idx] = val
        if len(kill_indices) > 1:
            kill_indices = np.concatenate([i for i in kill_indices])
        final_seq = "".join(np.delete(ar, kill_indices))
        if len(final_seq) > self.min_amp_len:
            if not row[1] in self.results.keys():
                self.results[row[1]] = []
            self.results[row[1]].append([f"{row[1]}_{row[0]}", final_seq])

    def stats(self, read_stats):
        '''Save CSV with details of all and unique reads'''
        dedupe_stats = {}
        for ref in self.results.keys():
            seen = set()
            for seq in self.results[ref]:
                seen.add(seq[1])
            dedupe_stats[ref] = len(seen)

        all_stats = {}
        # ref: [all reads, dedupe reads]
        for ref in read_stats.keys():
            all_stats[ref] = [read_stats[ref], dedupe_stats[ref]]
            loginfo(
                f" - Target: {ref}. Total reads: {read_stats[ref]}; Deduplicated reads: {dedupe_stats[ref]}")
        df = pd.DataFrame.from_dict(all_stats).T
        df = df.rename(columns={0: "total_reads", 1: "dedup_reads"})
        df.to_csv(f"{self.amp_folder}/read_statistics.csv")

    def save(self):
        '''Save DEDUPLICATED reads, separately grouped by target'''
        for ref in self.results.keys():
            with open(f"{self.amp_folder}/{ref}.fasta", "w") as f:
                seen = set()
                for idx, seq in enumerate(self.results[ref]):
                    if not seq[1] in seen:
                        f.write(
                            f">{self.a['ExpName']}_{idx} {seq[0]}\n{seq[1]}\n")
                    seen.add(seq[1])
            if self.do_aln_graphs:
                loginfo(
                    f"Generating alignment and graph for ref {ref}. This might take a few moments.")
                shell(
                    f"mafft --auto --thread {self.a['NThreads']} {self.amp_folder}/{ref}.fasta > {self.amp_folder}/{ref}.aln")
                # RM <TODO  add lib to requirements, swap CLI to API
                out = shell(
                    f"pymsaviz -i {self.amp_folder}/{ref}.aln -o {self.amp_folder}/{ref}.png --color_scheme Identity --show_consensus --show_grid", is_test=True)
                if "ValueError" in out:
                    logerr(
                        f"I couldn't produce an alignment plot for ref {ref}, this usually happens if you have so many unique amplicons to align that the plot would be insanely huge.")

    def clean(self):
        '''Delete TSV file when done'''
        shell(f"rm {self.bam_fname.replace('.bam','.tsv')}")

    def main(self):
        loginfo("Running amplicon analysis")
        seqs, read_stats = self.get_tsvs()
        self.crunch(seqs)
        self.stats(read_stats)
        self.save()
        self.clean()
