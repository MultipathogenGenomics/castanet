from __future__ import print_function
import sys
import re
import pandas as pd

from app.utils.argparsers import parse_args_bam_parse
from app.utils.error_handlers import error_handler_parse_bam_positions
from app.utils.shell_cmds import make_dir
from app.utils.utility_fns import get_gene_orgid, trim_long_fpaths


class Parse_bam_positions:
    '''
    Parse contents of bam file, via reading shell commands passed in.
    Count reads mapped to each target, generate read groupings for consensus calling.
    Should only get called by another Python script due to requirement for shell input.
    '''

    def __init__(self, argies) -> None:
        '''N.b. argies are not a dict, but a namespace!!'''
        self.min_match_length = 40
        self.argies = argies
        self.n = 10  # Min n reads to decide we want to make a consensus
        self.reads_by_hit = {}

    def getmatchsize(self, cigar):
        '''Find matches in cigar string with regex, return count'''
        matches = re.findall(r'([0-9]+)M', cigar)
        if not len(matches):
            return 0
        return sum(int(x) for x in matches)

    def build_target_dbs(self, ref, seq, id):
        if not ref in self.reads_by_hit.keys():
            self.reads_by_hit[ref] = [[id, seq]]
        else:
            self.reads_by_hit[ref].append([id, seq])

    def parse_bam_position(self, l):
        '''GENERATE COUNTS STAGE: For each line passed in, parse fields of interest; identify matches and print back to stdout.'''
        fields = l.split()
        id, ref, pos, ref2, tlen, cigar, seq = fields[0], fields[2], fields[3], fields[6], int(
            fields[8]), fields[5], fields[9]

        match = tlen >= self.min_match_length
        if not self.argies.SingleEnded == "True":
            improper_match = (tlen == 0) and (self.getmatchsize(
                cigar) >= self.min_match_length) and (get_gene_orgid(ref) == get_gene_orgid(ref2))
        else:
            '''RM < Experimental, for use with single ended sets (e.g. when Sequencer explodes mid-run)'''
            improper_match = (tlen == 0) and (self.getmatchsize(
                cigar) >= self.min_match_length)
            try:
                tlen = self.getmatchsize(cigar)
            except:
                tlen = 0

        if match or improper_match:
            '''Properly paired and match is of decent mapped length OR
            Improperly paired BUT same gene AND match is of decent mapped length (via CIGAR string lookup) AND RNAME ref organism is same to RNEXT ref org'''
            print(f'{ref},{pos},{tlen},{self.argies.SeqName}')
            self.build_target_dbs(ref, seq, id)
        else:
            return

    def filter_bam(self, l, reads_to_drop):
        '''POST FILTER STAGE'''
        if l.startswith('@'):
            sys.stdout.write(l)
            return
        if len(reads_to_drop):
            fields = l.split()
            ref = fields[2]
            pos = int(fields[3])
            tlen = int(fields[8])
            '''output only reads that are not found in the indexed READS_TO_DROP file'''
            if not (ref, pos, tlen) in reads_to_drop.index:
                sys.stdout.write(l)
            else:
                return

    def save_hit_dbs(self):
        grp_aln_f = f"experiments/{self.argies.ExpName}/grouped_reads/"
        make_dir(f"mkdir {grp_aln_f}")

        for key in self.reads_by_hit.keys():
            '''Save list of grouped read QNAME ids for calling consensuses, if more than n reads'''
            if len(self.reads_by_hit[key]) < self.n:
                '''Don't save grouped reads if less than n'''
                continue

            short_key = trim_long_fpaths(key)
            make_dir(f"mkdir {grp_aln_f}{short_key}")
            with open(f"{grp_aln_f}{short_key}/{short_key}.lst", "w") as file:
                [file.write(f"{self.reads_by_hit[key][i][0]}\n")
                    for i in range(len(self.reads_by_hit[key]))]

    def main(self):
        '''Entrypoint. Multi functional across generate counts and post filter.'''
        error_handler_parse_bam_positions(sys.argv)
        if self.argies.Mode == "filter":
            try:
                reads_to_drop = pd.read_csv(
                    self.argies.ExpDir + self.argies.FilterFile)
            except:
                raise SystemExit(
                    f"Couldn't find reads to drop file: {self.argies.FilterFile}. Did your run generate one?")

        with open(f"experiments/{self.argies.ExpName}/{self.argies.SeqName}_bamview.txt") as f:
            for l in f:
                if self.argies.Mode == "parse":
                    self.parse_bam_position(l)
                elif self.argies.Mode == "reparse":
                    # RM < TODO deprecate
                    self.parse_bam_position(l)
                else:
                    self.filter_bam(l, reads_to_drop)

        if len(self.reads_by_hit) == 0:
            '''No hits found between input BAM file and reference sequences'''
            return
        else:
            '''Save data for consensus call fns'''
            self.save_hit_dbs()


if __name__ == '__main__':
    cls = Parse_bam_positions(parse_args_bam_parse())
    cls.main()
