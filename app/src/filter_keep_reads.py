from __future__ import division
import os
import subprocess as sp
from collections import deque

from app.utils.shell_cmds import loginfo, read_line, shell, logerr
from app.utils.error_handlers import error_handler_filter_keep_reads
from app.utils.system_messages import end_sec_print
from app.utils.utility_fns import enumerate_read_files


class FilterKeepReads:
    '''
    Retain and exclude specific reads from raw seqs according to Kraken annotations
    and optionally via user-specified rules, through comparison with lineage file.
    '''

    def __init__(self, argies) -> None:
        '''Convert API arguments to format of argparser'''
        self.a = argies
        self.a["input_file"] = enumerate_read_files(argies["ExpDir"])
        self.a["kraken"] = f"{argies['SaveDir']}/{argies['ExpName']}/{argies['ExpName']}.kraken"
        '''Run error handler, build output fnames, extend retain/exclude IDs from names'''
        if self.a["DoKrakenPrefilter"]:
            self.a["o"], self.a["ExcludeIds"], self.a["RetainIds"] = error_handler_filter_keep_reads(
                self.a)
        '''Empty vars'''
        self.reads_to_exclude, self.reads_to_keep = deque(), deque()

    def cmd_string(self, inpath):
        '''Expand file shell command depending on input file extension'''
        return ['gunzip', "-c"] if inpath.endswith('.gz') else ['cat']

    def filter_Reads(self, inpath, out_h):
        '''
        Collect the four lines belonging to the next read and filter on readname.
        ASSUMES NO BLANK LINES IN INPUT.
        '''
        handle = sp.Popen([i for i in self.cmd_string(inpath)] + [inpath],
                          bufsize=8192, stdout=sp.PIPE).stdout
        line = read_line(handle)
        num_reads = 0
        while line:
            readheader, seq, qh, qual = line, read_line(
                handle), read_line(handle), read_line(handle)
            line = read_line(handle)
            readname = readheader[1:].split('/')[0].split()[0]
            to_keep = True if not self.reads_to_keep or (
                self.reads_to_keep and (readname in self.reads_to_keep)) else False
            to_exclude = True if (self.reads_to_exclude and (
                readname in self.reads_to_exclude)) else False
            if to_keep and not to_exclude:
                out_h.write(''.join((readheader, seq, qh, qual)))
                num_reads += 1
        return num_reads

    def main(self):
        '''Entrypoint'''
        if self.a["DoKrakenPrefilter"]:
            '''Read in Kraken file, populate exclude/retain ID loc sets'''
            end_sec_print("INFO: Filtering reads using Kraken2 annotations.")
            handle = sp.Popen([i for i in self.cmd_string(
                self.a['kraken'])] + [self.a['kraken']], bufsize=8192, stdout=sp.PIPE).stdout
            line = read_line(handle)
            while line:
                readname, taxid = line.split()[1:3]
                if taxid in self.a["ExcludeIds"]:
                    self.reads_to_exclude.append(readname)
                elif taxid in self.a["RetainIds"]:
                    self.reads_to_keep.append(readname)
                line = read_line(handle)
            self.reads_to_exclude, self.reads_to_keep = frozenset(
                self.reads_to_exclude), frozenset(self.reads_to_keep)
            loginfo(
                f'Excluding {len(self.reads_to_exclude)} reads ({len(self.a["ExcludeIds"])} taxa specified).')
            loginfo(
                f'Retaining {len(self.reads_to_keep) if self.reads_to_keep else "all other"} reads ({len(self.a["RetainIds"])} taxa specified).')

        else:
            '''If not doing kraken filter, just copy and decompress the files'''
            logerr("Skipping Kraken Prefilter as you specified to.")
            self.reads_to_exclude = self.reads_to_keep = ()
            self.a['o'] = [
                f'{self.a["SaveDir"]}/{self.a["ExpName"]}/{self.a["ExpName"]}_{i}_filt.fastq' for i in range(1, 3)]

        '''Iterate over input file, filter reads by excl and retain rules, save output'''
        for (inpath, outpath) in zip(self.a["input_file"], self.a["o"]):
            with open(outpath, 'w', buffering=8192) as out_h:
                num_reads = self.filter_Reads(inpath, out_h)
            loginfo(f'Wrote {num_reads} reads to {outpath}.')

        shell(
            f"rm {self.a['SaveDir']}/{self.a['ExpName']}/{self.a['ExpName']}.kraken")
        end_sec_print("INFO: Filter reads complete.")
