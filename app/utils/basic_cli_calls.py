import os
import pickle
import subprocess as sp

from app.utils.shell_cmds import shell, loginfo


def samtools_index(fpath):
    shell(f"samtools index {fpath}", "Samtools Index Call")


def get_read_num(args, bam_fname):
    '''Get read num from pre-mapping stage if possible, else default to bam'''
    pickles = [i for i in os.listdir(
        f"{args['SaveDir']}/{args['ExpName']}/") if i[-12:] == "rawreadnum.p"]
    if len(pickles) == 1:
        loginfo(f"Retrieving raw read numbers from trimmed fastq")
        read_num = pickle.load(
            open(f"{args['SaveDir']}/{args['ExpName']}/{pickles[0]}", "rb"))
    else:
        loginfo(f"Retrieving read numbers from bam (fastq was not found; this might be because of your choice of pipeline)")
        read_num = samtools_read_num(bam_fname)
    return read_num


def samtools_read_num(bam_name, flag=""):
    '''Retrieve total n reads from master BAM file'''
    # Add '-F 260' to switch to only primary aligned mapped reads
    res = sp.run(
        f"samtools view -c {flag} {bam_name}", shell=True, capture_output=True)
    return int(res.stdout.decode("utf-8").replace("\n", ""))


def bwa_index(fpath):
    shell(f"bwa-mem2 index {fpath}")


def bam_to_fastq(bam, fastq):
    shell(f"samtools fastq {bam} > {fastq}")


def find_and_delete(fpath, pat):
    shell(f"find {fpath} -name '{pat}' -delete")


def rm(delete, flag=""):
    shell(f"rm {flag} {delete}")
