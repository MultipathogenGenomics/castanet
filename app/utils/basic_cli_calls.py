import subprocess as sp

from app.utils.shell_cmds import shell


def samtools_index(fpath):
    shell(f"samtools index {fpath}", "Samtools Index Call")


def samtools_read_num(outdir, seqname, flag=""):
    '''Retrieve total n reads from master BAM file'''
    # Add '-F 260' to switch to only primary aligned mapped reads
    res = sp.run(
        f"samtools view -c {flag} {outdir}/{seqname}.bam", shell=True, capture_output=True)
    return int(res.stdout.decode("utf-8").replace("\n", ""))


def bwa_index(fpath):
    shell(f"./bwa-mem2-2.2.1_x64-linux/bwa-mem2 index {fpath}")


def bam_to_fastq(bam, fastq):
    shell(f"samtools fastq {bam} > {fastq}")


def find_and_delete(fpath, pat):
    shell(f"find {fpath} -name '{pat}' -delete")


def rm(delete, flag=""):
    shell(f"rm {flag} {delete}")
