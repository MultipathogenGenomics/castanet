from app.utils.shell_cmds import shell


def samtools_index(fpath):
    shell(f"samtools index {fpath}", "Samtools Index Call")


def bwa_index(fpath):
    shell(f"./bwa-mem2-2.2.1_x64-linux/bwa-mem2 index {fpath}")


def bam_to_fastq(bam, fastq):
    shell(f"samtools fastq {bam} > {fastq}")


def find_and_delete(fpath, pat):
    shell(f"find {fpath} -name '{pat}' -delete")


def rm(delete, flag=""):
    (f"rm {flag} {delete}")
