from app.utils.shell_cmds import shell


def samtools_index(fpath):
    shell(f"samtools index {fpath}", "Samtools Index Call")


def bwa_index(fpath):
    shell(f"./bwa-mem2-2.2.1_x64-linux/bwa-mem2 index {fpath}")
