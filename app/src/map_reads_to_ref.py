import os
from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print
from app.utils.shell_cmds import stoperr


def run_map(p):
    '''Use BWA and Samtools to map reads from each sample to targets'''
    p["ExpDir"] = f"{p['ExpDir']}/"
    end_sec_print(
        f"INFO: Beginning initial mapping using BWA\nThis may take a while for large files")
    shell(f"bwa-mem2 index {p['RefStem']}")
    shell(
        f"bwa-mem2 mem -t {os.cpu_count()} {p['RefStem']} experiments/{p['ExpName']}/{p['ExpName']}_[12]_clean.fastq | samtools view -F4 -Sb - | samtools sort - 1> experiments/{p['ExpName']}/{p['ExpName']}.bam")
    if os.stat(f"experiments/{p['ExpName']}/{p['ExpName']}.bam").st_size < 2:
        stoperr(f"BWA-MEM2 produced an empty BAM file. Check your BWA-MEM2 installation and that your input reads are of sufficent quality. This might also indicate an out of memory error, if you're crunching a huge dataset.")
    shell(f"rm experiments/{p['ExpName']}/{p['ExpName']}_[12]_clean.fastq")
    end_sec_print(f"INFO: BWA mapping complete")
