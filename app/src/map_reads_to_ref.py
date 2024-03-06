import os
from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print


def run_map(p):
    '''Use BWA and Samtools to map reads from each sample to targets'''
    p["ExpDir"] = f"{p['ExpDir']}/"
    end_sec_print(
        f"INFO: Beginning initial mapping using BWA\nThis may take a while for large files")
    shell(f"bwa-mem2 index {p['RefStem']}")
    shell(
        f"bwa-mem2 mem -t {os.cpu_count()} {p['RefStem']} experiments/{p['ExpName']}/{p['ExpName']}_[12]_clean.fastq | samtools view -F4 -Sb - | samtools sort - 1> experiments/{p['ExpName']}/{p['ExpName']}.bam")
    shell(f"rm experiments/{p['ExpName']}/{p['ExpName']}_[12]_clean.fastq")
    end_sec_print(f"INFO: BWA mapping complete")
