import os
from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print


def run_map(p, bwa_path='./bwa-mem2-2.2.1_x64-linux/bwa-mem2', api_entry=True):
    '''Use BWA and Samtools to map reads from each sample to targets'''
    if not api_entry:
        p = {
            "ExpDir": p.ExpDir,
            "SeqName": p.SeqName,
            "RefStem": p.RefStem
        }
    else:
        p["ExpDir"] = f"{p['ExpDir']}/"
    end_sec_print(
        f"INFO: Beginning initial mapping using BWA\nThis may take a while for large files")
    shell(f"{bwa_path} index {p['RefStem']}")
    shell(
        f"{bwa_path} mem -t {os.cpu_count()} {p['RefStem']} experiments/{p['ExpName']}/{p['SeqName']}_[12]_clean.fastq | samtools view -F4 -Sb - | samtools sort - 1> experiments/{p['ExpName']}/{p['SeqName']}.bam")
    shell(f"rm experiments/{p['ExpName']}/{p['SeqName']}_[12]_clean.fastq")
    end_sec_print(f"INFO: BWA mapping complete")
