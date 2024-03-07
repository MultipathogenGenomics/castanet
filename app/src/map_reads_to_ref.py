import os
from app.utils.utility_fns import enumerate_read_files
from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print
from app.utils.shell_cmds import stoperr


def run_map(p):
    '''Use BWA and Samtools to map reads from each sample to targets'''
    '''Default in_files are created by the trimming step'''
    in_files = [f"experiments/{p['ExpName']}/{p['ExpName']}_1_clean.fastq",
                f"experiments/{p['ExpName']}/{p['ExpName']}_2_clean.fastq"]
    CLEAN_UP = True
    for fn in in_files:
        '''If default infiles not present, look in ExpDir for user-specified ones'''
        if not os.path.exists(fn):
            in_files = enumerate_read_files(f"{p['ExpDir']}/")
            CLEAN_UP = False
    for fn in in_files:
        if os.stat(fn).st_size < 2:
            stoperr(
                f"Castanet found an input file: {fn}, but it's empty. Please check your input file have been processed appropriately for input to BWA-mem2.")

    end_sec_print(
        f"INFO: Beginning initial mapping using BWA\nThis may take a while for large files")
    shell(f"bwa-mem2 index {p['RefStem']}")
    shell(
        f"bwa-mem2 mem -t {os.cpu_count()} {p['RefStem']} {in_files[0]} {in_files[1]} | samtools view -F4 -Sb - | samtools sort - 1> experiments/{p['ExpName']}/{p['ExpName']}.bam")
    if os.stat(f"experiments/{p['ExpName']}/{p['ExpName']}.bam").st_size < 2:
        stoperr(f"BWA-MEM2 produced an empty BAM file. Check your BWA-MEM2 installation and that your input reads are of sufficent quality. This might also indicate an out of memory error, if you're crunching a huge dataset.")
    if CLEAN_UP:
        shell(f"rm experiments/{p['ExpName']}/{p['ExpName']}_[12]_clean.fastq")
    end_sec_print(f"INFO: BWA mapping complete")
