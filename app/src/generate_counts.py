import os
import subprocess as sp
from app.utils.utility_fns import enumerate_bam_files
from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print
from app.utils.shell_cmds import stoperr


def run_counts(p):
    '''Pipe mapped bam into parse functions, generate counts and consensus groupings'''
    '''Default BAM name, propagated in end to end function'''
    in_file = f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam"
    bamview_fname = f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_bamview.txt"
    if not os.path.exists(in_file):
        '''If firing in isolation, look in experiment directory for a bam file'''
        in_file = enumerate_bam_files(f'{p["ExpDir"]}/')
    if os.stat(in_file).st_size < 2:
        stoperr(
            f"Your input BAM file ({in_file}) is empty. Please review how it was generated, i.e. was your mapping successful.")

    '''Split samtools pipe to python with persistent file for ease of debugging'''
    end_sec_print("Info: Generating read counts ")
    shell(
        f"""samtools view -F2048 -F4 {in_file} > {bamview_fname}""")
    sp.run(  # Use subprocess run rather than Popen as complex call is a PITA
        f"python3 -m app.src.parse_bam -Mode parse -SeqName {p['ExpName']} -ExpDir {p['ExpDir']}/ -ExpName {p['ExpName']} -SingleEnded {p['SingleEndedReads']} -SaveDir {p['SaveDir']} -MatchLength {p['MatchLength']} | sort | uniq -c | sed s'/ /,/'g | sed s'/^[,]*//'g > {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_PosCounts.csv", shell=True)
    shell(f"rm {bamview_fname}")
    end_sec_print("INFO: Counts generated")
