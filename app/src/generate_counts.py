import subprocess as sp
from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print


def run_counts(p):
    '''Pipe mapped bam into parse functions, generate counts and consensus groupings'''
    p["ExpDir"] = f"{p['ExpDir']}/"
    bamview_fname = f"experiments/{p['ExpName']}/{p['ExpName']}_bamview.txt"
    single_ended = p["SingleEndedReads"]

    '''Split samtools pipe to python with persistent file for ease of debugging'''
    end_sec_print("Info: Generating read counts ")
    shell(
        f"""samtools view -F2048 -F4 experiments/{p['ExpName']}/{p['ExpName']}.bam > {bamview_fname}""")
    sp.run(  # Use subprocess run rather than Popen as complex call is a PITA
        f"python3 -m app.src.parse_bam -Mode parse -SeqName {p['ExpName']} -ExpDir {p['ExpDir']} -ExpName {p['ExpName']} -SingleEnded {single_ended} | sort | uniq -c | sed s'/ /,/'g | sed s'/^[,]*//'g > experiments/{p['ExpName']}/{p['ExpName']}_PosCounts.csv", shell=True)
    shell(f"rm {bamview_fname}")
    end_sec_print("INFO: Counts generated")
