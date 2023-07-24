import subprocess as sp
from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print


def run_counts(p, api_entry=True):
    '''Pipe mapped bam into parse functions, generate counts and consensus groupings'''
    if not api_entry:
        p = {
            "ExpDir": p.ExpDir,
            "SeqName": p.SeqName,
        }
    else:
        p["ExpDir"] = f"{p['ExpDir']}/"
    bamview_fname = f"experiments/{p['ExpName']}/{p['SeqName']}_bamview.txt"

    '''Split samtools pipe to python with persistent file for ease of debugging'''
    end_sec_print("Info: Generating read counts ")
    shell(
        f"""samtools view -F2048 -F4 experiments/{p['ExpName']}/{p['SeqName']}.bam > {bamview_fname}""")
    sp.run(  # Use subprocess run rather than Popen as complex call is a PITA
        f"python3 -m app.src.parse_bam -Mode parse -SeqName {p['SeqName']} -ExpDir {p['ExpDir']} -ExpName {p['ExpName']} | sort | uniq -c | sed s'/ /,/'g | sed s'/^[,]*//'g > experiments/{p['ExpName']}/{p['SeqName']}_PosCounts.csv", shell=True)
    shell(f"rm {bamview_fname}")
    end_sec_print("INFO: Counts generated")
