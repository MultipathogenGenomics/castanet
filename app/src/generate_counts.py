import subprocess as sp
from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print


def run_counts(p, api_entry=True):
    if not api_entry:
        p = {
            "ExpDir": p.ExpDir,
            "SeqName": p.SeqName,
        }
    else:
        p["ExpDir"] = f"{p['ExpDir']}/"

    bamview_fname = f"experiments/{p['ExpName']}/{p['SeqName']}_bamview.txt"
    shell(
        f"""samtools view -F2048 -F4 experiments/{p['ExpName']}/{p['SeqName']}.bam > {bamview_fname}""")

    # Easier to run call to Python script via sp run rather than utility fn shell()
    sp.run(
        f"python3 -m app.src.parse_bam -Mode parse -SeqName {p['SeqName']} -ExpDir {p['ExpDir']} -ExpName {p['ExpName']} | sort | uniq -c | sed s'/ /,/'g | sed s'/^[,]*//'g > experiments/{p['ExpName']}/{p['SeqName']}_PosCounts.csv", shell=True)
    end_sec_print("INFO: Counts generated")
