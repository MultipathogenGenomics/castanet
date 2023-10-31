from app.utils.shell_cmds import shell, stoperr
from app.utils.utility_fns import read_fa
from app.utils.system_messages import end_sec_print


def run_trim(p, trim_path='java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar', api_entry=True):
    '''Call BWA trim CLI tool; check it worked; remove interim files'''
    if not api_entry:
        p = {
            "ExpDir": p.ExpDir,
            "SeqName": p.SeqName,
            "NThreads": 4,
            "AdaptP": p.AdaptP
        }
    else:
        p["ExpDir"] = f"{p['ExpDir']}/"
    shell(f"{trim_path} PE -threads {p['NThreads']} experiments/{p['ExpName']}/{p['SeqName']}_1_filt.fastq experiments/{p['ExpName']}/{p['SeqName']}_2_filt.fastq experiments/{p['ExpName']}/{p['SeqName']}_1_clean.fastq experiments/{p['ExpName']}/{p['SeqName']}_1_trimmings.fq experiments/{p['ExpName']}/{p['SeqName']}_2_clean.fastq experiments/{p['ExpName']}/{p['SeqName']}_2_trimmings.fq ILLUMINACLIP:{p['AdaptP']}:2:10:7:1:true MINLEN:{p['TrimMinLen']}")  # was 80
    reads_1, reads_2 = read_fa(f"experiments/{p['ExpName']}/{p['SeqName']}_1_clean.fastq"), read_fa(
        f"experiments/{p['ExpName']}/{p['SeqName']}_2_clean.fastq")
    if len(reads_1) == 0 and len(reads_2) == 0:
        stoperr(f"Trimming produced empty files. Check your TrimMinLen parameter is not too short for your sequences and that Trimmomatic is isntalled.")
    shell(
        f"rm experiments/{p['ExpName']}/{p['SeqName']}_[12]_filt.fastq experiments/{p['ExpName']}/{p['SeqName']}_[12]_trimmings.fq")
    end_sec_print("INFO: Read Trimming complete.")
