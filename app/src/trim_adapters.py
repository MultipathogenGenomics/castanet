from app.utils.shell_cmds import shell, stoperr
from app.utils.utility_fns import read_fa
from app.utils.system_messages import end_sec_print


def run_trim(p, trim_path='java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar'):
    '''Call Trimmomatic trim CLI tool; check it worked; remove interim files'''
    p["ExpDir"] = f"{p['ExpDir']}/"
    files = {
        "in_files": [f"experiments/{p['ExpName']}/{p['ExpName']}_1_filt.fastq", f"experiments/{p['ExpName']}/{p['ExpName']}_2_filt.fastq"],
        "clean_files": [f"experiments/{p['ExpName']}/{p['ExpName']}_1_clean.fastq", f"experiments/{p['ExpName']}/{p['ExpName']}_2_clean.fastq"],
        "trim_files": [f"experiments/{p['ExpName']}/{p['ExpName']}_1_trimmings.fq", f"experiments/{p['ExpName']}/{p['ExpName']}_2_trimmings.fq"]
    }
    if p["DoTrimming"]:
        end_sec_print("INFO: Read Trimming beginning.")
        shell(f"{trim_path} PE -threads {p['NThreads']} {files['in_files'][0]} {files['in_files'][1]} {files['clean_files'][0]} {files['trim_files'][0]} {files['clean_files'][1]} {files['trim_files'][1]} ILLUMINACLIP:{p['AdaptP']}:2:10:7:1:true MINLEN:{p['TrimMinLen']}")

        '''Test for success'''
        try:
            reads_1, reads_2 = read_fa(files['clean_files'][0]), read_fa(
                files['clean_files'][1])
        except FileNotFoundError:
            stoperr(f"Trimming produced empty files. Check your TrimMinLen parameter is not too short for your sequences and that Trimmomatic is isntalled.")
        if len(reads_1) == 0 and len(reads_2) == 0:
            stoperr(f"Trimming produced empty files. Check your TrimMinLen parameter is not too short for your sequences and that Trimmomatic is isntalled.")
        for idx, fi in enumerate(files['in_files']):
            shell(
                f"rm {files['in_files'][idx]} {files['trim_files'][idx]}")
    else:
        '''If not trimming, just rename the filtered files to the trim output fnames'''
        for idx, fi in enumerate(files['in_files']):
            shell(f"mv {fi} {files['clean_files'][idx]}")
    end_sec_print("INFO: Read Trimming complete.")
