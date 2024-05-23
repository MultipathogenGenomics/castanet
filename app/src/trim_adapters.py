import os
from app.utils.shell_cmds import shell, stoperr
from app.utils.utility_fns import read_fa, enumerate_read_files
from app.utils.system_messages import end_sec_print


def run_trim(p, trim_path='java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar'):
    '''Call Trimmomatic trim CLI tool; check it worked; remove interim files'''
    p["ExpDir"] = f"{p['ExpDir']}/"
    CLEAN_UP = True
    files = {
        "in_files": [f"{p['ExpRoot']}/{p['ExpName']}/{p['ExpName']}_1_filt.fastq", f"{p['ExpRoot']}/{p['ExpName']}/{p['ExpName']}_2_filt.fastq"],
        "clean_files": [f"{p['ExpRoot']}/{p['ExpName']}/{p['ExpName']}_1_clean.fastq", f"{p['ExpRoot']}/{p['ExpName']}/{p['ExpName']}_2_clean.fastq"],
        "trim_files": [f"{p['ExpRoot']}/{p['ExpName']}/{p['ExpName']}_1_trimmings.fq", f"{p['ExpRoot']}/{p['ExpName']}/{p['ExpName']}_2_trimmings.fq"]
    }
    if not os.path.exists(files['in_files'][0]):
        '''If default filtered files not available, find input fastqs'''
        files['in_files'] = enumerate_read_files(p["ExpDir"])
        CLEAN_UP = False

    if p["DoTrimming"]:
        end_sec_print("INFO: Read Trimming beginning.")
        shell(f"{trim_path} PE -threads {p['NThreads']} {files['in_files'][0]} {files['in_files'][1]} {files['clean_files'][0]} {files['trim_files'][0]} {files['clean_files'][1]} {files['trim_files'][1]} ILLUMINACLIP:{p['AdaptP']}:2:10:7:1:true MINLEN:{p['TrimMinLen']}")

        '''Test for success'''
        if not os.path.exists(files['clean_files'][0]) or not os.path.exists(files['clean_files'][1]):
            stoperr(f"Trimming produced empty files. Check your TrimMinLen parameter is not too short for your sequences and that Trimmomatic is isntalled (you may use the dependency_check endpoint to check your installation).")
        if os.stat(files['clean_files'][0]).st_size < 2 or os.stat(files['clean_files'][1]).st_size < 2:
            stoperr(f"Trimming produced empty files. Check your TrimMinLen parameter is not too short for your sequences and that Trimmomatic is isntalled (you may use the dependency_check endpoint to check your installation).")
        if CLEAN_UP:
            for idx, fi in enumerate(files['in_files']):
                shell(
                    f"rm {files['in_files'][idx]} {files['trim_files'][idx]}")
    else:
        '''If not trimming, just rename the filtered files to the trim output fnames'''
        for idx, fi in enumerate(files['in_files']):
            shell(f"mv {fi} {files['clean_files'][idx]}")
    end_sec_print("INFO: Read Trimming complete.")
