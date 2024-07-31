import os
from app.utils.shell_cmds import shell, logerr
from app.utils.utility_fns import enumerate_read_files
from app.utils.system_messages import end_sec_print
from app.utils.error_handlers import error_handler_cli


# def run_trim(p, trim_path='java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar'): # RM < TODO test after switching to bioconda install
def run_trim(p, trim_path='trimmomatic'):
    '''Call Trimmomatic trim CLI tool; check it worked; remove interim files'''
    p["ExpDir"] = f"{p['ExpDir']}/"
    CLEAN_UP = True
    files = {
        "in_files": [f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_1_filt.fastq", f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_2_filt.fastq"],
        "clean_files": [f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_1_clean.fastq", f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_2_clean.fastq"],
        "trim_files": [f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_1_trimmings.fq", f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_2_trimmings.fq"]
    }
    if not os.path.exists(files['in_files'][0]):
        '''If default filtered files not available, find input fastqs'''
        files['in_files'] = enumerate_read_files(p["ExpDir"])
        CLEAN_UP = False

    if p["DoTrimming"]:
        end_sec_print("INFO: Read Trimming beginning.")
        out = shell(f"{trim_path} PE -threads {p['NThreads']} {files['in_files'][0]} {files['in_files'][1]} {files['clean_files'][0]} {files['trim_files'][0]} {files['clean_files'][1]} {files['trim_files'][1]} ILLUMINACLIP:{p['AdaptP']}:2:10:7:1:true MINLEN:{p['TrimMinLen']}", is_test=True)
        error_handler_cli(out, files['clean_files']
                          [0], "trimmomatic", test_f_size=True)

        if CLEAN_UP:
            for idx, fi in enumerate(files['in_files']):
                shell(
                    f"rm {files['in_files'][idx]} {files['trim_files'][idx]}")
    else:
        '''If not trimming, just rename the filtered files to the trim output fnames'''
        logerr(f"Skipping trimming as you specified to")
        for idx, fi in enumerate(files['in_files']):
            shell(f"mv {fi} {files['clean_files'][idx]}")
    end_sec_print("INFO: Read Trimming complete.")
