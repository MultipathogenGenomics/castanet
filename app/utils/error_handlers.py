import os
import pandas as pd
import subprocess as sp
from collections import deque
from termcolor import colored

from app.utils.shell_cmds import loginfo, stoperr, read_line


def error_handler_filter_keep_reads(argies):
    '''
    Print errors, and exit if necessary, on bad input data.
    Extend incl/exlc ID lists from input names.
    '''
    retain_list, exclude_list, argies["o"] = deque(), deque(), []
    '''Check in path, set up output path'''
    cnt = 1
    for inpath in argies["input_file"]:
        if not os.path.isfile(inpath):
            stoperr(f'Unable to open FastQ file {inpath}.')
        # outstem = os.path.basename(inpath.split('.gz')[0]) if inpath.endswith(
        #     '.gz') else os.path.basename(inpath) # TODO DEPRECATED
        # outpath = f'experiments/{argies["ExpName"]}/{os.path.splitext(outstem)[0]}_filt.fastq'
        '''Append suffix "filt" to output file'''
        argies["o"].append(
            f'{argies["SaveDir"]}/{argies["ExpName"]}/{argies["ExpName"]}_{cnt}_filt.fastq')
        cnt += 1
    loginfo(f'Output files {argies["o"]}')

    '''Check input files'''
    if len(argies["input_file"]) != len(argies["o"]):
        stoperr('Could not create output paths for all given input files.')
    if not os.path.isfile(argies["kraken"]):
        stoperr(
            f'Unable to open Kraken file {argies["kraken"]} for input {argies["input_file"]}. Have you run the preprocess command before this one to create the kraken file?')

    '''Check lineage file; iterate over linF, get ret/excl IDs from names'''
    if argies["RetainNames"] or argies["ExcludeNames"]:
        if not os.path.isfile(argies["LineageFile"]):
            stoperr(f'Unable to open lineage file {argies["LineageFile"]}.')
        else:
            taxa = dict.fromkeys(argies["RetainNames"].split(','), retain_list)
            '''If same taxid in both lists, exclusion takes precedence over retention'''
            taxa.update(dict.fromkeys(
                argies["ExcludeNames"].split(','), exclude_list))
            if '' in taxa:
                del taxa['']
            cmd = 'zcat' if argies["LineageFile"].endswith('.gz') else 'cat'
            handle = sp.Popen(
                (cmd, argies["LineageFile"]), bufsize=8192, stdout=sp.PIPE).stdout
            line = read_line(handle)
            while line:
                for taxname, taxlist in taxa.items():
                    if f',{taxname},' in line:
                        taxlist.append(line.split(',', 1)[0])
                line = read_line(handle)

    '''Check NCBI TaxIDs to retain/exclude'''
    if not type(argies["ExcludeIds"]) == frozenset:
        try:
            argies["ExcludeIds"] = (frozenset(argies["ExcludeIds"].split(
                ',')) | frozenset(exclude_list)) - frozenset([''])
        except:
            stoperr(f'TaxID(s) {argies["ExcludeIds"]} invalid.')
        try:
            argies["retain"] = (frozenset(argies["RetainIds"].split(
                ',')) | frozenset(retain_list)) - frozenset([''])
        except:
            stoperr(f'TaxID(s) {argies["RetainIds"]} invalid.')
    if not argies["ExcludeIds"] and not argies["RetainIds"]:
        stoperr('User opted to use Kraken filtering but no parameters were provided to exclude reads: please check your ExcludeIds and RetainIds arguments or set DoKrakenPreprocess to false.')

    return argies["o"], argies["ExcludeIds"], argies["RetainIds"]


def error_handler_parse_bam_positions(argvs) -> None:
    if len(argvs) < 2 or '-h' in argvs:
        stoperr(f'Usage: samtools view MyBamFile | {argvs[0]} \n\n')


def error_handler_analysis(argies) -> pd.DataFrame:
    '''Print errors, and exit if necessary, on bad input data. Make outdir'''
    '''Validate main infput file (dataframe from processed BAM files for this pool)'''
    if not os.path.isfile(argies["input_file"]):
        stoperr('Unable to open input file {0}.'.format(argies["input_file"]))

    '''Validate sample info file'''
    if argies["Clin"] != "":
        if not os.path.isfile(argies["Clin"]):
            stoperr(
                f'Unable to open clinical data from input file {argies["Clin"]}.')
        '''Validate clinical info file'''
        with open(argies["Clin"]) as clin_inf:
            clin_header_check = clin_inf.readline().split(',')
            if ('pt' not in clin_header_check):
                stoperr(
                    f'{argies["Clin"]} must contain at least the following columns: pt, clin_int')

    '''Open data frame'''
    try:
        df = pd.read_csv(argies["input_file"], compression=('gzip' if argies["input_file"].endswith('.gz') else None), header=None,
                         names=['n', 'target_id', 'startpos',
                                'maplen', 'sampleid'],
                         on_bad_lines="skip")
    except (IOError, TypeError, pd.errors.ParserError) as e:
        stoperr(f'Failed to read dataframe from {df} : {e}')

    if df.empty:
        stoperr(f"Your Positions Count file is empty, meaning that Castanet didn't detect any significant hits in your input sample. This can sometimes mask an upstream problem, but may also mean that your sample is low quality and/or genuinely has nothing that maps to your mapping reference.")

    if argies["DepthInf"] and not os.path.isfile(argies["DepthInf"]):
        stoperr(f'Unable to open precomputed depth file {argies["DepthInf"]}.')
    return df


def error_handler_consensus_ref_corrected(a, tar_name) -> bool:
    '''Don't construct a ref corrected genome if conditions met'''
    if not "GtOrg" in a.keys() and not "GtFile" in a.keys():
        print("WARNING: Not calling reference corrected consensus as no evaluation arguments were specified (GtOrg, GtFile)")
        return True
    if a["GtOrg"] == "" and a["GtFile"] == "":
        print("WARNING: Not calling reference corrected consensus as no evaluation arguments were specified (GtOrg, GtFile)")
        return True
    elif bool(a["GtOrg"] == "") ^ bool(a["GtFile"] == ""):
        print("WARNING: Not calling reference corrected consensus, both a GtOrg and GtFile need to be specified")
        return True
    if tar_name != a['GtOrg']:
        print(
            f"INFO: Target {tar_name} is not the GT organism, so ref-adjusted consensus not being built for it")
        return True
    return False


def error_handler_api(ex):
    import traceback
    import logging
    import re
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    err = traceback.format_exc()
    if type(ex) == SystemError:
        err_short = ansi_escape.sub("", str(ex))
    else:
        print(colored(f"Unclassified Castanet error:", 'red'))
        err_short = "Unclassified Castanet error: " + err.split('\n')[-2]
    logging.error(err)
    return f"Castanet run failed, please see error message and terminal for more details: {err_short}"


def check_readf_ext(fnames):
    '''DEPRECATED'''
    if "fq" in fnames[0]:
        ext = "fq"
    elif "fastq" in fnames[0]:
        ext = "fastq"
    else:
        stoperr(f"Input files ({fnames}) don't appear to be of type .fq or .fastq (regardless of gz compression). Please check your files are in the correct format and that there's nothing else in the folder.")
    if "gz" in fnames[0]:
        ext += ".gz"
    return ext


def error_handler_cli(out, out_fname, tool, test_out_f=True, test_f_size=False):
    cli_specific_errs = get_cli_tool_errors(tool)
    default_guidance = ""
    if "guidance" in cli_specific_errs.keys():
        default_guidance = cli_specific_errs["guidance"]
    if "command not found" in out.lower() or "segmentation fault" in out.lower() or not cli_specific_errs["healthy_msg"] in out.lower():
        if out == "" and test_f_size:
            # Don't fail here as it represents casting to an out file, which would have no response
            ...
        else:
            stoperr(f"{tool} doesn't seem to be installed or threw an error not recognised by the Castanet test suite. Please check the Castanet readme for installation instructions. {default_guidance}")
    if test_out_f:
        if not os.path.exists(out_fname):
            stoperr(
                f"{tool} seems to be installed but didn't produce output. Please check your installation manually.")
    if test_f_size:
        if os.stat(out_fname).st_size < 2:
            stoperr(f"{tool} produced an empty output file."
                    f"{cli_specific_errs['guidance']}")


def get_cli_tool_errors(cli_tool):
    err_objs = {
        "kraken": {"healthy_msg": "loading database information",
                   "guidance": "Ensure it's installed and that you've got a Kraken2 database in the place Castanet expects it (via KrakenDbDir argument)."},
        "java": {"healthy_msg": "usage: java"},
        "trimmomatic": {"healthy_msg": "trimmomaticpe: started with arguments",
                        "guidance": "Trimming produced empty files. Check your TrimMinLen parameter is not too short for your sequences and that Trimmomatic is isntalled (you may use the dependency_check endpoint to check your installation)."},
        "bwa-mem2": {"healthy_msg": "looking to launch executable",
                     "guidance": f"BWA-MEM2 produced an empty BAM file. Check your BWA-MEM2 installation and that your input reads are of sufficent quality. This might also indicate an out of memory error, if you're crunching a huge dataset: if so, rerun the experiment with less cores (NThreads parameter)."},
        "samtools": {"healthy_msg": "program: samtools"},
        "mafft": {"healthy_msg": "generating a scoring matrix for nucleotide"},
        "viral_consensus": {"healthy_msg": "",
                            "guidance": "ViralConsensus did not function properly when called. Check it's installed"},
    }
    return err_objs[cli_tool]
