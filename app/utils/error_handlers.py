import os
import pandas as pd
import subprocess as sp
from collections import deque

from app.utils.shell_cmds import loginfo, stoperr, read_line


def error_handler_filter_keep_reads(argies):
    '''
    Print errors, and exit if necessary, on bad input data.
    Extend incl/exlc ID lists from input names.
    '''
    retain_list, exclude_list, argies["o"] = deque(), deque(), []
    '''Check in path, set up output path'''
    for inpath in argies["input_file"]:
        if not os.path.isfile(inpath):
            stoperr(f'Unable to open FastQ file {inpath}.')
        outstem = os.path.basename(inpath.split('.gz')[0]) if inpath.endswith(
            '.gz') else os.path.basename(inpath)
        '''Append suffix "filt" to output file'''
        outpath = f'experiments/{argies["ExpName"]}/{os.path.splitext(outstem)[0]}_filt.fastq'
        argies["o"].append(outpath)
    loginfo(f'Output files {argies["o"]}')

    '''Check input files'''
    if len(argies["input_file"]) != len(argies["o"]):
        stoperr('Could not create output paths for all given input files.')
    if not os.path.isfile(argies["kraken"]):
        stoperr(
            f'Unable to open Kraken file {argies["kraken"]} for input {argies["input_file"]}.')

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
        stoperr('Nothing to do. Exiting.')

    return argies["o"], argies["ExcludeIds"], argies["RetainIds"]


def error_handler_parse_bam_positions(argvs):
    if len(argvs) < 2 or '-h' in argvs:
        stoperr(f'Usage: samtools view MyBamFile | {argvs[0]} \n\n')


def error_handler_analysis(argies):
    '''Print errors, and exit if necessary, on bad input data. Make outdir'''
    '''Validate main infput file (dataframe from processed BAM files for this pool)'''
    if not os.path.isfile(argies["input_file"]):
        stoperr('Unable to open input file {0}.'.format(argies["input_file"]))

    '''Validate sample info file'''
    # RM < TODO Samples input disabled.
    # if argies["Samples"] and not os.path.isfile(argies["Samples"]):
    #     stoperr(
    #         f'Unable to open sample data from input file {argies["Samples"]}.')
    # with open(argies["Samples"]) as samples_inf:
    #     samples_header_check = samples_inf.readline()
    #     if ('sampleid' not in samples_header_check) or ('pt' not in samples_header_check) or ('rawreadnum' not in samples_header_check):
    #         stoperr(
    #             f'{argies["Samples"]} must contain at least the following columns: pt, sampleid, rawreadnum')
    if argies["Clin"]:
        if not os.path.isfile(argies["Clin"]):
            stoperr(
                f'Unable to open clinical data from input file {argies["Clin"]}.')
        '''Validate clinical info file'''
        with open(argies["Clin"]) as clin_inf:
            clin_header_check = clin_inf.readline().split(',')
            if ('pt' not in clin_header_check):
                stoperr(
                    f'{argies["Clin"]} must contain at least the following columns: pt, clin_int')

    '''Validate batch name'''
    if not argies["ExpName"].replace('_', '').isalnum():
        stoperr(
            f'{argies["ExpName"]} is not a valid batch name. Must be alphanumeric.')

    '''Open data frame'''
    try:
        df = pd.read_csv(argies["input_file"], compression=('gzip' if argies["input_file"].endswith('.gz') else None), header=None,
                         names=['n', 'target_id', 'startpos', 'maplen', 'sampleid'])
    except (IOError, TypeError, pd.errors.ParserError) as e:
        stoperr(f'Failed to read dataframe from {df} : {e}')
    if argies["DepthInf"] and not os.path.isfile(argies["DepthInf"]):
        stoperr(f'Unable to open precomputed depth file {argies["DepthInf"]}.')
    return df


def error_handler_consensus_ref_corrected(a, tar_name):
    '''Don't construct a ref corrected genome if conditions met'''
    if a["GtOrg"] == "" and a["GtFile"] == "":
        print("WARNING: Not calling reference corrected consensus as no evaluation arguments were specified (GtOrg, GtFile)")
        return True
    elif bool(a["GtOrg"] == "") ^ bool(a["GtFile"] == ""):
        print("WARNING: Not calling reference corrected consensus, both a GtOrg and GtFile need to be specified")
        return True
    if tar_name != a['GtOrg']:
        print(f"Target {tar_name} is not the GT organism")
        return True
    return False
