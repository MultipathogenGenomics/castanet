import argparse
import os
from app.api import run_end_to_end, process_payload, do_batch
from app.utils.system_messages import end_sec_print
from app.utils.error_handlers import error_handler_api, stoperr


def defaults():
    return {
        "ExpDir": "./data/eval/",
        "ExpName": "CastanetTest",
        "SaveDir": "./experiments",
        "RefStem": "data/eval/ref.fa",
        "SingleEndedReads": False,
        "MatchLength": 40,
        "DoTrimming": True,
        "TrimMinLen": 36,
        "DoKrakenPrefilter": True,
        "LineageFile": "data/ncbi_lineages_2023-06-15.csv.gz",
        "ExcludeIds": "9606",
        "RetainIds": "",
        "RetainNames": "",
        "ExcludeNames": "Homo",
        "ConsensusMinD": 10,
        "ConsensusCoverage": 30,
        "ConsensusMapQ": 1,
        "ConsensusCleanFiles": True,
        "GtFile": "",
        "GtOrg": "",
        "KrakenDbDir": "kraken2_human_db/",
        "KeepDups": True,
        "Clin": "",
        "DepthInf": "",
        "SamplesFile": "",
        "PostFilt": False,
        "AdaptP": "data/all_adapters.fa",
        "NThreads": "auto"
    }


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Castanet Lite (Beta)"
    )
    '''N.b. Argparse is SO UNBELIEVABLY FUCKING SHIT that it can't evaluate booleans on optional arguments with a default, so we need to eval() strings passed to it later'''
    parser.add_argument('-Batch', required=False, default=False, type=str,
                        help="If True, conduct a batch run analysing multiple datasets; expects your ExpDir folder to contain sub-folders, each containing two (paired) read files.")
    parser.add_argument('-BAM', required=False, default=False, type=str,
                        help="If True, launch Castanet in BAM process mode, in which the software will skip initial processing and look in your input folder/s for BAM files rather than .fastq.gz")
    parser.add_argument('-ExpDir', required=True, type=str,
                        help="Folder containing your two paired read files, OR if batch = True, folder containing your experiment sub folders.")
    parser.add_argument('-ExpName', required=True, type=str,
                        help="Name for your experiment data.")
    parser.add_argument('-SaveDir', required=True, type=str,
                        help="Folder to save your experiment data.")
    parser.add_argument('-RefStem', required=True, type=str,
                        help="File location for mapping reference (fasta).")
    parser.add_argument('-DoKrakenPrefilter', required=False, default=True, type=str,
                        help="If True, do an initial Kraken pre-filter, using database specifid in -KrakenDbDir and list of NCBI TaxID(s) to exclude from -ExcludeIds.")
    parser.add_argument('-KrakenDbDir', required=False, default="kraken2_human_db", type=str,
                        help="If -DoKrakenPrefilter = True, path to Kraken2 database to do filtering.")
    parser.add_argument('-ExcludeIds', required=False, default="9606", type=str,
                        help="If -DoKrakenPrefilter = True, filter this/these NBCI TaxIDs (list of integers, separated by commas with no spaces)")
    parser.add_argument('-DoTrimming', required=False, default=True, type=str,
                        help="If True, use Trimmomatic to remove adapters and low quality sequences.")
    parser.add_argument('-DoConsensus', required=False, default=True, type=str,
                        help="If True, run the Castanet consensus sequence pipeline stage.")
    bool_fields = ["Batch", "BAM", "DoKrakenPrefilter",
                   "DoTrimming", "DoConsensus"]
    return parser, bool_fields


def populate_request(payload, parser, bool_vals):
    '''Overwrite default request object with argparser args'''
    for key, val in parser.__dict__.items():
        if key in bool_vals and type(val) == str:
            '''If argparse argument is string and should be bool. Seriously, fuck argument parsers, let's get out of 1994 people.'''
            val = eval(val)
        payload[key] = val
    return payload


def tests(payload):
    '''ExpDir'''
    if not os.path.exists(payload['ExpDir']):
        os.mkdir(payload['ExpDir'])
    '''ExpName'''
    if not payload['ExpName'].isprintable() or " " in payload['ExpName']:
        stoperr(f"Your experiment name (ExpName) isn't valid because it contains unprintable characters and/or spaces.")
    '''SaveDir'''
    if not os.path.exists(payload['SaveDir']):
        os.mkdir(payload['SaveDir'])
    '''RefStem'''
    if not os.path.isfile(payload['RefStem']):
        stoperr(f"Your mapping reference file (RefStem) doesn't exist: check spelling.")


def main():
    parser, bool_fields = parse_arguments()
    payload = populate_request(defaults(), parser.parse_args(), bool_fields)
    tests(payload)
    end_sec_print(f"Calling Castanet Lite with following arguments: {payload}")
    try:
        payload = process_payload(payload)
        end_sec_print(
            f"INFO: Starting run, saving results to {payload['ExpName']}.")
        if payload["BAM"]:
            '''Do analyse_my_bam pipeline'''
            if not payload["Batch"]:
                msg = run_end_to_end(payload, start_with_bam=True)
            else:
                payload["DataFolder"] = payload["ExpDir"]
                del payload["ExpDir"]
                msg = do_batch(payload, start_with_bam=True)
        else:
            '''Do end_to_end pipeline'''
            if not payload["Batch"]:
                '''Single end to end'''
                msg = run_end_to_end(payload)
            else:
                payload["DataFolder"] = payload["ExpDir"]
                del payload["ExpDir"]
                msg = do_batch(payload)
            return msg
    except Exception as ex:
        return error_handler_api(ex)


if __name__ == "__main__":
    main()
