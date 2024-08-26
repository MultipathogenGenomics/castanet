import os
from app.api import run_end_to_end, process_payload, do_batch
from app.utils.system_messages import end_sec_print
from app.utils.error_handlers import error_handler_api, stoperr
from app.utils.argparsers import parse_arguments_lite


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
    parser, bool_fields = parse_arguments_lite()
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
