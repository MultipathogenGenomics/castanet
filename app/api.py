import os
import re
import time
from fastapi import FastAPI
from fastapi.encoders import jsonable_encoder

from app.utils.shell_cmds import stoperr
from app.utils.timer import timing
from app.utils.system_messages import banner, end_sec_print
from app.utils.utility_fns import make_exp_dir, enumerate_read_files, read_fa
from app.utils.write_logs import write_input_params
from app.utils.eval import Evaluate
from app.utils.error_handlers import error_handler_api
from app.utils.generate_probe_files import ProbeFileGen
from app.utils.combine_batch_output import combine_output_csvs
from app.utils.dependency_check import Dependencies
from app.src.preprocess import run_kraken
from app.src.filter_keep_reads import FilterKeepReads
from app.src.trim_adapters import run_trim
from app.src.map_reads_to_ref import run_map
from app.src.generate_counts import run_counts
from app.src.consensus import Consensus
from app.src.analysis import Analysis
from app.src.amplicons import Amplicons
from app.src.post_filter import run_post_filter
from app.utils.test_imports import import_test
from app.utils.api_classes import (Batch_eval_data, E2e_eval_data, E2e_data, Preprocess_data, Filter_keep_reads_data, Amp_e2e_data,
                                   Trim_data, Mapping_data, Count_map_data, Analysis_data, Dep_check_data, Amplicon_data,
                                   Post_filter_data, Consensus_data, Eval_data, Convert_probe_data, Bam_workflow_data)

import_test()

description = """
CASTANET is software for analysis of targeted & metagenomic sequencing data, originally by tgolubch (https://github.com/tgolubch) and refactored to Python3 by mayne941 (https://github.com/Mayne941).
"""

tags_metadata = [
    {
        "name": "End to end pipelines",
        "description": "Run end-to-end Castanet jobs",
    },
    {
        "name": "Individual pipeline functions",
        "description": "Run individual functions from the Castanet pipeline",
    },
    {
        "name": "Convenience functions",
        "description": "Supplementary functions for data processing",
    },
    {
        "name": "Dev endpoints",
        "description": "Developer tools"
    }
]

banner()

app = FastAPI(
    title="Castanet",
    version="0.3",
    description=description,
    contact={
        "name": "Nuffield Department of Medicine, University of Oxford",
        "url": "https://www.ndm.ox.ac.uk/",
    },
    license_info={
        "name": "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)",
        "url": "https://creativecommons.org/licenses/by-nc/4.0/",
    },
    openapi_tags=tags_metadata
)

'''Utility Functions'''


def process_payload(payload) -> dict:
    payload = jsonable_encoder(payload)

    if "NThreads" in payload.keys():
        if type(payload["NThreads"]) == str:
            if payload["NThreads"] == "auto":
                import psutil
                import sys
                import numpy
                n_cpus = os.cpu_count()
                vmem = psutil.virtual_memory()[1]
                refseqs = read_fa(payload["RefStem"])
                sizes = []
                for i in refseqs:
                    sizes.append(sys.getsizeof(i[1]) * 28)
                size_vs_vmem = sum(sizes) * n_cpus
                if size_vs_vmem >= vmem:
                    payload["NThreads"] = numpy.ceil(
                        vmem/(max(sizes) * os.cpu_count()))
                else:
                    payload["NThreads"] = n_cpus
                end_sec_print(
                    f"Setting AUTO NThreads to: {payload['NThreads']}")
            elif payload["NThreads"] == "hpc":
                payload["NThreads"] == 1
            else:
                stoperr(
                    f"NThreads parameter should either be an integer, or 'auto' or 'hpc'.")
        elif type(payload["NThreads"]) == int:
            pass
        else:
            stoperr(
                f"NThreads parameter should either be an integer, or 'auto' or 'hpc'.")

    if "ConsensusMinD" in payload.keys():
        if payload["ConsensusMinD"] <= 2:
            stoperr(f"Consuensus min depth must exceed 2, otherwise you would inherit sections of reference sequence in the final remapped consensus.")

    write_input_params(payload)
    return payload


'''Dev Endpoints'''


@app.get("/", tags=["Dev endpoints"])
async def read_root() -> dict:
    return {"response": "API is healthy. Append the current URL to include '/docs/' at the end to visit the GUI."}


@app.post("/check_dependencies/", tags=["Convenience functions"])
async def check_deps(payload: Dep_check_data) -> str:
    try:
        clf = Dependencies(jsonable_encoder(payload))
        return clf.main()
    except Exception as ex:
        return error_handler_api(ex)


@app.post("/batch/", tags=["End to end pipelines"])
async def batch(payload: Batch_eval_data) -> str:
    payload = process_payload(payload)
    msg = do_batch(payload)
    return msg


def do_batch(payload, start_with_bam=False):
    st = time.time()
    payload["BatchName"] = payload["DataFolder"]
    payload["StartTime"] = st
    agg_analysis_csvs, agg_analysis_name = [], f'{payload["ExpName"]}.csv'
    errs = []

    if not start_with_bam:
        '''Standard end to end pipelines'''
        SeqNamesList = [enumerate_read_files(
            folder, payload["BatchName"]) for folder in sorted(os.listdir(payload["BatchName"])) if not folder == "__pycache__"]
        SeqNamesList = [i for i in SeqNamesList if not i == []]
    else:
        '''BAM only pipelines'''
        SeqNamesList = [os.path.join(dp, f) for dp, dn, fn in os.walk(
            os.path.expanduser(payload["BatchName"])) for f in fn]

    if len(SeqNamesList) == 0:
        stoperr(f"No files could be detected to analyse. "
                f"Castanet expects batch runs to point towards either: a data folder containing sub-folders for each sample, which should contain only 2 read files each (end to end batch pipelines); "
                f"or a folder containing multiple bam files (analyse bam batch pipelines). "
                f"Please refer to Castanet's readme for more details.")

    for SeqNames in SeqNamesList:
        try:
            if not start_with_bam:
                '''End to end pipelines'''
                exp_name = SeqNames[0].split(
                    "/")[-3]  # RM < TODO TEST THIS IS ROBUST WITH DIFFERENT FOL STRUCTURES
                payload["SeqNames"] = SeqNames
                payload["ExpDir"] = "/".join(SeqNames[0].split("/")[:-1])
                payload["ExpName"] = exp_name
            else:
                '''BAM only pipelines'''
                exp_name = SeqNames.split("/")[-2]
                payload["ExpDir"] = "/".join(SeqNames.split("/")[:-1])
                payload["ExpName"] = exp_name
            agg_analysis_csvs.append(
                f"experiments/{payload['ExpName']}/{exp_name}_depth_with_clin.csv")
            run_end_to_end(payload, start_with_bam)
            do_eval(payload)
        except Exception as ex:
            err = error_handler_api(ex)
            errs.append(exp_name)
            end_sec_print(
                f"REGISTERED ERROR {exp_name} WITH EXCEPTION: {err}")

    msg = combine_output_csvs(agg_analysis_csvs, agg_analysis_name)
    end_sec_print(msg)
    if len(errs) < 1:
        return "f***\nBatch complete. Time to complete: {time.time() - st} ({(time.time() - st)/len(SeqNames)} per sample)\n{msg}\nFailed to process following samples: {errs}***"
    else:
        return "Batch process task completed with errors. See terminal output for details."


@app.post("/end_to_end_eval/", tags=["Dev endpoints"])
async def end_to_end_eval(payload: E2e_eval_data) -> None:
    try:
        payload = process_payload(payload)
        payload["StartTime"] = time.time()
        end_sec_print(
            f"INFO: Starting run, saving results to {payload['ExpName']}.")
        msg = run_end_to_end(payload)
        do_eval(payload)
        return msg
    except Exception as ex:
        return f"Castanet run failed, please see error message and terminal for more details: {ex}"


@app.post("/evaulate/", tags=["Dev endpoints"])
async def evaluate(payload: Eval_data) -> None:
    payload = process_payload(payload)
    payload["StartTime"] = time.time()
    end_sec_print(
        f"INFO: Starting run evaluation.")
    do_eval(payload)


def do_eval(payload) -> None:
    clf = Evaluate(payload)
    clf.main()


'''Consumer endpoints'''


@app.post("/end_to_end/", tags=["End to end pipelines"])
async def end_to_end(payload: E2e_data) -> None:
    try:
        payload = process_payload(payload)
        end_sec_print(
            f"INFO: Starting run, saving results to {payload['ExpName']}.")
        msg = run_end_to_end(payload)
        return msg
    except Exception as ex:
        return error_handler_api(ex)


@app.post("/analyse_my_bam/", tags=["End to end pipelines"])
async def end_to_end(payload: Bam_workflow_data) -> None:
    try:
        payload = process_payload(payload)
        end_sec_print(
            f"INFO: Starting run, saving results to {payload['ExpName']}.")
        msg = run_end_to_end(payload, start_with_bam=True)
        return msg
    except Exception as ex:
        return error_handler_api(ex)


@app.post("/analyse_my_amplicons/", tags=["End to end pipelines"])
async def end_to_end(payload: Amp_e2e_data) -> None:
    try:
        payload = process_payload(payload)
        end_sec_print(
            f"INFO: Starting run, saving results to {payload['ExpName']}.")
        msg = run_amp_end_to_end(payload, start_with_bam=False)
        return msg
    except Exception as ex:
        return error_handler_api(ex)


@timing
def run_end_to_end(payload, start_with_bam=False) -> str:
    end_sec_print(f"INFO: Starting run, experiment: {payload['ExpName']}")
    make_exp_dir(f'{payload["SaveDir"]}/{payload["ExpName"]}')
    if not start_with_bam:
        if payload["DoKrakenPrefilter"]:
            run_kraken(payload)
        do_filter_keep_reads(payload)
        run_trim(payload)
        run_map(payload)
    run_counts(payload, start_with_bam)
    run_analysis(payload, start_with_bam)
    if payload["PostFilt"]:
        run_post_filter(payload)
    if payload["DoConsensus"]:
        do_consensus(payload, start_with_bam)
    return "Task complete. See terminal output for details."


@timing
def run_amp_end_to_end(payload, start_with_bam=False) -> str:
    end_sec_print(f"INFO: Starting run, experiment: {payload['ExpName']}")
    make_exp_dir(f'{payload["SaveDir"]}/{payload["ExpName"]}')
    if not start_with_bam:
        if payload["DoKrakenPrefilter"]:
            run_kraken(payload)
        do_filter_keep_reads(payload)
        run_trim(payload)
        run_map(payload)
    run_amplicons(payload)
    return "Task complete. See terminal output for details."


@app.post("/preprocess/", tags=["Individual pipeline functions"])
async def preprocess(payload: Preprocess_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(f'{payload["SaveDir"]}/{payload["ExpName"]}')
    run_kraken(payload)
    return "Task complete. See terminal output for details."


@app.post("/filter_keep_reads/", tags=["Individual pipeline functions"])
async def filter_keep_reads(payload: Filter_keep_reads_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(f'{payload["SaveDir"]}/{payload["ExpName"]}')
    do_filter_keep_reads(payload)
    return "Task complete. See terminal output for details."


def do_filter_keep_reads(payload) -> None:
    cls = FilterKeepReads(payload)
    cls.main()


@app.post("/trim_data/", tags=["Individual pipeline functions"])
async def trim_data(payload: Trim_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(f'{payload["SaveDir"]}/{payload["ExpName"]}')
    run_trim(payload)
    return "Task complete. See terminal output for details."


@app.post("/mapping/", tags=["Individual pipeline functions"])
async def mapping(payload: Mapping_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(f'{payload["SaveDir"]}/{payload["ExpName"]}')
    run_map(payload)
    return "Task complete. See terminal output for details."


@app.post("/generate_counts/", tags=["Individual pipeline functions"])
async def count_mapped(payload: Count_map_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(f'{payload["SaveDir"]}/{payload["ExpName"]}')
    run_counts(payload)
    return "Task complete. See terminal output for details."


@app.post("/analysis/", tags=["Individual pipeline functions"])
async def analysis(payload: Analysis_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(f'{payload["SaveDir"]}/{payload["ExpName"]}')
    run_analysis(payload)
    return "Task complete. See terminal output for details."


def run_analysis(payload, start_with_bam=False) -> None:
    cls = Analysis(payload, start_with_bam)
    cls.main()


@app.post("/consensus/", tags=["Individual pipeline functions"])
async def consensus(payload: Consensus_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(f'{payload["SaveDir"]}/{payload["ExpName"]}')
    do_consensus(payload)
    return "Task complete. See terminal output for details."


def do_consensus(payload, start_with_bam=False):
    clf = Consensus(payload, start_with_bam)
    clf.main()


@app.post("/post_filter/", tags=["Individual pipeline functions"])
async def post_filter(payload: Post_filter_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(f'{payload["SaveDir"]}/{payload["ExpName"]}')
    run_post_filter(payload)
    return "Task complete. See terminal output for details."


@app.post("/amplicons/", tags=["Individual pipeline functions"])
async def analysis(payload: Amplicon_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(f'{payload["SaveDir"]}/{payload["ExpName"]}')
    run_amplicons(payload)
    return "Task complete. See terminal output for details."


def run_amplicons(payload) -> None:
    cls = Amplicons(payload)
    cls.main()


'''Convenience functions'''


@app.post("/convert_mapping_reference/", tags=["Convenience functions"])
async def convertprobes(payload: Convert_probe_data) -> str:
    payload = jsonable_encoder(payload)
    clf = ProbeFileGen(payload)
    clf.main()
    return f"Task complete. Output saved to: {payload['OutFolder']}/{payload['OutFileName']}.fasta / .csv."
