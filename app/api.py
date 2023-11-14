import os
import re
import time
from fastapi import FastAPI, BackgroundTasks
from fastapi.encoders import jsonable_encoder

from app.utils.shell_cmds import stoperr, shell
from app.utils.timer import timing
from app.utils.system_messages import banner, end_sec_print
from app.utils.utility_fns import make_exp_dir
from app.utils.write_logs import write_input_params
from app.utils.eval import Evaluate
from app.utils.generate_probe_files import ProbeFileGen
from app.utils.combine_batch_output import combine_output_csvs
from app.src.preprocess import run_kraken
from app.src.filter_keep_reads import FilterKeepReads
from app.src.trim_adapters import run_trim
from app.src.map_reads_to_ref import run_map
from app.src.generate_counts import run_counts
from app.src.consensus import Consensus
from app.src.analysis import Analysis
from app.src.post_filter import run_post_filter
from app.utils.api_classes import (Batch_eval_data, E2e_eval_data, E2e_data, Preprocess_data, Filter_keep_reads_data,
                                   Trim_data, Mapping_data, Count_map_data, Analysis_data,
                                   Post_filter_data, Consensus_data, Eval_data, Convert_probe_data)

description = """
CASTANET is software for analysis of targeted metagenomics sequencing data, originally by tgolubch (https://github.com/tgolubch) and refactored to Python3 by mayne941 (https://github.com/Mayne941).
"""

tags_metadata = [
    {
        "name": "End to end pipeline",
        "description": "Run an end-to-end Castanet job",
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
    payload["NThreads"] = os.cpu_count()
    if "ConsensusMinD" in payload.keys():
        if payload["ConsensusMinD"] <= 2:
            stoperr(f"Consuensus min depth must exceed 2, otherwise you would inherit sections of reference sequence in the final remapped consensus.")
    write_input_params(payload)
    return payload


'''Dev Endpoints'''


@app.get("/", tags=["Dev endpoints"])
async def read_root() -> dict:
    return {"response": "API is healthy. Append the current URL to include '/docs/' at the end to visit the GUI."}


@app.post("/batch_eval/", tags=["Dev endpoints"])
async def batch(payload: Batch_eval_data) -> str:
    st = time.time()
    payload = process_payload(payload)
    payload["StartTime"] = time.time()
    SeqNames = get_batch_seqnames(payload["BatchName"])
    agg_analysis_csvs, agg_analysis_name = [], f'{payload["ExpName"]}.csv'
    errs = []
    for i in SeqNames:
        try:
            payload["ExpDir"] = "/".join(i[1][1].split("/")[:-1])
            payload["ExpName"] = payload["SeqName"] = i[0]
            agg_analysis_csvs.append(
                f"experiments/{payload['ExpName']}/{payload['SeqName']}_depth_with_clin.csv")
            run_end_to_end(payload)
            do_eval(payload)
        except Exception as ex:
            errs.append(i[0])
            end_sec_print(
                f"REGISTERED ERROR {payload['SeqName']} WITH EXCEPTION: {ex}")
    msg = combine_output_csvs(agg_analysis_csvs, agg_analysis_name)
    print(
        f"***\nBatch complete. Time to complete: {time.time() - st} ({(time.time() - st)/len(SeqNames)} per sample)\n{msg}\nFailed to process following samples: {errs}***")
    return "Task complete. See terminal output for details."


def get_batch_seqnames(batch_name) -> list:
    fstems = []
    folders = sorted(os.listdir(batch_name))
    for folder in folders:
        f_full = [f'{batch_name}/{folder}/{"_".join(i.split("_")[:-1])}' for i in sorted(
            os.listdir(f"{batch_name}/{folder}")) if re.match(r"[\s\S]*?\.fastq.gz", i)]
        if "_R1" in f_full[0] or "_R2" in f_full[0]:
            temp = []
            raw_names = [f'{batch_name}/{folder}/{"_".join(i.split("_"))}' for i in sorted(
                os.listdir(f"{batch_name}/{folder}")) if re.match(r"[\s\S]*?\.fastq.gz", i)]
            for i in range(0, 2):
                temp.append(
                    f'{raw_names[i].replace(f"_R1", "").replace("_R2","").split(".fastq.gz")[0]}_{i+1}.fastq.gz')
                shell(f"mv {raw_names[i]} {temp[i]}")
        else:
            temp = [i for i in os.listdir(f"{batch_name}/{folder}")]
        f = ["_".join(i.split("_")[:-1]).split("/")[-1]
             for i in temp if re.match(r"[\s\S]*?\.fastq.gz", i)]
        assert len(
            f) == 2,  "Incorrect number of files in directory, please ensure your experiment folder contains only two fasta.gz files."
        assert f[0] == f[1], "Inconsistent naming between paired read files, please revise your naming conventions."
        fstems.append([list(set(f))[0], f_full])
    return fstems


@app.post("/end_to_end/", tags=["Dev endpoints"])
async def end_to_end(payload: E2e_data) -> None:
    payload = process_payload(payload)
    end_sec_print(
        f"INFO: Starting run, saving results to {payload['ExpName']}.")
    run_end_to_end(payload)


@app.post("/end_to_end_eval/", tags=["Dev endpoints"])
async def end_to_end_eval(payload: E2e_eval_data) -> None:
    payload = process_payload(payload)
    payload["StartTime"] = time.time()
    end_sec_print(
        f"INFO: Starting run, saving results to {payload['ExpName']}.")
    run_end_to_end(payload)
    do_eval(payload)


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


@app.post("/end_to_end/", tags=["End to end pipeline"])
async def end_to_end(payload: E2e_data) -> None:
    payload = process_payload(payload)
    end_sec_print(
        f"INFO: Starting run, saving results to {payload['ExpName']}.")
    run_end_to_end(payload)


@timing
def run_end_to_end(payload) -> str:
    end_sec_print(f"INFO: Starting run, experiment: {payload['ExpName']}")
    make_exp_dir(payload["ExpName"])
    run_kraken(payload)
    do_filter_keep_reads(payload)
    run_trim(payload)
    run_map(payload)
    run_counts(payload)
    run_analysis(payload)
    if payload["PostFilt"]:
        run_post_filter(payload)
    do_consensus(payload)
    return "Task complete. See terminal output for details."


@app.post("/preprocess/", tags=["Individual pipeline functions"])
async def preprocess(payload: Preprocess_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(payload["ExpName"])
    run_kraken(payload)
    return "Task complete. See terminal output for details."


@app.post("/filter_keep_reads/", tags=["Individual pipeline functions"])
async def filter_keep_reads(payload: Filter_keep_reads_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(payload["ExpName"])
    do_filter_keep_reads(payload)
    return "Task complete. See terminal output for details."


def do_filter_keep_reads(payload) -> None:
    cls = FilterKeepReads(payload)
    cls.main()


@app.post("/trim_data/", tags=["Individual pipeline functions"])
async def trim_data(payload: Trim_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(payload["ExpName"])
    run_trim(payload)
    return "Task complete. See terminal output for details."


@app.post("/mapping/", tags=["Individual pipeline functions"])
async def mapping(payload: Mapping_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(payload["ExpName"])
    run_map(payload)
    return "Task complete. See terminal output for details."


@app.post("/generate_counts/", tags=["Individual pipeline functions"])
async def count_mapped(payload: Count_map_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(payload["ExpName"])
    run_counts(payload)
    return "Task complete. See terminal output for details."


@app.post("/consensus/", tags=["Individual pipeline functions"])
async def consensus(payload: Consensus_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(payload["ExpName"])
    do_consensus(payload)
    return "Task complete. See terminal output for details."


def do_consensus(payload):
    clf = Consensus(payload)
    clf.main()


@app.post("/analysis/", tags=["Individual pipeline functions"])
async def analysis(payload: Analysis_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(payload["ExpName"])
    run_analysis(payload)
    return "Task complete. See terminal output for details."


def run_analysis(payload) -> None:
    cls = Analysis(payload)
    cls.main()


@app.post("/post_filter/", tags=["Individual pipeline functions"])
async def post_filter(payload: Post_filter_data) -> str:
    payload = process_payload(payload)
    make_exp_dir(payload["ExpName"])
    run_post_filter(payload)
    return "Task complete. See terminal output for details."


'''Convenience functions'''


@app.post("/convert_probes/", tags=["Convenience functions"])
async def convertprobes(payload: Convert_probe_data) -> str:
    payload = jsonable_encoder(payload)
    clf = ProbeFileGen(payload)
    clf.main()
    return f"Task complete. Output saved to: {payload['OutFolder']}/{payload['OutFileName']}.fasta / .csv."
