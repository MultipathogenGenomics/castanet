import os
import re
from fastapi import FastAPI, BackgroundTasks
from fastapi.encoders import jsonable_encoder

from app.utils.timer import timing
from app.utils.system_messages import banner, end_sec_print
from app.utils.utility_fns import make_exp_dir
from app.utils.write_logs import write_input_params
from app.src.preprocess import run_kraken
from app.src.filter_keep_reads import FilterKeepReads
from app.src.trim_adapters import run_trim
from app.src.map_reads_to_ref import run_map
from app.src.generate_counts import run_counts
from app.src.consensus import Consensus
from app.src.analysis import Analysis
from app.src.post_filter import run_post_filter
from app.utils.api_classes import (Batch_data, E2e_data, Preprocess_data, Filter_keep_reads_data,
                                   Trim_data, Mapping_data, Count_map_data, Analysis_data,
                                   Post_filter_data, Consensus_data)

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
    license_info={  # RM < TODO CHECK LICENSE
        "name": "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)",
        "url": "https://creativecommons.org/licenses/by-nc/4.0/",
    },
    openapi_tags=tags_metadata
)

'''Utility Functions'''


def process_payload(payload) -> dict:
    payload = jsonable_encoder(payload)
    write_input_params(payload)
    return payload


'''Dev Endpoints'''


@app.get("/", tags=["Dev endpoints"])
async def read_root() -> dict:
    return {"response": "API is healthy. Append the current URL to include '/docs/' at the end to visit the GUI."}


'''Consumer endpoints'''


@app.post("/batch/", tags=["Dev endpoints"])
async def batch(payload: Batch_data) -> str:
    payload = process_payload(payload)
    SeqNames = get_batch_seqnames(payload["BatchName"])
    for i in SeqNames:
        payload["ExpDir"] = "/".join(i[1][1].split("/")[:-1])
        payload["ExpName"] = payload["SeqName"] = i[0]
        run_end_to_end(payload)
    return "Task complete. See terminal output for details."


def get_batch_seqnames(batch_name) -> list:
    fstems = []
    for folder in os.listdir(batch_name):
        f_full = [f'{batch_name}/{folder}/{"_".join(i.split("_")[:-1])}' for i in os.listdir(
            f"{batch_name}/{folder}") if re.match(r"[\s\S]*?\.fastq.gz", i)]
        f = ["_".join(i.split("_")[:-1]) for i in os.listdir(
            f"{batch_name}/{folder}") if re.match(r"[\s\S]*?\.fastq.gz", i)]
        assert len(
            f) == 2,  "Incorrect number of files in directory, please ensure your experiment folder contains only two fasta.gz files."
        assert f[0] == f[1], "Inconsistent naming between paired read files, please revise your naming conventions."
        fstems.append([list(set(f))[0], f_full])
    return fstems


@app.post("/end_to_end/", tags=["End to end pipeline"])
async def end_to_end(payload: E2e_data) -> None:
    payload = process_payload(payload)
    end_sec_print(
        f"INFO: Starting run, saving results to {payload['ExpName']}.")
    run_end_to_end(payload)


@timing
def run_end_to_end(payload) -> str:
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
