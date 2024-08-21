import pytest
import os
import mock
from app.utils.shell_cmds import shell
from app.src.preprocess import run_kraken, rm_existing_kraken


def get_default_args():
    return {
        "ExpDir": "./data/",
        "ExpName": "test/experiments/test/",
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
        "DoConsensus": True,
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
        "NThreads": "auto",
        "SeqNames": ["data/eval/sim_reads_1.fastq.gz", "data/eval/sim_reads_2.fastq.gz"]
    }


@mock.patch(
    "app.src.preprocess.enumerate_read_files",
    mock.MagicMock(return_value=["fucked.fastq.gz"])
)
def test_run_kraken():
    p = get_default_args()
    with pytest.raises(SystemError):
        run_kraken(p)


def test_rm_existing_kraken():
    p = get_default_args()
    out_fnames = ["test/experiments/test/kraken.kraken",
                  "test/experiments/test/kraken_report.tsv"]
    shell(f"mkdir test/experiments/test/")
    shell(f"touch {out_fnames[0]} {out_fnames[1]}")
    rm_existing_kraken(p)
    assert not os.path.exists("test/experiments/test/kraken.kraken") or not os.path.exists(
        "test/experiments/test/kraken_report.tsv")


# if __name__ == "__main__":
#     test_rm_existing_kraken()
