import pytest
import os
import mock

from app.utils.shell_cmds import shell
from app.src.preprocess import run_kraken, rm_existing_kraken
from test.utils import get_default_args


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
