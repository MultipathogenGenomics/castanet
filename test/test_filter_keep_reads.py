import pytest
import os
import shutil

from test.utils import get_random_str, make_rand_dir, create_test_file, get_default_args
from app.src.filter_keep_reads import FilterKeepReads
from app.src.preprocess import run_kraken
from app.utils.shell_cmds import shell
from app.utils.utility_fns import enumerate_read_files


def run_prefilter(prefilter):
    '''Initialise default args and directories'''
    p = get_default_args()
    p["DoKrakenPrefilter"], fstem = prefilter, f'{p["SaveDir"]}/{p["ExpName"]}'
    if os.path.exists(fstem):
        shutil.rmtree(fstem)
    os.mkdir(fstem)
    '''If doing prefilter, call Kraken'''
    if prefilter:
        run_kraken(p)
    '''Call filter keep reads'''
    clf = FilterKeepReads(p)
    clf.main()
    out_fs = enumerate_read_files(fstem)
    '''Check out files created with size > 2b'''
    for out_f in out_fs:
        assert os.stat(out_f).st_size > 2
    n_reads_orig = int(shell(
        f"samtools view -c {enumerate_read_files(p['ExpDir'])[0]}", ret_output=True).decode())
    n_reads_filt = int(shell(
        f'samtools view -c {p["SaveDir"]}/{p["ExpName"]}/{p["ExpName"]}_{1}_filt.fastq', ret_output=True).decode())
    '''Check outfile n reads is expected: n if no prefilter, n - 1 if prefilter'''
    if prefilter:
        expected = n_reads_orig - 1
    else:
        expected = n_reads_orig
    assert n_reads_filt == expected
    shutil.rmtree(fstem)


def test_filter_keep_reads_badinput():
    p = get_default_args()
    '''Test empty initialisation'''
    with pytest.raises(TypeError):
        _ = FilterKeepReads()
    '''Test breaks with empty/absent kraken file'''
    with pytest.raises(SystemError):
        _ = FilterKeepReads(p)


def test_filter_keep_reads_noprefilter():
    '''Test the `DoKrakenPrefilter` switch: should create filt files with same len as input'''
    run_prefilter(prefilter=False)


def test_filter_keep_reads_prefilter():
    '''Should remove 1 human synthetic read from dummy data'''
    run_prefilter(prefilter=True)


if __name__ == "__main__":
    test_filter_keep_reads_badinput()
    test_filter_keep_reads_noprefilter()
    test_filter_keep_reads_prefilter()
