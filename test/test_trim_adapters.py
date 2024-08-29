import pytest
import os
import shutil

from test.utils import get_random_str, make_rand_dir, create_test_file, get_default_args
from app.src.trim_adapters import run_trim
from app.src.filter_keep_reads import FilterKeepReads
from app.src.preprocess import run_kraken


def init_kraken(p):
    run_kraken(p)
    clf = FilterKeepReads(p)
    clf.main()


def init_do_trim(do_trim, entry_via_pipelne):
    p = get_default_args()
    fstem = f"{p['SaveDir']}/{p['ExpName']}/"
    os.mkdir(fstem)
    if entry_via_pipelne:
        init_kraken(p)
    '''If not trimming, adjust payload'''
    outf_ext = "" if (do_trim or entry_via_pipelne) else ".gz"
    p["DoTrimming"] = do_trim
    '''Run fn to be tested'''
    run_trim(p)
    '''Check output has size > 2b'''
    trim_out = f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_1_clean.fastq{outf_ext}"
    assert os.stat(trim_out).st_size > 2
    shutil.rmtree(fstem)

# RM < TODO Test with dead input


def test_trim_adapters_trimentry_dotrim():
    '''When just trimming is called, relies on enumerating fresh files'''
    init_do_trim(do_trim=True, entry_via_pipelne=False)


def test_trim_adapters_trimentry_donttrim():
    '''When just trimming is called, relies on enumerating fresh files'''
    init_do_trim(do_trim=False, entry_via_pipelne=False)


def test_trim_adapters_pipelinentry_dotrim():
    '''With prior kraken run, but trimming set to false'''
    init_do_trim(do_trim=True, entry_via_pipelne=True)


def test_trim_adapters_pipelinentry_donttrim():
    '''With prior kraken run, with trimming set to true'''
    init_do_trim(do_trim=False, entry_via_pipelne=True)
