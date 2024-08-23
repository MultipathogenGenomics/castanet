import pytest
import os
import shutil

from test.utils import get_random_str, make_rand_dir, create_test_file, get_default_args
from app.src.trim_adapters import run_trim


def test_trim_adapters_trimentry_dotrim():
    '''When just trimming is called, relies on enumerating fresh files'''
    p = get_default_args()
    fstem = f"{p['SaveDir']}/{p['ExpName']}/"
    os.mkdir(fstem)
    run_trim(p)
    '''Check output has size > 2b'''
    trim_out = f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_1_clean.fastq"
    assert os.stat(trim_out).st_size > 2
    shutil.rmtree(fstem)


def test_trim_adapters_trimentry_donttrim():
    '''When just trimming is called, relies on enumerating fresh files'''
    p = get_default_args()
    fstem = f"{p['SaveDir']}/{p['ExpName']}/"
    p["DoTrimming"] = False
    os.mkdir(fstem)
    run_trim(p)
    '''Check output has size > 2b'''
    trim_out = f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_1_clean.fastq"
    assert os.stat(trim_out).st_size > 2
    shutil.rmtree(fstem)


def test_trim_adapters_notrim():
    '''With prior kraken run, but trimming set to false'''
    p = get_default_args()
    p["DoTrimming"] = False


if __name__ == "__main__":
    test_trim_adapters_trimentry_donttrim()
