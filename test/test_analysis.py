import pytest
import os
import shutil

from test.utils import get_random_str, make_rand_dir, get_default_args
from app.src.analysis import Analysis
from app.src.generate_counts import run_counts
from app.src.map_reads_to_ref import run_map


def init_map(p):
    run_map(p)


def init_generate_counts(p):
    run_counts(p)


def init_analysis(start_with_bam=False):
    '''Init start files'''
    p = get_default_args()
    fstem = f"{p['SaveDir']}/{p['ExpName']}/"
    if os.path.exists(fstem):
        shutil.rmtree(fstem)
    os.mkdir(fstem)

    init_map(p)
    run_counts(p)

    clf = Analysis(p, start_with_bam)
    clf.main()
    assert os.stat(
        f'{p["SaveDir"]}/{p["ExpName"]}/{p["ExpName"]}_fullself.df.csv.gz').st_size > 2

    shutil.rmtree(fstem)


def test_analysis():
    init_analysis()


test_analysis()
