import pytest
import os
import shutil

from test.utils import get_random_str, make_rand_dir, get_default_args
from app.src.generate_counts import run_counts
from app.src.map_reads_to_ref import run_map
from app.utils.shell_cmds import shell


def init_map(p):
    run_map(p)


def init_generate_counts(start_with_bam, fake_files=False):
    '''Init start files'''
    p = get_default_args()
    fstem = f"{p['SaveDir']}/{p['ExpName']}/"
    if os.path.exists(fstem):
        shutil.rmtree(fstem)
    os.mkdir(fstem)

    '''Prepare input files'''
    if fake_files:
        p["ExpDir"] = get_random_str()
        with pytest.raises(SystemError):
            run_counts(p, start_with_bam)

    else:
        run_map(p)
        if start_with_bam:
            '''Test entrypoint where premade bam lives in ExpDir'''
            expdir_fstem = make_rand_dir()
            p['ExpDir'] = expdir_fstem
            shell(f"cp {fstem}test.bam {expdir_fstem}.bam")

        '''Call script to be tested'''
        run_counts(p, start_with_bam)
        '''Check main output form generate_counts (position counts csv) exists and is correct size'''
        assert os.stat(
            f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_PosCounts.csv").st_size > 2

    '''Clean up'''
    shutil.rmtree(fstem)
    if start_with_bam and not fake_files:
        shutil.rmtree(expdir_fstem)


def test_generate_counts_bamstart():
    '''Start with a BAM file (i.e. no preprocess or map)'''
    init_generate_counts(start_with_bam=True)


def test_generate_counts_fastqstart():
    '''Start mid-pipeline'''
    init_generate_counts(start_with_bam=False)


def test_generate_counts_wrong_n_bams():
    '''Test with multiple BAM files in indir or none'''
    init_generate_counts(start_with_bam=True, fake_files=True)
