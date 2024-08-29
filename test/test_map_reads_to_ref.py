import pytest
import os
import shutil

from test.utils import get_random_str, make_rand_dir, create_test_file, get_default_args
from app.src.map_reads_to_ref import run_map
from app.src.trim_adapters import run_trim


def init_preprocessing(p):
    '''Call preprocessing (trim)'''
    p["DoTrimming"] = True
    run_trim(p)


def init_run_map(do_preprocessing, fake_files=False):
    '''Init'''
    p = get_default_args()
    fstem = f"{p['SaveDir']}/{p['ExpName']}/"
    os.mkdir(fstem)
    if do_preprocessing:
        '''Run preprocessing'''
        init_preprocessing(p)
    if fake_files:
        p['ExpDir'] = get_random_str()
    '''Call script to test'''
    if fake_files:
        '''Fake input files should raise stoperr'''
        with pytest.raises(SystemError):
            run_map(p)
    else:
        '''Real files should lead to a real BAM file output'''
        run_map(p)
        assert os.stat(
            f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam").st_size > 2
        if do_preprocessing:
            '''Check deleted the temp files if did preprocessing'''
            assert not os.path.exists(
                f"rm {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_1_clean.fastq")
    shutil.rmtree(fstem)


def test_map_no_preprocessing():
    '''Test without preprocessing, i.e. map test data files'''
    init_run_map(do_preprocessing=False)


def test_map_with_preprocessing():
    '''Test with preprocessing, i.e. reads trimmed files'''
    init_run_map(do_preprocessing=True)


def test_map_no_input_fs():
    '''Test with non-existnt input files'''
    init_run_map(do_preprocessing=False, fake_files=True)
