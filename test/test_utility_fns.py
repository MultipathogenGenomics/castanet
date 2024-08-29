import pytest
import os
import shutil

from test.utils import get_random_str, make_rand_dir, create_test_file
from app.utils.utility_fns import (make_exp_dir, get_gene_orgid, read_fa, save_fa,
                                   trim_long_fpaths, enumerate_bam_files, enumerate_read_files)


def test_make_exp_dir():
    dir_name = f"./test/{get_random_str()}"
    make_exp_dir(dir_name)
    assert os.path.exists(dir_name)
    os.rmdir(dir_name)


def test_get_gene_orgid():
    target_id = "hepatitiserat_1678145_5_cluster10"
    out = get_gene_orgid(target_id)
    assert out[0] == "hepatitiserat" and out[0] == out[1]
    target_id = "bact000001_leptospira_genericstrai"
    out = get_gene_orgid(target_id)
    assert out[0] != out[1]


def test_read_fa():
    fa = read_fa("./data/eval/ref.fa")
    assert len(fa) == 1 and fa[0][0][0] == ">" and fa[0][1][0] == "G"


def test_save_fa():
    fa = [[">seq1", "ATCG"], [">seq2", "ATCG"]]
    fa_fname = f"./test/{get_random_str()}.fasta"
    with open(fa_fname, "w") as f:
        [f.write(f"{i[0]}\n{i[1]}\n") for i in fa]
    assert os.path.exists(fa_fname)
    fa_read = read_fa(fa_fname)
    assert len(
        fa_read) == 2 and fa_read[0][0][0] == ">" and fa_read[0][1] == "ATCG"
    os.remove(fa_fname)


def test_trim_long_fpaths():
    max_len = 100
    key_lens = [max_len, max_len+1]
    for key_len in key_lens:
        key = ''.join(get_random_str(key_len))
        trimmed_key = trim_long_fpaths(key, 100)
        if key_len > max_len:
            assert len(trimmed_key) == max_len
        else:
            assert len(trimmed_key) == len(key)


def test_enumerate_read_files():
    '''Correct usage scenario'''
    fstem = make_rand_dir()
    correct_f_structure = [
        f"{fstem}/read_1.fastq.gz", f"{fstem}/read_2.fastq.gz"]
    for f in correct_f_structure:
        create_test_file(f)
    fs = enumerate_read_files(fstem)
    assert len(fs) == 2
    shutil.rmtree(fstem)

    '''Non-existing folder'''
    nonexistant_fstem = get_random_str()
    with pytest.raises(SystemError):
        enumerate_read_files(nonexistant_fstem)

    '''Extra files ending in fq in directory'''
    fstem = make_rand_dir()
    incorrect_f_structure = [f"{fstem}/read_1.fastq.gz",
                             f"{fstem}/read_2.fastq.gz", f"{fstem}/read_3.fastq.gz"]
    for f in incorrect_f_structure:
        create_test_file(f)
    with pytest.raises(SystemError):
        fs = enumerate_read_files(fstem)
    os.remove(incorrect_f_structure[2])

    '''Extra files not of fq type in directory'''
    non_fa_file = incorrect_f_structure[0].replace(".fastq.gz", ".foo")
    create_test_file(non_fa_file)
    fs = enumerate_read_files(fstem)
    assert len(fs) == 2
    shutil.rmtree(fstem)


def test_enumerate_bam_files():
    '''Correct usage'''
    fstem = make_rand_dir()
    bam_fname = f"{fstem}/foo.bam"
    create_test_file(bam_fname)
    fs = enumerate_bam_files(fstem)
    assert fs == bam_fname

    '''Two BAMs'''
    extra_bam_fname = bam_fname.replace(".bam", "2.bam")
    create_test_file(extra_bam_fname)
    with pytest.raises(AssertionError):
        enumerate_bam_files(fstem)
    os.remove(extra_bam_fname)

    '''Additional non-bam file'''
    junk_fname = bam_fname.replace(".bam", ".foo")
    create_test_file(junk_fname)
    fs = enumerate_bam_files(fstem)
    assert fs == bam_fname
    shutil.rmtree(fstem)
