import pytest
import os
import shutil
import random
import pickle

from test.utils import get_default_args, get_random_str, make_rand_dir, create_test_file
from app.utils.hash_files import hash_me, check_infile_hashes


def test_hash_me():
    '''Normal use'''
    fstem = make_rand_dir()
    fname = f"{fstem}/{get_random_str}.foo"
    create_test_file(fname)
    hash = hash_me(fname)
    assert len(hash) == 16
    shutil.rmtree(fstem)
    '''Non-existent file'''
    with pytest.raises(SystemError):
        hash_me("nonexistent_file")


def test_check_infile_hashes():
    ''''First time - should create 2 new hashes'''
    payload = get_default_args()
    exp_dir = f'{payload["SaveDir"]}/{payload["ExpName"]}'
    payload_new = check_infile_hashes(payload, exp_dir)
    for fname in payload["SeqNames"]:
        assert os.path.exists(f"{exp_dir}/hashes/{fname.split('/')[-1]}.p")
    assert payload == payload_new
    '''Second time - should run through the hash comparison'''
    payload = check_infile_hashes(payload, exp_dir)
    for fname in payload["SeqNames"]:
        assert os.path.exists(f"{exp_dir}/hashes/{fname.split('/')[-1]}.p")
    '''Modify a hash - should error when mismatch detected'''
    rand_hash = random.getrandbits(16)
    pickle.dump(rand_hash, open(
        f"{exp_dir}/hashes/{payload['SeqNames'][0].split('/')[-1]}.p", "wb"))
    with pytest.raises(SystemError):
        payload = check_infile_hashes(payload, exp_dir)
    shutil.rmtree(f'{payload["SaveDir"]}/{payload["ExpName"]}')


if __name__ == "__main__":
    test_hash_me()
    test_check_infile_hashes()
