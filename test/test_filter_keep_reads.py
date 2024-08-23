import pytest
import os
import shutil

from test.utils import get_random_str, make_rand_dir, create_test_file, get_default_args
from app.src.filter_keep_reads import FilterKeepReads
from app.src.preprocess import run_kraken

def test_filter_keep_reads():
    p = get_default_args()
    '''Test empty initialisation'''
    with pytest.raises(TypeError):
        clf = FilterKeepReads()
    '''Test breaks with empty/absent kraken file'''
    with pytest.raises(SystemError):
        clf = FilterKeepReads(p)
    '''Test the `DoKrakenPrefilter` switch'''
    p["DoKrakenPrefilter"] = False
    if os.path.exists(f'{p["SaveDir"]}/{p["ExpName"]}'):
        shutil.rmtree(f'{p["SaveDir"]}/{p["ExpName"]}')
    os.mkdir(f'{p["SaveDir"]}/{p["ExpName"]}')
    run_kraken(p)
    assert not os.path.exists(f'{p["SaveDir"]}/{p["ExpName"]}/{p["ExpName"]}.kraken')



if __name__ == "__main__":
    test_filter_keep_reads()