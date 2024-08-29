import pytest
import os
import shutil
import json

from test.utils import get_default_args
from app.utils.write_logs import write_input_params


def test_write_input_params():
    payload = get_default_args()
    fstem = f"{payload['SaveDir']}/{payload['ExpName']}/"
    if os.path.exists(fstem):
        shutil.rmtree(fstem)
    os.mkdir(fstem)
    fname = f"{fstem}/run_parameters.json"
    write_input_params(payload)
    assert os.path.exists(fname)
    with open(fname, "r") as f:
        fs = json.load(f)
    assert fs == payload
    shutil.rmtree(fstem)


if __name__ == "__main__":
    test_write_input_params()
