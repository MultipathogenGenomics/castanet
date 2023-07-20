import json
from app.utils.shell_cmds import make_dir


def write_input_params(payload):
    make_dir(f"./experiments/{payload['ExpName']}/")
    with open(f"./experiments/{payload['ExpName']}/run_parameters.json", "w") as f:
        f.write(json.dumps(payload))
