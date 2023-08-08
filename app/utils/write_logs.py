import json


def write_input_params(payload):
    with open(f"./experiments/{payload['ExpName']}/run_parameters.json", "w") as f:
        f.write(json.dumps(payload))
