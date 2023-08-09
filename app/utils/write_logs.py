import json


def write_input_params(payload):
    try:
        with open(f"./experiments/{payload['ExpName']}/run_parameters.json", "w") as f:
            f.write(json.dumps(payload))
    except FileNotFoundError as e:
        print("INFO: Not logging input params as operating in batch mode")
