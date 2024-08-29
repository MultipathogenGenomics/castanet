import json


def write_input_params(payload):
    try:
        payload_for_sav = payload.copy()
        for key, val in payload_for_sav.items():
            '''Frozensets don't serialise to JSON'''
            if type(val) == frozenset:
                payload_for_sav[key] = list(val)
        with open(f"{payload['SaveDir']}/{payload['ExpName']}/run_parameters.json", "w") as f:
            f.write(json.dumps(payload_for_sav))
    except FileNotFoundError as e:
        print("INFO: Not logging input params as operating in batch mode")
