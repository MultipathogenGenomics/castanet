import os
from fastapi.encoders import jsonable_encoder
from app.utils.system_messages import end_sec_print
from app.utils.error_handlers import error_handler_api, stoperr
from app.utils.mapping_ref_convert import MappingRefConverter
from app.utils.argparsers import parse_arguments_cmr


def defaults():
    return {
        "InFile": "./my_mapping_ref.fasta",
        "OutFile": "./my_mapping_ref_table.csv"
    }


def populate_request(payload, parser):
    '''Overwrite default request object with argparser args'''
    for key, val in parser.__dict__.items():
        payload[key] = val
    return payload


def tests(payload):
    '''KrakenDbDir'''
    if not os.path.isfile(payload['InFile']):
        stoperr(f"Your InFile doesn't seem to exist. Check spelling and path.")


def main():
    parser = parse_arguments_cmr()
    payload = populate_request(defaults(), parser.parse_args())
    tests(payload)
    end_sec_print(f"Calling Castanet Lite with following arguments: {payload}")
    try:
        clf = MappingRefConverter(payload, sneaky_mode=False)
        return clf.main()
    except Exception as ex:
        return error_handler_api(ex)


if __name__ == "__main__":
    main()
