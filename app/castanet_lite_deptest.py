import os
from fastapi.encoders import jsonable_encoder
from app.utils.system_messages import end_sec_print
from app.utils.error_handlers import error_handler_api, stoperr
from app.utils.dependency_check import Dependencies
from app.utils.argparsers import parse_arguments_deptest


def defaults():
    return {
        "AdaptP": "data/all_adapters.fa",
        "KrakenDbDir": "kraken2_human_db/"
    }


def populate_request(payload, parser):
    '''Overwrite default request object with argparser args'''
    for key, val in parser.__dict__.items():
        payload[key] = val
    return payload


def tests(payload):
    '''KrakenDbDir'''
    if not os.path.exists(payload['KrakenDbDir']):
        stoperr(f"Your KrakenDbDir doesn't seem to exist. Check spelling and path.")


def main():
    parser = parse_arguments_deptest()
    payload = populate_request(defaults(), parser.parse_args())
    tests(payload)
    end_sec_print(f"Calling Castanet Lite with following arguments: {payload}")
    try:
        clf = Dependencies(jsonable_encoder(payload))
        return clf.main()
    except Exception as ex:
        return error_handler_api(ex)


if __name__ == "__main__":
    main()
