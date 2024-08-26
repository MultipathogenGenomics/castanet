import argparse


def parse_args_bam_parse():
    '''Parse args for filter and parse BAM functions. Programmatic only - no user interaction should be necessary.'''
    parser = argparse.ArgumentParser(
        description=__doc__, epilog='Do you want to parse or filter?')
    parser.add_argument('-ExpDir', required=False,
                        default="parse", help="Pass API arg in via shell.")
    parser.add_argument('-SingleEnded', required=True,
                        default="parse", help="Pass API arg in via shell.")
    parser.add_argument('-SeqName', required=True,
                        default="parse", help="Pass API arg in via shell.")
    parser.add_argument('-ExpName', required=False,
                        default="parse", help="Pass API arg in via shell.")
    parser.add_argument('-FilterFile', required=False,
                        default="parse", help="Pass API arg in via shell.")
    parser.add_argument('-Mode', required=True,
                        default="parse", help="parse or filter.")
    parser.add_argument('-SaveDir', required=True,
                        default="./experiments", help="Pass API arg in via shell.")
    parser.add_argument('-MatchLength', required=True,
                        default="./experiments", help="Pass API arg in via shell.")
    return parser.parse_args()
