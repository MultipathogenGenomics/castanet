import random
import string
import os


def get_random_str(n=8):
    '''Return random string of len n'''
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(n))


def make_rand_dir():
    '''Makes directory in test with random name'''
    fstem = f"./test/{get_random_str()}/"
    os.mkdir(fstem)
    return fstem


def create_test_file(fname):
    with open(fname, 'a'):
        os.utime(fname, None)


def get_default_args():
    return {
        "ExpDir": "./data/eval/",
        "ExpName": "test",
        "SaveDir": "./experiments",
        "RefStem": "data/eval/ref.fa",
        "SingleEndedReads": False,
        "MatchLength": 40,
        "DoTrimming": True,
        "TrimMinLen": 36,
        "DoKrakenPrefilter": True,
        "LineageFile": "data/ncbi_lineages_2023-06-15.csv.gz",
        "ExcludeIds": "9606",
        "RetainIds": "",
        "RetainNames": "",
        "ExcludeNames": "Homo",
        "DoConsensus": True,
        "ConsensusMinD": 10,
        "ConsensusCoverage": 30,
        "ConsensusMapQ": 1,
        "ConsensusCleanFiles": True,
        "GtFile": "",
        "GtOrg": "",
        "KrakenDbDir": "kraken2_human_db/",
        "KeepDups": True,
        "Clin": "",
        "DepthInf": "",
        "SamplesFile": "",
        "PostFilt": False,
        "AdaptP": "data/all_adapters.fa",
        "NThreads": os.cpu_count() - 1,
        "SeqNames": ["data/eval/sim_reads_1.fastq.gz", "data/eval/sim_reads_2.fastq.gz"]
    }
