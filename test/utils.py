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
