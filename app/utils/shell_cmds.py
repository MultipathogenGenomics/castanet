import subprocess as sp
import sys


def shell(args, calling_fn="Misc shell function"):
    '''Call Bash shell with input string as argument'''
    whitelist = ["mafft"]
    _ = sp.Popen(args, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = _.communicate()
    if not any(x.lower() in whitelist for x in whitelist):
        sp_error_handler(
            out, err, f"ERROR: Something went wrong with CLI call: {calling_fn}, dumping STDOUT/STDERR to shell.")


def sp_error_handler(out, err, name):
    '''Kill program if error found; used in combo with POpen commands.'''
    if err != b"":
        raise SystemExit(f"***\n{name}:\nOut: {out}\nErr: {err}\n***")


def loginfo(s):
    '''Log Info statements to stderr'''
    sys.stderr.write(f'  INFO: {s}\n')


def logerr(s):
    '''Log Warning/Error statements to stderr'''
    sys.stderr.write(f'  Warning: {s}\n')


def stoperr(s, errcode=1):
    '''Call sys exit on finish or err'''
    status = 'Finished' if not errcode else 'Error'
    sys.stderr.write(f'  {status}: {s}\n')
    raise SystemError(f'  {status}: {s}\nErrcode: {errcode}')


def read_line(h):
    '''Parse line from input file'''
    return h.readline().decode("utf-8")


def make_dir(name):
    shell(f"mkdir -p {name}")
