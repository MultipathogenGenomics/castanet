import os
import hashlib
import pickle

from app.utils.utility_fns import enumerate_read_files, make_exp_dir
from app.utils.shell_cmds import stoperr


def hash_me(fname):
    md5 = hashlib.md5()
    try:
        with open(fname, 'rb') as f:
            for chunk in iter(lambda: f.read(8192 * md5.block_size), b''):
                md5.update(chunk)
        return md5.digest()
    except FileNotFoundError:
        stoperr(f"Input file not found when attempting to hash")


def check_infile_hashes(payload, exp_dir):
    '''If an ExpDir already exists, check if hashes exist and gate if input data changed.'''
    payload["SeqNames"] = enumerate_read_files(payload["ExpDir"])
    existing_hashes = {}
    if os.path.exists(exp_dir):
        try:
            for fname in payload["SeqNames"]:
                existing_hashes[fname] = pickle.load(
                    open(f"{exp_dir}/hashes/{fname.split('/')[-1]}.p", "rb"))
        except:
            stoperr(f"You're trying to run an experiment in a directory that already exists, but has no data hashes to compare against. "
                    f"This is a safety feature to stop you from overwriting data. "
                    f"Either change your save directory or delete the existing file.")

    '''Make experiment directory and hash new files'''
    make_exp_dir(exp_dir)

    new_hashes, failstate = {}, False
    for fname in payload["SeqNames"]:
        new_hashes[fname] = hash_me(fname)

    if len(existing_hashes) > 0:
        '''If hashes exist, compare to new ones'''
        status = "Previous hashes existed, but they're identical to the new ones, meaning data versions are the same"
        for fname in existing_hashes.keys():
            if not fname in new_hashes.keys():
                failstate = True
                continue
            if not existing_hashes[fname] == new_hashes[fname]:
                failstate = True
        if failstate:
            stoperr(f"You've tried to create a Castanet run to a folder (SaveDir) that already exists."
                    f"You can only do this if the input data (in DataFolder) are IDENTICAL to the previous run of with the same SaveDir"
                    f"This is a safety feature to prevent you from overwriting experimental data. Try changing the SaveDir and re-running.")
    else:
        '''Write new hashes if none existed before'''
        status = "No previous hashes existed, so Castanet has generated some"
        os.mkdir(f"{exp_dir}/hashes/")
        for fname in new_hashes.keys():
            pickle.dump(new_hashes[fname], open(
                f"{exp_dir}/hashes/{fname.split('/')[-1]}.p", "wb"))
    print(f"Completed experiment data hash check ({status}).")
    return payload


if __name__ == "__main__":
    print(hash_me(input("Type relative path of file to hash\n")))
