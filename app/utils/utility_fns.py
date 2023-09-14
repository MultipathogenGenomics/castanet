import os
import pandas as pd
from app.utils.get_genbank import DownloadGenBankFile


def make_exp_dir(ExpName):
    '''Checks if dir exists; creates if not'''
    if not os.path.exists(f"experiments/{ExpName}"):
        os.makedirs(f"experiments/{ExpName}")


def get_gene_orgid(target_id):
    '''Find gene and orgid at specific ref; return (gene, orgid).'''
    parts = target_id.split('_')
    return (parts[0], parts[-1] if parts[0].startswith('BACT') else parts[0])


def read_fa(fpath):
    '''Read fasta file to list of lists in format [[>name, seq], [...]]'''
    seqs = []
    try:
        with open(fpath, "r") as f:
            for l in f:
                if l[0] == ">":
                    seqs.append(f"?{l}")
                else:
                    seqs.append(l.replace("\n", ""))
    except:
        '''Current release of ViralConsensus will sometimes dump binary into the output file - this is a temp fix'''
        print(
            f"*** WARNING: Encoding on your input fasta file ({fpath}) is messed up. Attempting to rescue it...")
        seqs = []  # RESET SEQS as it doesn't get wiped by transition to exception
        bases = ["A", "T", "C", "G", "N", "-"]
        with open(fpath, "r", encoding="ISO-8859-1") as f:
            for l in f.readlines():
                if l[0] not in bases:
                    seqs.append(f"?{l}")
                else:
                    seqs.append(l)

        for i in range(len(seqs)):
            if seqs[i][0] == "?":
                continue
            else:
                seqs[i] = "".join([j for j in seqs[i] if j in bases])

    seqs_split = [i.split("\n") for i in "".join(seqs).split("?")]
    return [i for i in seqs_split if not i == [""]]


def save_fa(fpath, pat):
    with open(fpath, "w") as f:
        f.write(pat)


def get_reference_org(gt_file, seq_name, folder_stem):
    gt_table = pd.read_csv(gt_file)
    acc_id = gt_table[gt_table["Primary_accession"] ==
                      seq_name]["GenBank_accession"].item()
    if type(acc_id) != str:
        print(f"WARNING: Reference has no ground truth genome sequence!")
        return [f">NO REFERENCE AVAILABLE", "AAAAA"]
    ref_gb = DownloadGenBankFile(
        f"{folder_stem}consensus_data/GROUND_TRUTH_{gt_file}.gb", acc_id, "test@test.com")
    return [f">{acc_id}", str(ref_gb[acc_id].seq)]


def trim_long_fpaths(key):
    '''Curtail very long probe names that can't be used as folder names'''
    if len(key) > 100:
        return key[0:100]
    else:
        return key
