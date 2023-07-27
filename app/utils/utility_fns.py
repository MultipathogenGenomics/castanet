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
    with open(fpath, "r") as f:
        for l in f:
            if l[0] == ">":
                seqs.append(f"?{l}")
            else:
                seqs.append(l.replace("\n", ""))
    seqs_split = [i.split("\n") for i in "".join(seqs).split("?")]
    return [i for i in seqs_split if not i == [""]]


def save_fa(fpath, pat):
    with open(fpath, "w") as f:
        f.write(pat)


def get_reference_org(gt_file, seq_name, folder_stem):
    gt_table = pd.read_csv(gt_file)
    acc_id = gt_table[gt_table["Primary_accession"] ==
                      seq_name]["GenBank_accession"].item()
    ref_gb = DownloadGenBankFile(
        f"{folder_stem}consensus_data/GROUND_TRUTH_{gt_file}.gb", acc_id, "test@test.com")
    return [f">{acc_id}", str(ref_gb[acc_id].seq)]
