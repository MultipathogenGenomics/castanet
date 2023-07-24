import os


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
