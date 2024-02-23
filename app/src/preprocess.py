from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print


def run_kraken(p, api_entry=True):
    '''Call Kraken2 to remove unwanted reads'''
    if not api_entry:
        p = {
            "KrakenDbDir": "kraken2_human_db/",
            "NThreads": 4,
            "ExpDir": p.ExpDir,
            "SeqName": p.SeqName
        }
    shell(f'kraken2 --db {p["KrakenDbDir"]} --threads {p["NThreads"]} {p["ExpDir"]}/{p["SeqName"]}_1.fq.gz > experiments/{p["ExpName"]}/{p["SeqName"]}_1.kraken')
    end_sec_print(f"Kraken2 annotations complete")
