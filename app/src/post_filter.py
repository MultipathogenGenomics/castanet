from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print


def run_post_filter(p):
    p["ExpDir"] = f"{p['ExpDir']}/"
    shell(f"samtools view -h {p['ExpDir']}{p['ExpDir']}.bam | python3 -m app.src.parse_bam -Mode filter -ExpDir {p['ExpDir']} -SeqName {p['ExpDir']} -FilterFile {p['ExpName']}_reads_to_drop.csv")
    end_sec_print(f"INFO: Post filter complete")
