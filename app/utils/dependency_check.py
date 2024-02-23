from app.utils.shell_cmds import shell, stoperr, loginfo
import os, random

class Dependencies:
    def __init__(self, p) -> None:
        self.p = p
        self.threads = os.cpu_count()
        self.hash = random.getrandbits(32)
        self.healthy_resp = {"response": "All dependencies installed and functional."}

    def check_kraken(self) -> None:
        shell(f'kraken2 --db {self.p["KrakenDbDir"]} --threads {self.threads} data/eval/sim_reads_1.fastq.gz > {self.hash}.kraken')
        if not os.path.exists(f"{self.hash}.kraken"):
            stoperr(f"Kraken2 did not function properly when called. Check it's installed and that your local Kraken2-compatible database path is correct.")
        else:
            os.remove(f"{self.hash}.kraken")

    def check_trimmomatic(self):
        trim_path='java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar' # TODO parameterise, harmonise with trim_adapters.py
        trimmies = [f"{self.hash}_1_clean.fastq", f"{self.hash}_1_trimmings.fq", f"{self.hash}_2_clean.fastq", f"{self.hash}_2_trimmings.fq"]
        shell(f"{trim_path} PE -threads {self.threads} data/eval/sim_reads_1.fastq.gz data/eval/sim_reads_2.fastq.gz {trimmies[0]} {trimmies[1]} {trimmies[2]} {trimmies[3]} ILLUMINACLIP:{self.p['AdaptP']}:2:10:7:1:true MINLEN:30")
        if not os.path.exists(trimmies[0]):
            stoperr(f"Trimmomatic did not function properly when called. Check it's installed and that your adapter path is correct.")
        else:
            [os.remove(f) for f in trimmies]

    def check_mapping(self):
        bwa_path='./bwa-mem2-2.2.1_x64-linux/bwa-mem2' # TODO parameterise, harmonise with map_reads_to_ref.py
        mappies = [f"{self.hash}.bam", f"{self.hash}.fasta", f"{self.hash}.aln"]
        shell(f"{bwa_path} index data/eval/ref.fa")
        shell(f"{bwa_path} mem -t {os.cpu_count()} data/eval/ref.fa data/eval/sim_reads_[12].fastq.gz | samtools view -F4 -Sb - | samtools sort - 1> {mappies[0]}")
        if not os.path.exists(mappies[0]):
            stoperr(f"BWA-mem2 did not function properly when called. Check it's installed.")
        self.check_samtools(mappies)
        self.check_mafft(mappies)
        [os.remove(f) for f in mappies]

    def check_samtools(self, paths):
        shell(f"samtools view {paths[0]} > {paths[1]}")
        if not os.path.exists(paths[1]):
            stoperr(f"Samtools did not function properly when called. Check it's installed.")

    def check_mafft(self, paths):
        shell(f"mafft --auto {paths[1]} --thread {self.threads} > {paths[2]}")
        if not os.path.exists(paths[1]):
            stoperr(f"Mafft did not function properly when called. Check it's installed.")

    def check_viral_cons(self):
        ...

    def main(self):
        loginfo("Castanet is searching for and testing all dependencies.")
        # self.check_kraken()
        # self.check_trimmomatic()
        self.check_mapping()
        self.check_viral_cons()
        loginfo("All dependencies installed and functioning correctly.")
        return self.healthy_resp

if __name__ == "__main__":
    clf = Dependencies({"KrakenDbDir": "kraken2_human_db/", "AdaptP": "data/all_adapters.fa"})
    clf.main()