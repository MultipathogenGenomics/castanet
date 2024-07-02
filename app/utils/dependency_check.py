from app.utils.shell_cmds import shell, loginfo
from app.utils.error_handlers import error_handler_cli
import os
import random


class Dependencies:
    def __init__(self, p) -> None:
        self.p = p
        self.fol = "test_files/"
        self.threads = os.cpu_count()
        self.hash = random.getrandbits(32)
        self.healthy_resp = "All dependencies installed and functional."

    def check_kraken(self) -> None:
        loginfo("Testing Kraken2")
        out = shell(
            f'kraken2 --db {self.p["KrakenDbDir"]} --paired --threads {self.threads} --output {self.fol}{self.hash}.kraken data/eval/sim_reads_[12].fastq.gz', is_test=True)
        error_handler_cli(out, f"{self.fol}{self.hash}.kraken", "kraken")
        os.remove(f"{self.fol}{self.hash}.kraken")

    def check_trimmomatic(self):
        loginfo("Testing Trimmomatic")
        out = shell("java", is_test=True)
        error_handler_cli(out, "dummy_outf", "java", test_out_f=False)

        trim_path = 'java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar'
        trimmies = [f"{self.fol}{self.hash}_1_clean.fastq", f"{self.fol}{self.hash}_1_trimmings.fq",
                    f"{self.fol}{self.hash}_2_clean.fastq", f"{self.fol}{self.hash}_2_trimmings.fq"]
        out = shell(
            f"{trim_path} PE -threads {self.threads} data/eval/sim_reads_1.fastq.gz data/eval/sim_reads_2.fastq.gz {trimmies[0]} {trimmies[1]} {trimmies[2]} {trimmies[3]} ILLUMINACLIP:{self.p['AdaptP']}:2:10:7:1:true MINLEN:30", is_test=True)
        error_handler_cli(out, trimmies[0], "trimmomatic", test_f_size=True)
        [os.remove(f) for f in trimmies]

    def check_mapping(self):
        loginfo("Testing BWA-mem2")
        mappies = [f"{self.fol}{self.hash}.bam", f"{self.fol}{self.hash}.fasta",
                   f"{self.fol}{self.hash}.aln", f"{self.fol}{self.hash}_cons.fasta"]
        out = shell(f"bwa-mem2 index data/eval/ref.fa", is_test=True)
        shell(
            f"bwa-mem2 mem -t {os.cpu_count()} data/eval/ref.fa data/eval/sim_reads_[12].fastq.gz | samtools view -F4 -Sb - | samtools sort - 1> {mappies[0]}")
        error_handler_cli(out, mappies[0], "bwa-mem2", test_f_size=True)
        self.check_samtools(mappies)
        self.check_mafft(mappies)
        self.check_viral_cons(mappies)
        [os.remove(f) for f in mappies]

    def check_samtools(self, paths):
        loginfo("Testing Samtools")
        out = shell("samtools", is_test=True)
        shell(f"samtools fasta {paths[0]} > {paths[1]}")
        error_handler_cli(out, paths[1], "samtools")

    def check_mafft(self, paths):
        loginfo("Testing Mafft")
        out = shell(f"mafft --auto {paths[1]} > {paths[2]}", is_test=True)
        error_handler_cli(out, paths[2], "mafft")

    def check_viral_cons(self, paths):
        loginfo("Testing ViralConsensus")
        out = shell(
            f"viral_consensus -i {paths[0]} -r data/eval/ref.fa -o {paths[3]} --min_depth 1", is_test=True)
        error_handler_cli(out, paths[3], "viral_consensus", test_f_size=True)

    def main(self):
        loginfo("Castanet is searching for and testing all dependencies.")
        loginfo(
            f"Dependencies will be tested by processing dummy data: output will be saved to {self.fol}, then cleaned if all functions as expected.")
        shell(f"mkdir {self.fol}")
        self.check_kraken()
        self.check_trimmomatic()
        self.check_mapping()
        loginfo("All dependencies installed and functioning correctly.")
        shell(f"rmdir {self.fol}")
        return self.healthy_resp


if __name__ == "__main__":
    clf = Dependencies({"KrakenDbDir": "kraken2_human_db/",
                       "AdaptP": "data/all_adapters.fa"})
    clf.main()
