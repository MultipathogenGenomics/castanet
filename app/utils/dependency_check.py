from app.utils.shell_cmds import shell, stoperr, loginfo
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
        if not "Loading database information" in out:
            stoperr(
                f"Kraken2 doesn't seem to be installed. Please check the Castanet readme for installation instructions.")
        if not os.path.exists(f"{self.fol}{self.hash}.kraken"):
            stoperr(f"Kraken2 seems to be installed but didn't produce output. Check the Kraken2 database downloaded and unpacked correctly.")
        else:
            os.remove(f"{self.fol}{self.hash}.kraken")

    def check_trimmomatic(self):
        loginfo("Testing Trimmomatic")
        out = shell("java", is_test=True)
        if not "Usage: java [-options] class" in out:
            stoperr(
                "Java doesn't seem to be installed. Please check the Castanet readme for installation instructions.")

        trim_path = 'java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar'
        trimmies = [f"{self.fol}{self.hash}_1_clean.fastq", f"{self.fol}{self.hash}_1_trimmings.fq",
                    f"{self.fol}{self.hash}_2_clean.fastq", f"{self.fol}{self.hash}_2_trimmings.fq"]
        out = shell(
            f"{trim_path} PE -threads {self.threads} data/eval/sim_reads_1.fastq.gz data/eval/sim_reads_2.fastq.gz {trimmies[0]} {trimmies[1]} {trimmies[2]} {trimmies[3]} ILLUMINACLIP:{self.p['AdaptP']}:2:10:7:1:true MINLEN:30", is_test=True)
        if not "TrimmomaticPE: Started with arguments:" in out:
            stoperr(
                f"Trimmomatic doesn't seem to be installed. Please check the Castanet readme for installation instructions")
        if not os.path.exists(trimmies[0]):
            stoperr(
                f"Trimmomatic seems to be installed but didn't produce output. Please check your installation manually.")
        else:
            [os.remove(f) for f in trimmies]

    def check_mapping(self):
        loginfo("Testing BWA-mem2")
        # TODO parameterise, harmonise with map_reads_to_ref.py
        mappies = [f"{self.fol}{self.hash}.bam", f"{self.fol}{self.hash}.fasta",
                   f"{self.fol}{self.hash}.aln", f"{self.fol}{self.hash}_cons.fasta"]
        out = shell(f"bwa-mem2 index data/eval/ref.fa", is_test=True)
        if not "Looking to launch executable" in out:
            stoperr(
                f"BWA-mem2 doesn't seem to be installed. Please check the Castanet readme for installation instructions.")

        shell(
            f"bwa-mem2 mem -t {os.cpu_count()} data/eval/ref.fa data/eval/sim_reads_[12].fastq.gz | samtools view -F4 -Sb - | samtools sort - 1> {mappies[0]}")
        if not os.path.exists(mappies[0]):
            stoperr(
                f"BWA-mem2 seems to be installed but didn't produce output. Please check your installation manually.")
        self.check_samtools(mappies)
        self.check_mafft(mappies)
        self.check_viral_cons(mappies)
        [os.remove(f) for f in mappies]

    def check_samtools(self, paths):
        loginfo("Testing Samtools")
        out = shell("samtools", is_test=True)
        if not "Program: samtools" in out:
            stoperr(
                f"Samtools doesn't seem to be installed. Please check the Castanet readme for installation instructions.")
        shell(f"samtools fasta {paths[0]} > {paths[1]}")
        if not os.path.exists(paths[1]):
            stoperr(
                f"Samtools did not function properly when called. Please check your installation manually.")

    def check_mafft(self, paths):
        loginfo("Testing Mafft")
        out = shell(f"mafft --auto {paths[1]} > {paths[2]}", is_test=True)
        if not "generating a scoring matrix for nucleotide" in out:
            stoperr(
                f"Mafft doesn't seem to be installed. Please check the Castanet readme for installation instructions.")
        if not os.path.exists(paths[1]):
            stoperr(
                f"Mafft did not function properly when called. Please check your installation manually.")

    def check_viral_cons(self, paths):
        loginfo("Testing ViralConsensus")
        out = shell(f"./ViralConsensus/viral_consensus", is_test=True)
        if not "Missing required argument:" in out:
            stoperr(
                f"ViralConsensus doesn't seem to be installed. Please check the Castanet readme for installation instructions.")
        shell(
            f"./ViralConsensus/viral_consensus -i {paths[0]} -r data/eval/ref.fa -o {paths[3]} --min_depth 1")
        if not os.path.exists(paths[3]):
            stoperr(
                f"ViralConsensus did not function properly when called. Check it's installed.")

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
