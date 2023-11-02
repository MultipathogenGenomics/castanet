import os
import re
import pandas as pd


class ProbeFileGen:
    '''
    Reads in all FASTA files in a folder, cleans headers, joins seqs and QC's
    file. Generates probe length file, then saves unified probes. Header format:
    >Species_[some unique string, no whitespace]-[…description, optional…]
    '''

    def __init__(self, a) -> None:
        self.fstem = a["OutFolder"]
        self.in_stem = a['InputFolder']
        self.out_fname = f"{self.fstem}{a['OutFileName']}.fasta"
        self.plen_fname = f"{self.fstem}{a['OutFileName']}.csv"
        self.files = os.listdir(self.in_stem)
        self.all_seqs = []
        self.master_seq_counter = 0
        self.stop_words = [[".mafft_consensus", ""], ["mafft", ""], ["consensus", ""], ['"', ""], ["E.coli", "Escherichia-coli"],
                           [",", ""], [" ", "-"], [".mafft", ""], ["_TRUE",
                                                                   ""], [".fst", ""], ["__", ""], ["--", "-"], ["/", "-"],
                           ["(", ""], [")", ""], [":", "-"], [".", "-"]]
        self.split_names = [
            "enterovirus", "influenza"
        ]

    def header_cleaner(self, header) -> str:
        '''Remove stop words, re-organise headers with rMLST gene first with primary org'''

        for i in self.stop_words:
            '''Stop word removal'''
            header = header.replace(i[0], i[1])

        if header == "":
            '''Blanks will get removed later'''
            return ""

        if header[0:5].lower() == ">bact":
            '''If rMLST gene first, reorganise...'''
            species_match = re.findall(
                r"_[A-Z]{1}[a-z]*-[0-9]*_[a-z]*", header.replace("|", "_"))
            if species_match:
                leading_string = re.sub(
                    r'[0-9]', '', species_match[0]).replace('_', '')
                if len([i for i in leading_string.split('-') if not i == '']) == 1:
                    leading_string = f'{leading_string.replace("_","")}GenericStrain'
                header = f">{leading_string.replace('-','_')}_{header.replace('>', '')}"
            else:
                raise ValueError(f"I COULDN'T PROCESS HEADER: {header}")

        elif "bact0" in header.lower() and not header[0:5].lower() == ">bact":
            header = header.replace("-", "_")

        for j in self.split_names:
            '''SPECIFIC TO DEVELOPER'S SET - Convert EnterovirusX/InfluenzaY to _X'''
            if j in header and not "bact" in header.lower():
                header = header.lower()
                header = f'>{j}_{header.split(j)[1].replace("_","").replace("-","")}_{j}{"".join(header.split(j)[2:])}'

        while header[-1] == "_":
            '''Fix random amount of trailing _'s'''
            header = header[:-1]

        return header.lower().replace("--", "-").replace("|", "-")

    def qc(self) -> None:
        '''Raise error if duplicate entries or gaps found'''
        self.all_seqs = [i for i in self.all_seqs if not i == ["", ""]]
        seq_names, seqs = [i[0] for i in self.all_seqs], [i[1]
                                                          for i in self.all_seqs]
        deduped_names = set(seq_names)
        if len(deduped_names) != len(seq_names):
            import collections
            raise ValueError(
                f"I found duplicate probe name/s: {[i for i, c in collections.Counter(seq_names).items() if c > 1]}. Please fix before re-trying.")
        if len(seq_names) != len(seqs):
            raise ValueError(
                f"Unequal number of sequences and headers. Please fix and run again.")

    def generate_probe_lengths(self) -> None:
        '''Generate CSV file with probe ids/lens'''
        plens = [{"target_id": i[0].replace(
            ">", ""), "target_len": len(i[1])} for i in self.all_seqs]
        df = pd.DataFrame(plens)
        df = df.sort_values(by="target_id")
        df.to_csv(self.plen_fname, index=False)

    def main(self) -> None:
        '''Read files in, the call cleaning, QC and output generators'''
        for file in self.files:
            if file == self.out_fname.split("/")[-1]:
                '''Don't include previous output'''
                continue

            seqs = []
            with open(f"{self.in_stem}{file}", "r") as f:
                [seqs.append(i.replace("\n", "")) for i in f]

            is_first_code_block = False
            current_seq_block, current_header = "", ""
            clean_seqs = []

            '''Iterate over seqs, sort and concat incomplete sequences. Call text cleaner on header lines'''
            for seq in seqs:
                if len(seq) == 0:
                    """Empty line, ignore"""
                    continue

                if seq[0] == ">":
                    """line = header"""
                    self.master_seq_counter += 1
                    clean_seqs.append([self.header_cleaner(
                        current_header), current_seq_block.upper()])
                    current_header = seq
                    is_first_code_block = True
                    current_seq_block = ""
                    continue

                elif is_first_code_block:
                    """line = initial chunk of sequence"""
                    is_first_code_block = False
                    current_seq_block = seq
                    continue
                else:
                    """line = not intitial chunk of sequence"""
                    current_seq_block = current_seq_block + seq

            print(f"Reading in {len(clean_seqs)} probes from file: {file}")
            self.all_seqs.append(clean_seqs)

        self.all_seqs = [
            i for sublist in self.all_seqs for i in sublist if not i == ""]
        self.qc()
        self.generate_probe_lengths()

        print(
            f"Aggregated {len(self.all_seqs)} probes into output file: {self.out_fname}.")
        with open(self.out_fname, "w") as f:
            [f.write(f"{i[0]}\n{i[1]}\n")
             for i in self.all_seqs if not i[0] == ""]


if __name__ == "__main__":
    '''Provide input seq name'''
    clf = ProbeFileGen("data/2023_panel/",
                       "data/2023_panel/raw_seqs/", "all_seqs")
    clf.main()
