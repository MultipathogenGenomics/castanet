import os
import sys
import re
import random
import pandas as pd
from app.utils.utility_fns import read_fa
from app.utils.shell_cmds import loginfo, logerr, stoperr


class MappingRefConverter:
    def __init__(self, payload, sneaky_mode=True):
        self.payload = payload
        self.default_aggregation_val = "unaggregated"
        if not sneaky_mode:
            self.in_file = payload["InFile"]
            self.out_file = payload["OutFile"]
        else:
            self.in_file = self.payload["RefStem"]
            self.out_file = f"{payload['SaveDir']}/{payload['ExpName']}/ref.fa"
            self.payload["RefStem"] = self.out_file
        self.sneaky_mode = sneaky_mode

    def make_csv(self) -> pd.DataFrame:
        fixed_fasta = False
        if not self.in_file.endswith((".fa", ".fasta", ".fna")):
            stoperr(
                f"Input mapping reference file {self.in_file} is not a fasta file.")
        try:
            fastas = read_fa(self.in_file)
        except Exception as e:
            raise ValueError(
                f"Error reading fasta file {self.in_file}. Please ensure it is in valid fasta format.") from e

        for fasta in fastas:
            '''Check for leading ">"'''
            if not fasta[0].startswith(">"):
                if self.payload["FixFasta"]:
                    logerr(
                        f"Fasta entry {fasta[0]} does not have a header starting with '>'. I'm adding a header of '>{fasta[0].split()[0]}', which will be saved to your new fasta file.")
                    fasta[0] = f">{fasta[0]}"
                    fixed_fasta = True
                else:
                    stoperr(
                        f"Fasta entry {fasta[0]} does not have a header starting with '>'. Please check your file is a valid FASTA, or run with FixFasta=True to automatically replace invalid characters with 'n' and add missing '>' to headers.")
            '''Check for invalid characters in sequence'''
            tmp = set([i for i in fasta[1].lower().replace("a", "").replace(
                "t", "").replace("c", "").replace("g", "").replace("n", "")])
            if len(tmp) > 0:
                if self.payload["FixFasta"]:
                    logerr(
                        f"Fasta entry {fasta[0]} has non-ATCGN character/s in its sequence: {tmp}. I'm replacing these with 'n', which will be saved to your new fasta file.")
                    fasta[1] = re.sub(r'[^ATCGNatcgn]', 'n', fasta[1])
                    fixed_fasta = True
                else:
                    stoperr(
                        f"Fasta entry {fasta[0]} has non-ATCGN character/s in its sequence: {tmp}. Please check your file is a valid FASTA, or run with FixFasta=True to automatically replace invalid characters with 'n' and add missing '>' to headers.")

        if fixed_fasta:
            saved_file = self.in_file.split(".")[0] + "_fixed.fa"
            with open(saved_file, "w") as f:
                [f.write(f"{i[0]}\n{i[1]}\n") for i in fastas]
            loginfo(
                f"Some formatting issues were found in your input fasta file {self.in_file} and have been automatically fixed. The converted mapping reference has been saved to {saved_file}. Please check this file to ensure the formatting is correct, then rerun this function with the edited fasta file if you are happy to proceed.")
            return pd.DataFrame(), True

        agg_headers, descriptions, seqs, organisms, rmlst = [], [], [], [], []
        for fasta in fastas:
            if "bact0" in fasta[0].lower():  # Aggregaton to key with "bact0"
                probe_regexes = [
                    re.compile(r'bact[0-9]+_([A-Za-z]+)-[0-9]+[-_]([A-Za-z]+)'),
                    re.compile(r'bact[0-9]+_[0-9]+_([A-Za-z]+_[A-Za-z_]+)'),
                    re.compile(r'bact[0-9]+_([a-z]+_[a-z_]+)'),
                    re.compile(r'bact[0-9]+_([A-Za-z]+)-[0-9]+')
                ]

                def _pat_search(s):
                    '''Private function to return empty string instead of error when pattern is not matched.'''
                    try:
                        res = probe_regexes[0].findall(s)
                        if not res:
                            res = (probe_regexes[1].findall(s),)
                            if not res[0]:
                                has_cluster = re.search(r'_cluster_[0-9]+', s)
                                if has_cluster:
                                    pat = has_cluster[0]
                                    s = f"{s.replace(pat, '')}"

                                res = (probe_regexes[2].findall(s),)
                                if not res[0]:
                                    res = (probe_regexes[3].findall(s),)
                        if not res[0]:
                            return ''
                        name = '-'.join(res[0])

                        if name[-1] == "_":
                            # Fix for old probe set with random trailing _'s
                            name = name[:-1]

                    except Exception as e:
                        logerr(
                            f"Castanet couldn't parse one or more of your probe names. Please ensure you've converted it to Castanet format with the /convert_mapping_reference/ endpoint and that input format was consistent with the format expected (see documentation).\n{s}\n{e}")
                        return s
                    return name

                s=_pat_search(fasta[0][1:].lower())

                if any([s.startswith(x.lower()) for x in ["escherichia","klebsiella","enterobacter","shigella","serratia"]]):
                    s="enterobacteraciae"
                elif s in ['Streptococcus-pyogenes', 'Streptococcus-agalactiae']:
                    s = "streptococcus-agalactiae-pyogenes"
                elif s in ['streptococcus-pneumoniae','streptococcus-pseudopneumoniae', 'streptococcus-mitis', 'streptococcus-oralis']:
                    s = 'streptococcus-mitisgroup'

                match = re.findall(r"bact[0-9]*", fasta[0].lower())
                rmlst.append(match[0])
                agg_headers.append("rmlst-"+s)
                descriptions.append(
                    "_".join(fasta[0].split("_")[1:]).replace(",", ""))
                organisms.append("rmlst-"+s)

            elif "-segment" in fasta[0].lower():
                org = fasta[0].lower().split("-segment")[0][1:]
                seg = fasta[0].lower().split("-segment")[-1].split("_")[0]
                rmlstname = "segment-"+seg
                rmlst.append(rmlstname)
                agg_headers.append(org)
                organisms.append(org)
                descriptions.append(
                    "_".join(fasta[0].split("_")[1:]).replace(",", ""))
            else:
                if "segment" in fasta[0].lower() and "-segment" not in fasta[0].lower():
                    loginfo("If you would like segments of a virus to be aggregated please format target file fasta header as virusA-segmentB")
                rmlst.append("")
                if len(fasta[0].split("_")) < 2:
                    logerr(f"Mapping reference {fasta[0]} has no underscores, so will not aggregate with any other references! Please refer to documentation. "
                           f"I'm setting this to '{self.default_aggregation_val}'.")
                    agg_headers.append(self.default_aggregation_val)
                    descriptions.append(fasta[0].replace(">", ""))
                else:
                    agg_headers.append(fasta[0].split("_")[0].replace(">", ""))
                    descriptions.append(
                        "_".join(fasta[0].split("_")[1:]).replace(",", ""))
                organisms.append(fasta[0].split(
                    "_")[0].replace(">", "").split("-")[0])
            try:
                seqs.append(fasta[1])
            except IndexError as e:
                stoperr(
                    f"Fasta entry {fasta[0]} has no sequence associated with it. Please check your file is a valid FASTA.")
        df = pd.DataFrame({"organism": organisms, "probetype": agg_headers,
                           "description": descriptions, "sequence": seqs, "rmlst": rmlst})

        return df, False

    def input_checks(self, df) -> pd.DataFrame:
        '''Scans header organism and probetype values for disallowed characters. Stop if found and report to user.'''
        aggregation_headers = df["organism"].unique(
        ).tolist() + df["probetype"].unique().tolist()
        disallowed_chars = [" ", "/", "\\", ":", "@", "(", ")", "]", "[", ";", "#", "$", "%", "^", "&",
                            "*", "?", "\"", "<", ">", ","]
        errors = []
        for header in aggregation_headers:
            for char in disallowed_chars:
                if char in header:
                    errors.append([char, header])

        if len(errors) > 0:
            stoperr(f"Disallowed characters '{[i[0] for i in errors]}' found in probetype value '{[i[1] for i in errors]}'"
                    f" Please remove or replace, then try again.")

        return df

    def generate_hash(self, df) -> pd.DataFrame:
        df["key"] = df.apply(
            lambda x: f"{random.getrandbits(128)}", axis=1)
        return df

    def generate_fasta(self, df) -> list:
        fasta = []
        for _, row in df.iterrows():
            fasta.append(
                [f">{row['probetype']}_{row['key']}", row['sequence']])
        return fasta

    def save_output(self, df, fasta) -> None:
        if not self.sneaky_mode:
            '''Save CSV'''
            df[["organism", "probetype", "description", "key", "rmlst"]].to_csv(
                f"{self.out_file}", index=False)
        else:
            if not os.path.exists(f"{self.payload['SaveDir']}/{self.payload['ExpName']}/"):
                os.makedirs(
                    f"{self.payload['SaveDir']}/{self.payload['ExpName']}/")
            self.payload["MappingRefTable"] = f"{self.payload['SaveDir']}/{self.payload['ExpName']}/MappingRefTable.csv"
            df = df.applymap(lambda s: s.lower() if type(s) == str else s)
            df[["organism", "probetype", "description", "key", "rmlst"]].to_csv(
                self.payload["MappingRefTable"], index=False)
            with open(self.out_file, "w") as f:
                for header, seq in fasta:
                    f.write(f"{header}\n{seq}\n")

    def validate_user_csv(self, df):
        try:
            if df[["organism", "probetype", "key"]].isnull().values.any():
                stoperr(f"Your input MappingRefTable has empty values in the probetype and/or description columns. "
                        f"Castanet can't proceed as it needs names for each target we map to. "
                        f"Please manually edit these, or re-generate the mapping reference with the /convert_mapping_ref/ function.")
        except:
            stoperr(f"Castanet couldn't read your mapping reference table. Please ensure it is a CSV file with columns 'organism', 'probetype', and 'key'.")

    def join_seqs_to_df(self, df, users_df):
        out_df = users_df.copy()
        out_df["sequence"] = df["sequence"]
        return out_df

    def main(self):
        '''Convert an input CSV or FASTA mapping reference description file to a Castanet-compatible RefStem'''
        df, res = self.make_csv()
        if res:
            return "Please restart following correction of input fasta file as described in the log message above."

        if not self.sneaky_mode:
            loginfo(f"Converting mapping reference file: {self.in_file}")
            df = self.input_checks(df)
            df = self.generate_hash(df)
            fasta = ""

        if self.sneaky_mode:
            if os.path.isfile(self.payload["MappingRefTable"]):
                loginfo(
                    f'Parsing user supplied MappingRefTable: {self.payload["MappingRefTable"]}.')
                users_df = pd.read_csv(
                    self.payload["MappingRefTable"], index_col=None)
                self.validate_user_csv(users_df)
                df = self.join_seqs_to_df(df, users_df)
            else:
                loginfo(
                    f"No MappingRefTable supplied. Generating one automatically.")
                df = self.input_checks(df)
                df = self.generate_hash(df)
            fasta = self.generate_fasta(df)

        self.save_output(df, fasta)
        complete_msg = f"Conversion complete! Output saved to: {self.out_file}"
        if not self.sneaky_mode:
            loginfo(complete_msg)
            return complete_msg
        else:
            return self.payload
