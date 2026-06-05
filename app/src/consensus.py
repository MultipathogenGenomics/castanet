import os
import io
import numpy as np
import pandas as pd
import pickle as p
from Bio import AlignIO
from collections import Counter
import plotly.express as px
import re

from app.utils.timer import timing
from app.utils.shell_cmds import shell, make_dir, loginfo, stoperr, logerr
from app.utils.utility_fns import read_fa, save_fa
from app.utils.fnames import get_consensus_fnames
from app.utils.system_messages import end_sec_print
from app.utils.basic_cli_calls import (
    samtools_index, bwa_index, find_and_delete, rm, samtools_read_num)
from app.utils.error_handlers import error_handler_cli
from app.utils.similarity_graph import call_graph

import warnings
# Pandas zero div errors
warnings.simplefilter(action='ignore', category=RuntimeWarning)


class Consensus:
    '''Take all targets in one probetype/species aggregation, call consensus for each,
    flatten consensuses into single sequence.'''

    def __init__(self, payload, start_with_bam) -> None:
        self.a = payload
        self.a['ConsensusCoverage'] = 10  # Hard-coded from V9.0
        self.a["folder_stem"] = f"{self.a['SaveDir']}/{self.a['ExpName']}/"
        self.target_consensuses = {}
        self.insufficient_coverage_orgs = []
        self.refs = [[i[0][0:101], i[1]] for i in read_fa(self.a["RefStem"])]
        self.probe_names = pd.read_csv(
            f"{self.a['SaveDir']}/{self.a['ExpName']}/probe_aggregation.csv")
        self.fnames = get_consensus_fnames(self.a)
        try:
            self.grouped_reads = p.load(
                open(self.fnames["grouped_reads"], "rb"))
        except:
            stoperr(f"Failed to load groupedreads.p. To save space, this file is DELTED after running the consensus algorithm if debugmode = False. Please re-run the pipeline, selecting debug mode on if you wish to run the consensus algorithm several times.")
        self.grouped_reads_keep = {}
        self.subconsensuses = {}
        if start_with_bam:
            self.fnames['master_bam'] = f"{self.a['ExpDir']}/{[i for i in os.listdir(self.a['ExpDir']) if i[-4:] == '.bam'][0]}"
        self.eval_stats, self.naive_consensuses, self.coverage = {}, {}, None
        self.lut = pd.read_csv(self.a["MappingRefTable"], index_col=False)
        make_dir(f"{self.a['folder_stem']}consensus_data/")
        make_dir(f"{self.a['folder_stem']}consensus_sequences/")

    def filter_bam(self, tar_name) -> None:
        '''Filter bam to specific target, call consensus sequence for sam alignment records, grouped by target'''
        try:
            '''Don't give the user the hashed header name, it will only upset them'''
            probename="_".join(tar_name.split("_")[:-1]).lower()
            subset = self.lut.loc[self.lut["probetype"] == probename,"rmlst"]
            if not subset.empty:
                hgroup = subset.iloc[0]
                if pd.notna(hgroup):
                    sneaky_name = f'{tar_name.split("_")[0]}_{hgroup}'
                else:
                    sneaky_name = f'{tar_name.split("_")[0]}_{self.lut[self.lut["key"] == tar_name.split("_")[-1]]["description"].item()[0:100]}'
        except TypeError:
            '''If user has somehow broken fasta header'''
            sneaky_name = f'{tar_name.split("_")[0]}'
        except ValueError:
            stoperr(f"Castanet doesn't yet support calling a consensus sequence without you supplying a valid mapping reference table.")
        loginfo(f"Calling subconsensus for target: {sneaky_name}")
        tar_name = tar_name.lower()
        is_regex = False
        if len(tar_name) > 99:
            is_regex = True
            prefix = re.escape(tar_name[:100])
            match_name = f"(^|\s){prefix}(\d*)($|\s)"
        else:
            match_name = tar_name

        coverage = self.coverage[self.coverage['#rname'].str.lower().str.contains(
            match_name, regex=is_regex)]  # Needs to be a regex with separate match in case mismatch on length

        if coverage.shape[0] > 1:
            coverage = coverage[coverage["#rname"] == tar_name[0:100]]

        try:
            if tar_name != coverage.iloc[0][0]:
                self.naive_consensuses[tar_name] = self.naive_consensuses[coverage.iloc[0][0]]
                self.naive_consensuses.pop(coverage.iloc[0][0])
        except:
            stoperr(
                f"Couldn't match consensus {sneaky_name} to read library. Have you tried to run this on a pre-existing data folder? Are you sure your mapping reference names are compatible with Castanet?")

        if coverage.empty:
            raise loginfo(
                f"Could not generate a consensus for target: {sneaky_name}\nCheck your probe names correspond to the target names in your bam file. It's possible your probe naming scheme is incompatible with castanet.")

        if float(coverage['meanmapq'].item()) < self.a["ConsensusMapQ"]:
            '''If coverage/depth don't surpass threshold, delete grouped reads dir'''
            loginfo(
                f"Not adding subconsensus for {sneaky_name} to consensus for organism, Map Q was under minimmum threshold you set ({coverage['meanmapq'].item()})")
            return
        else:
            '''Else, call consensus on this target'''
            if not tar_name in self.subconsensuses.keys():
                self.subconsensuses[tar_name] = []
            if tar_name in self.naive_consensuses.keys():
                self.subconsensuses[tar_name].append(
                    f"{self.naive_consensuses[tar_name]}")
            elif tar_name[:100] in self.naive_consensuses.keys():
                self.subconsensuses[tar_name].append(
                    f"{self.naive_consensuses[tar_name[:100]]}")

    def collate_consensus_seqs(self, tar_name) -> None:
        '''Read and collate consensus seqs from per target to per organism'''
        try:
            seqs_and_refs = [[self.aggregate_to_probename(
                tar_name), tar_name, i] for i in self.subconsensuses[tar_name]]    # aggn, refn, seq
        except Exception as ex:
            print(f"***{ex}")
            return

        if len(seqs_and_refs) == 0:
            return

        consensus_org = seqs_and_refs[0][0]
        if not consensus_org in self.target_consensuses.keys():
            self.target_consensuses[consensus_org] = []

        for i in seqs_and_refs:
            self.target_consensuses[consensus_org].append({
                "tar_name": i[1],
                "consensus_seq":  ', '.join([i[2] for i in seqs_and_refs])
            })

    def aggregate_to_probename(self, ref) -> str:
        '''Group targets to organism via probe name (compiled in analysis.py)'''
        match = self.probe_names.iloc[np.where(
            np.isin(self.probe_names["orig_target_id"].str.lower(), ref[0:100].replace(">", "")))[0]]  # TODO < Curtailment
        if match.empty:
            print(
                f"WARNING: Couldn't match reads to probe name: {self.probe_names['target_id']}")
            return f"Unmatched"
        elif pd.notna(match["rmlst"].item()):
            return f"{match['probetype'].item()}-{match['rmlst'].item()}"
        else:
            return f"{match['probetype'].item()}"

    def call_flat_consensus(self, org_name) -> None:
        '''Create consensus sequences'''
        '''Make folder and dictionary key for supplementary stats'''
        self.eval_stats[org_name] = {}
        organism_consensus_dir = f"{self.a['folder_stem']}consensus_data/{org_name}/"
        if not os.path.isdir(organism_consensus_dir):
            shell(f"mkdir {organism_consensus_dir}")

        '''Filter bam to organism-specific targets, further filter by coverage %'''
        coverage_filter, agg_names = self.filter_bam_to_organism(org_name)

        if len(coverage_filter) == 0:
            loginfo(
                f"No remapped consensus will be generated for {org_name} as coverage was too low on all target consensues")
            self.insufficient_coverage_orgs.append(org_name)
            return

        '''Filter tar consensuses on coverage, re-make target alignment and consensus to filtered list, save'''
        self.filter_tar_consensuses(org_name, coverage_filter)
        if len(self.target_consensuses[org_name]) == 0:
            '''If all target consensuses are removed by the filter, don't make a consensus for this org'''
            loginfo(
                f"Not proceeding to generate consensus for {org_name}, as coverage across flat consensus was too low.")
            shell(f"rm -rf {organism_consensus_dir}")
            return

        # TEST -- AGGREGATE
        if agg_names.shape[0] > 1:
            # 0 Get list of targets per gene
            agg_refseq = ""
            all_subconsensuses = []
            genes = agg_names["AGGREGATE"].value_counts().index.tolist()
            genes.sort()
            for gene in genes:
                gene_targets = agg_names[agg_names["AGGREGATE"]
                                         == gene]["target_id"].tolist()
                for target in self.target_consensuses[org_name]:
                    if target["tar_name"] in gene_targets:
                        ref_names = list(set([i["tar_name"].replace(">", "").lower(
                        ) for i in self.target_consensuses[org_name] if i["tar_name"] in target["tar_name"]]))
                        tmp = target
                        tmp["refseq"] = [ref for ref in self.refs if ref[0].replace(
                            ">", "").split(" ")[0].lower() in ref_names][0][1]
                        all_subconsensuses.append(
                            [f">{target['tar_name']}", target["consensus_seq"]])

                # 1. Get each refseqs
                agg_refseq += tmp['refseq']

            # TODO << DAMN THIS IS MESSY. CLEANER TO MAKE LOTS OF FLAT CONSENSUSES THEN STICK THEM TOGETHER AFTER?

            # 2. Save concatenated refseq and cleaned subconsensuses for MSA
            # TODO < HARMONISE WITH NAMES IN USUAL PATHWAY
            ref_aln_fnme = f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_concat_refseq.fasta"
            flat_cons_seqs = f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_concat_subconsensuses.fasta"

            with open(f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_concat_refseq.fasta", "w") as f:
                [f.write(f">{gene}_concat_refseq\n{agg_refseq}\n")]
            with open(f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_concat_subconsensuses.fasta", "w") as f:
                [f.write(f"{i[0]}\n{i[1]}\n") for i in all_subconsensuses]

            # 3 make alignment of concat refseq with
            out = shell(f"mafft --thread {self.a['NThreads']} --auto --addfragments {flat_cons_seqs} {ref_aln_fnme}"
                        f"> {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln", is_test=True)
            error_handler_cli(
                out, f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln", "mafft", test_f_size=True)

            # 4 get flat consensus as per usual pathway
            flat_consensus = self.dumb_consensus_AGGREGATE(
                f"{self.a['folder_stem']}consensus_data/{org_name}/", org_name)

        else:
            self.build_msa_requisites(org_name)

            # TODO << harmonie flatten_consensus input to work with aggregate mode too
            flat_consensus = self.flatten_consensus(org_name)

        save_fa(f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_flat_consensus_sequence.fasta",
                f">{org_name}_consensus\n{flat_consensus}")

        '''Remap to re-made flat consensus, to make `re-mapped consensus`'''
        self.remap_flat_consensus(org_name)

        '''Dump any additional stats to pickle'''
        self.dump_stats(org_name)

    def dump_stats(self, org_name) -> None:
        with open(f"{self.a['folder_stem']}consensus_data/{org_name}/supplementary_stats.p", 'wb') as f:
            p.dump(self.eval_stats[org_name], f, protocol=p.HIGHEST_PROTOCOL)

    def build_msa_requisites(self, org_name) -> None:
        '''Create fasta files containing target reference seqs and consensus seqs, for downstream MSA'''
        ref_seq_names = list(set([i["tar_name"].replace(">", "").lower()[:100]
                                  for i in self.target_consensuses[org_name]]))
        ref_seqs = [ref for ref in self.refs if ref[0].replace(
            ">", "").split(" ")[0].lower() in ref_seq_names]

        assert len(
            ref_seqs) > 0, f"Couldn't match ref sequences to target name for {org_name}"

        with open(self.fnames['flat_cons_refs'], "w") as f:
            [f.write(f"{i[0]}\n{i[1]}\n") for i in ref_seqs]
        with open(self.fnames['flat_cons_seqs'], "w") as f:
            [f.write(f">{i['tar_name']}_CONS\n{i['consensus_seq']}\n")
             for i in self.target_consensuses[org_name]]
        shell(
            f"cat {self.fnames['flat_cons_seqs']} {self.fnames['flat_cons_refs']} > {self.a['folder_stem']}/consensus_data/unaligned_consensuses_and_refs.fna")

    def flatten_consensus(self, org_name) -> str:
        '''Make MSA of references, then add fragments from target consensuses'''
        loginfo(f"making consensus alignments for target group: {org_name}")
        ref_aln_fname = f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_ref_alignment.aln"
        ref_aln = read_fa(self.fnames['flat_cons_refs'])
        assert len(
            ref_aln) > 0, f"Reference alignment for {org_name} doesn't exist."

        if len(ref_aln) > 1:
            '''Align flat consensus references'''
            out = shell(
                f"mafft --thread {self.a['NThreads']} --auto {self.fnames['flat_cons_refs']} > {ref_aln_fname}", is_test=True)

            error_handler_cli(out, ref_aln_fname, "mafft")
        else:
            '''If only 1 reference seq, the alignment wouldn't have worked - defer to temp refs file in these cases'''
            ref_aln_fname = self.fnames['flat_cons_refs']

        ref_aln_with_reads_fname = f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln"
        out = shell(f"mafft --thread {self.a['NThreads']} --auto --addfragments {self.fnames['flat_cons_seqs']} {ref_aln_fname}"
                    f"> {ref_aln_with_reads_fname}", is_test=True)

        error_handler_cli(out, ref_aln_with_reads_fname,
                          "mafft", test_f_size=True)

        try:
            if self.a["DebugMode"]:
                call_graph(self.a["ExpName"], org_name, f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_alignment.aln",
                           f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_target_consensus_alignment", self.a["SaveDir"], is_eval=False)
        except FileNotFoundError:
            raise SystemError(
                "Castanet couldn't construct a consensus alignment graph")

        '''Return flat consensus'''
        return self.dumb_consensus(f"{self.a['folder_stem']}consensus_data/{org_name}/", org_name)

    def dumb_consensus_AGGREGATE(self, alnfpath, org_name) -> list:
        '''Produce an un-referenced/`flat` consensus sequence for file of target and target ref seqs'''
        # TODO < Disgusting, harmonise with regular dumb_consensus.
        # (15/04/26) Bodged like this as concat alignments are necessarily sparse. Unclear whether to concat smaller subalignments or keep as-is.
        from Bio.Align.AlignInfo import SummaryInfo
        aln = AlignIO.read(
            f"{alnfpath}{org_name}_consensus_alignment.aln", 'fasta')
        return str(SummaryInfo(aln).dumb_consensus())

    def dumb_consensus(self, alnfpath, org_name) -> list:
        '''Produce an un-referenced/`flat` consensus sequence for file of target and target ref seqs'''
        def base_cons(s):
            ''' Return strict consensus for a set of bases (eg. column in alignment), ignoring gaps. '''
            if not s or (len(s.replace("-", "")) <= 0.1 * len(s)):
                return ('', np.nan)

            try:
                just_measured_bases = s[-int(len(s))                                        :].lower().replace("-", "").replace("n", "")
                consbase, consnum = Counter(
                    just_measured_bases).most_common()[0]

            except IndexError:  # TODO < MESSY
                return ('', np.nan)

            if float(consnum)/len(s) < 0.1:
                consbase, consnum = Counter(
                    s.replace("-", "")).most_common()[0]

            return consbase, float(consnum)/len(s)

        aln = AlignIO.read(
            f"{alnfpath}{org_name}_consensus_alignment.aln", 'fasta')
        cluster_cons = pd.DataFrame(
            pd.Series(base_cons(aln[:, i])) for i in range(len(aln[0]))).dropna()

        '''Plot identity for QC'''
        cluster_cons.columns = ['cons', 'ident']
        cluster_cons["ident"].rolling(120).mean().plot()
        if self.a["DebugMode"]:
            fig = px.line(x=cluster_cons.index,
                          y=cluster_cons["ident"], title="Flat consensus identity")
            fig.write_image(
                f"{alnfpath}{org_name}_flat_consensus_identity.png")
        return "".join(cluster_cons["cons"].tolist())

    def filter_bam_to_organism(self, org_name) -> list:
        '''Output coverage stats for target consensuses'''
        if self.probe_names["probetype"].eq(org_name).any():
            probels = self.probe_names[self.probe_names['probetype']
                                       == org_name]['orig_target_id'].str.lower().tolist()
        else:
            probels = self.probe_names[self.probe_names['genename']
                                       == org_name]['orig_target_id'].str.lower().tolist()
        coverage_df = self.coverage[
            self.coverage["#rname"].isin(probels) | self.coverage["#rname"].str.startswith(tuple(probels), na=False)]
        assert not coverage_df.empty, f"Call to samtools coverage returned empty output. Check that your bam file is indexed and that the path to it is correct."

        '''Get coverage for each consensus, filter collated bam by consensus coverage and map q'''
        coverage_df = coverage_df[(coverage_df["coverage"] >= self.a['ConsensusCoverage']) & (
            coverage_df["meanmapq"] >= self.a['ConsensusMapQ'])]
        coverage_df.to_csv(
            f"{self.a['folder_stem']}consensus_data/{org_name}/target_consensus_coverage.csv", index=False, header=False)
        coverage_filter = coverage_df["#rname"].tolist()

        if len(coverage_filter) == 0:
            return [], pd.DataFrame()
        else:
            '''Filter master bam for targets with sufficient coverage'''
            probels = [f'{i}' for i in coverage_filter]
            # Split into 2 due to backslash in f string
            probels = [i.replace("|", "\|") for i in probels]
            cmd = f"samtools view -@ {self.a['NThreads']} -b {self.fnames['master_bam']} {' '.join(probels)} > {self.a['folder_stem']}consensus_data/{org_name}/collated_reads.bam"
            shell(cmd)
            # TODO test output

            '''Estimate number of mapped reads in the final alignment (get just primary mapped reads, div 2 to average F & R strands)'''
            self.eval_stats[org_name]["filtered_collated_read_num"] = round(samtools_read_num(
                f"{self.a['folder_stem']}consensus_data/{org_name}/collated_reads.bam", f'-F 0x904 -q {self.a["ConsensusMapQ"]}') / 2)

            ##### TEST -- bact aggregates ####
            agg_names = self.probe_names[self.probe_names["target_id"].isin(
                coverage_filter)]
            ########

            return coverage_filter, agg_names

    def filter_tar_consensuses(self, org_name, filter) -> None:
        '''Purge target consensus from master list if coverage was lower than threshold (aln is consequently remade)'''
        to_del = [i for i in range(len(self.target_consensuses[org_name]))
                  if not self.target_consensuses[org_name][i]["tar_name"].replace(">", "") in filter]
        if not to_del:
            return
        else:
            to_del.reverse()
            for i in to_del:
                '''Not done in loop and in reverse to not break the iterator'''
                del self.target_consensuses[org_name][i]
            return

    # @timing
    def remap_flat_consensus(self, org_name) -> None:
        '''Remap reads to flattened consensus, save, call stats, remove raw fastas'''
        flat_cons_fname = f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_remapped_consensus_sequence.fasta"

        if self.a["Mapper"] == "bwa":
            bwa_index(
                f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_flat_consensus_sequence.fasta")
            shell(f"samtools fastq -@ {self.a['NThreads']} {self.fnames['master_bam']} |"
                  f"bwa-mem2 mem -t {self.a['NThreads']} {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_flat_consensus_sequence.fasta - | "
                  f"viral_consensus -i - -r {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_flat_consensus_sequence.fasta -o {flat_cons_fname} --min_depth {self.a['ConsensusMinD']} --out_pos_counts {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_pos_counts.csv")

        # TODO < tidy common elements # TODO < delete btl intermediate files
        elif self.a["Mapper"] == "bowtie2":
            shell(
                f"bowtie2-build {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_flat_consensus_sequence.fasta {self.a['folder_stem']}consensus_data/{org_name}/reference_indices", is_test=True)
            shell(f"samtools fastq -@ {self.a['NThreads']} {self.fnames['master_bam']} |"
                  f"bowtie2 -x {self.a['folder_stem']}consensus_data/{org_name}/reference_indices -U - -p {self.a['NThreads']} --local -I 50 --maxins 2000 --no-unal |"
                  f"viral_consensus -i - -r {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_flat_consensus_sequence.fasta -o {flat_cons_fname} --min_depth {self.a['ConsensusMinD']} --out_pos_counts {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_pos_counts.csv")

        elif self.a["Mapper"] == "minimap2":
            shell(f"samtools fastq -@ {self.a['NThreads']} {self.fnames['master_bam']} |"
                  f"minimap2 -ax map-ont {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_flat_consensus_sequence.fasta - |"
                  f"viral_consensus -i - -r {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_flat_consensus_sequence.fasta -o {flat_cons_fname} --min_depth {self.a['ConsensusMinD']} --out_pos_counts {self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_pos_counts.csv")

        try:
            error_handler_cli("", flat_cons_fname,
                              "viral_consensus", test_f_size=True)
        except Exception as ex:
            logerr(f"Consnesus sequence generation failed for {org_name} while attempting to remap to the flat consensus. "
                   f"Skipping remapped consensus generation. Error details: {ex}")
            return

        try:
            self.fix_terminal_gaps(f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_consensus_pos_counts.csv",
                                   f"{self.a['folder_stem']}consensus_data/{org_name}/{org_name}_remapped_consensus_sequence.fasta")
        except FileNotFoundError as ex:
            stoperr(f"Castanet call to ViralConsensus produced empty output. Check that it's installed using the dependency_check endpoint (see readme)")
        find_and_delete(
            f"{self.a['folder_stem']}consensus_data/{org_name}/", f"*.fasta.*")

    def fix_terminal_gaps(self, in_fname, out_fname) -> None:
        '''Trim terminal gaps'''
        cons = pd.read_csv(in_fname, sep="\t")
        n_pos = cons.shape[0]
        if self.a["ConsensusTrimTerminals"]:
            '''If total reads at pos x < threshold AND in leading/trailing 5% of reads, mark for deletion'''
            cons["del"] = cons.apply(lambda x: np.where(x["Total"] < self.a["ConsensusMinD"] and (
                x["Pos"] < n_pos * 0.05 or x["Pos"] > n_pos * 0.95), 1, 0), axis=1)
        else:
            cons["del"] = 0
        '''Rm terminal gaps, re-index, re-call consensus.'''
        cons = cons[cons["del"] == 0]
        cons = cons.drop(columns=["Pos", "Total", "del"])
        cons["con"] = cons.apply(lambda x: x.idxmax(), axis=1)
        '''Re-do index and totals'''
        cons["Pos"] = np.arange(1, cons.shape[0] + 1)
        cons["Total"] = cons.apply(
            lambda x: x["A"] + x["T"] + x["C"] + x["G"] + x["-"], axis=1)
        '''Fix artificial adnylation where totals = 0 (pd idxmax annoyingly picks first col in this case)'''
        cons["con"] = cons.apply(
            lambda x: x["con"] if not x["Total"] == 0 else "-", axis=1)
        cons["-"] = cons.apply(lambda x: 1 if x["Total"] == 0 else 0, axis=1)
        cons["con"] = cons["con"].astype(str)
        cons.to_csv(in_fname)
        # RM < TODO Add in deduplicated depth :S
        fasta_header = f"{self.a['ExpName']}_{in_fname.split('/')[-1].split('_')[0]}_consensus_MinDepth{self.a['ConsensusMinD']}"
        save_fa(
            out_fname, f">{fasta_header}\n{''.join(cons['con'].tolist())}")
        save_fa(
            f"{self.a['folder_stem']}/consensus_sequences/{out_fname.split('/')[-1]}", f">{fasta_header}\n{''.join(cons['con'].tolist())}")
        loginfo(
            f"INFO: Consensus sequence saved to {self.a['folder_stem']}/consensus_sequences/")

    def clean_incomplete_consensus(self) -> None:
        '''If we had insufficient coverage for organism x, clean it from self vars and folder tree'''
        self.subconsensuses.clear()
        for org_name in self.insufficient_coverage_orgs:
            del self.target_consensuses[org_name]
            shell(f"rm -r {self.a['folder_stem']}/consensus_data/{org_name}/")
        rm(f"{self.fnames['flat_cons_seqs']} {self.fnames['flat_cons_refs']}")
        rm(f"{self.a['folder_stem']}/consensus_data/unaligned_consensuses_and_refs.fna")

    def tidy(self) -> None:
        '''Remove intermediate files to save disc space'''
        rm(f"{self.fnames['collated_reads_fastq']}")
        rm(f"{self.fnames['temp_folder']}", "-r")
        if not self.a["DebugMode"]:
            rm(f"{self.fnames['grouped_reads']}")

    def generate_summary(self, org) -> None:
        try:
            dfpath = f"{self.a['folder_stem']}/consensus_seq_stats.csv"
            cols = ["target", "n_bases", "n_reads",
                    "gc_pc", "missing_bs", "ambig_bs", "cov"]
            df = pd.DataFrame(columns=cols)

            if os.path.isfile(dfpath):
                os.remove(dfpath)

            '''remapped cons stats'''
            c_df = pd.read_csv(
                f"{self.a['folder_stem']}/consensus_data/{org}/{org}_consensus_pos_counts.csv")
            gc = round((c_df["G"].sum() + c_df["C"].sum()) /
                       c_df["Total"].sum() * 100, 2)

            missing = c_df[c_df["Total"] == 0].shape[0]
            ambigs = c_df[(c_df["-"] != 0) & (c_df["A"] == 0) & (c_df["C"]
                                                                 == 0) & (c_df["T"] == 0) & (c_df["G"] == 0)].shape[0]
            coverage = round(
                (1 - ((missing + ambigs) / c_df["Total"].sum())) * 100, 2)

            '''Get additional stats on consensus remapping'''
            try:
                additional_stats = {}
                with open(f"{self.a['folder_stem']}/consensus_data/{org}/supplementary_stats.p", 'rb') as f:
                    additional_stats["n_remapped_seqs"] = p.load(
                        f)["filtered_collated_read_num"]

                c_stats = pd.DataFrame([[org, c_df["Total"].sum(
                ), additional_stats['n_remapped_seqs'], gc, missing, ambigs, coverage]], columns=cols)
                df = pd.concat([df, c_stats], axis=0, ignore_index=True)
                df.to_csv(dfpath)
            except FileNotFoundError:
                logerr(f"Couldn't find supplementary stats for {org}, skipping addition of remapped read count to summary csv."
                       f"This can happen if individual pipeline stages are run out-of-synch with each other, or if you're generating a consensus with horizontal (rMLST) aggregation.")

            if self.a["DebugMode"]:
                '''Plot consensus coverage'''
                px.line(c_df, x="Pos", y="Total", title=f"Consensus coverage, {org} ({self.a['ExpName']})",
                        labels={"Pos": "Position", "Total": "Num Reads"}).write_image(
                    f"{self.a['folder_stem']}/consensus_data/{org}/{org}_consensus_coverage.png")
        except Exception as ex:
            logerr(
                f"Couldn't generate summary for {org}. This usually happens if a consensus sequence failed to generate. Error details: {ex}")

    def main(self) -> None:
        '''Entrypoint. Index main bam, filter it, make target consensuses, then create flattened consensus'''
        end_sec_print(
            "INFO: Calling consensus sequences\nThis may take a little while...")
        samtools_index(f"{self.fnames['master_bam']}")
        group_consensus_fname = f"{self.a['folder_stem']}consensus_temp.fasta"

        '''Get consensus and coverage for each target, memoize'''
        out = shell("samtools", is_test=True)
        shell(  # Quicker all in one than out of loop
            f"""samtools consensus -@ {self.a['NThreads']} --min-depth {self.a["ConsensusMinD"]} -f fasta '{self.fnames['master_bam']}' -o '{group_consensus_fname}'""")
        error_handler_cli(out, group_consensus_fname,
                          "samtools", test_out_f=True, test_f_size=True)

        naive_consensuses_raw = {i[0].replace(
            ">", ""): i[1] for i in read_fa(f"{group_consensus_fname}")}
        for key in naive_consensuses_raw.keys():
            if len(key) > 100:
                self.naive_consensuses[key[0:100].lower(
                )] = naive_consensuses_raw[key]
            else:
                self.naive_consensuses[key.lower(
                )] = naive_consensuses_raw[key]

        end_sec_print("INFO: Calling coverage across all targets")
        self.coverage = pd.read_csv(io.StringIO(shell(f"samtools coverage '{self.fnames['master_bam']}'", "Coverage, consensus filter bam", ret_output=True).decode(
        )), sep="\t")
        assert not self.coverage.empty, "Call to samtools coverage returned empty output. Check that your bam file is indexed and that the path to it is correct."
        self.coverage["#rname"] = self.coverage["#rname"].str.lower()

        for key in self.grouped_reads.keys():
            self.filter_bam(key)
        self.grouped_reads.clear()

        '''Consensus for each thing target group'''
        [self.collate_consensus_seqs(tar_name)
            for tar_name in self.subconsensuses.keys()]
        [self.call_flat_consensus(
            i) for i in self.target_consensuses.keys() if i != "Unmatched"]

        '''Tidy up'''
        self.clean_incomplete_consensus()
        self.tidy()

        '''Call CSV summary generator'''
        [self.generate_summary(i) for i in os.listdir(
            f"{self.a['folder_stem']}/consensus_data/") if not "GROUND_TRUTH" in i and not ".fna" in i and not i.startswith(".")]

        '''Final tidy up'''
        shell(f"rm {group_consensus_fname}")
        find_and_delete(
            f"{self.a['folder_stem']}consensus_data/", "*.p")
        find_and_delete(
            f"{self.a['folder_stem']}consensus_data/", "*.bam")

        end_sec_print("INFO: Consensus calling complete")


if __name__ == "__main__":
    cls = Consensus()
    cls.main()
