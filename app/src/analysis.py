from __future__ import division
import os
import re
import numpy as np
import pandas as pd
import plotly.express as px

from app.utils.system_messages import end_sec_print
from app.utils.shell_cmds import loginfo, stoperr, logerr, shell
from app.utils.error_handlers import error_handler_analysis
from app.utils.basic_cli_calls import get_read_num
from app.utils.utility_fns import trim_long_fpaths, read_fa, enumerate_bam_files


class Analysis:
    def __init__(self, argies, start_with_bam, api_entry=True) -> None:
        self.a = argies
        self.output_dir = f"{self.a['SaveDir']}/{self.a['ExpName']}/"
        self.bam_fname = f"{self.a['SaveDir']}/{self.a['ExpName']}/{self.a['ExpName']}.bam" if not start_with_bam else f"{argies['ExpDir']}/{[i for i in os.listdir(argies['ExpDir']) if i[-4:] == '.bam'][0]}"
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        if not os.path.exists(self.bam_fname):
            '''If entry from analyse endpoint, cp bam file from input folder to experiment folder'''
            shell(
                f"cp {enumerate_bam_files(self.a['ExpDir'])} {self.bam_fname}")
        if api_entry:
            self.a["input_file"] = f"{self.output_dir}/{self.a['ExpName']}_PosCounts.csv"
        self.df = error_handler_analysis(self.a)
        self.df["target_id"] = self.df["target_id"].str.lower()
        self.lut = pd.read_csv(self.a["MappingRefTable"], index_col=False)
        self.probe_regexes = [
            re.compile(r'bact[0-9]+\_[\s\S]*'),
            # re.compile(r'bact[0-9]+_([A-Za-z]+)-[0-9]+[-_]([A-Za-z]+)'), # TODO < DEPRECATED AS OF 9.3
            # re.compile(r'bact[0-9]+_[0-9]+_([A-Za-z]+_[A-Za-z_]+)'),
            # re.compile(r'bact[0-9]+_([a-z]+_[a-z_]+)'),
            # re.compile(r'bact[0-9]+_([A-Za-z]+)-[0-9]+')
        ]

    def add_probelength(self):
        '''Add length of target_id to each row of master df after splitting probelength data.'''
        loginfo(f"Generating probe lengths from input probes file (RefStem)")
        try:
            plens = [{"target_id": i[0].replace(
                ">", ""), "target_len": len(i[1])} for i in read_fa(self.a["RefStem"])]
            probelengths = pd.DataFrame(plens)
            probelengths = probelengths.sort_values(by="target_id")
            probelengths.to_csv(
                f"{self.output_dir}/probe_lengths.csv", index=False)
        except:
            stoperr(
                f'Failed to read probe information. Is {self.a["RefStem"]} a valid multifasta file?')

        probelengths_mod = self.add_probetype(probelengths)
        self.df = self.df.merge(
            probelengths_mod, left_on='target_id', right_on='target_id', how='left')
        return probelengths_mod

    def add_probetype(self, pdf):
        ''' For probe aggregation, determine organism and gene ID.
        Splits genename and probetype into separate cols, then does manual adjustments
        '''
        loginfo('Aggregating by organism and gene name.')

        '''Apply normalisation to both probe and master dataframes to allow for different probe name conventions'''
        pdf['orig_target_id'] = pdf['target_id'].copy()
        pdf['orig_target_id'] = pdf.apply(
            lambda x: trim_long_fpaths(x["orig_target_id"]), axis=1)
        pdf['target_id'] = pdf['target_id'].str.lower()

        pdf["key"] = pdf["target_id"].str.split("_").str[-1]
        self.lut["key"] = self.lut["key"].astype(str)

        pdf["organism"] = pdf.apply(
            lambda x: self.lut[self.lut["key"] == x["key"]]["organism"].iloc[0], axis=1)
        pdf["rmlst"] = pdf.apply(lambda x: self.lut[self.lut["key"] == x["target_id"].split(
            "_")[-1]]["rmlst"].iloc[0], axis=1)
        pdf['genename'] = pdf.target_id.str.lower().apply(
            lambda x: x.split("_")[0])

        pdf["genename"] = pdf.apply(lambda x: x["genename"].lower(), axis=1)
        '''More precise definition for the different virus types'''
        pdf.loc[pdf.target_id == 'roseolovirus_allrecords_cluster_1',
                'genename'] = 'HHV7_roseolovirus_allrecords_cluster_1'
        pdf.loc[pdf.target_id == 'roseolovirus_allrecords_cluster_2',
                'genename'] = 'HHV6_roseolovirus_allrecords_cluster_2'
        pdf.loc[pdf.genename == 'enterovirus', 'genename'] = pdf.loc[pdf.genename ==
                                                                     'enterovirus'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))
        pdf.loc[pdf.genename == 'coronaviridae', 'genename'] = pdf.loc[pdf.genename ==
                                                                       'coronaviridae'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:3]))
        pdf.loc[pdf.genename == 'adenoviridae', 'genename'] = pdf.loc[pdf.genename ==
                                                                      'adenoviridae'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))
        pdf.loc[pdf.genename == 'flaviviridae', 'genename'] = pdf.loc[pdf.genename ==
                                                                      'flaviviridae'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))
        pdf.loc[pdf.genename == 'influenza', 'genename'] = pdf.loc[pdf.genename ==
                                                                   'influenza'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))
        pdf.loc[pdf.genename == 'paramyxoviridae', 'genename'] = pdf.loc[pdf.genename ==
                                                                         'paramyxoviridae'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))
        pdf.loc[pdf.genename == 'parvoviridae', 'genename'] = pdf.loc[pdf.genename ==
                                                                      'parvoviridae'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))

        '''Append and apply horizontal aggregation keys (i.e. rmlst)'''
        pdf['probetype'] = pdf.genename.str.lower()
        pdf["genename"] = pdf.apply(lambda x: x["organism"] + "-" + x["rmlst"] if pd.notna(x["rmlst"]) else x["organism"], axis=1)
        pdf["AGGREGATE"] = pdf.apply(lambda x: f'{x["organism"]}' if str(
            x["rmlst"]).startswith("rmlst") else x["target_id"].split("_")[0], axis=1)


        loginfo(
            f'Organism and gene summary: {pdf.organism.nunique()} organisms, up to {pdf.groupby("probetype").probetype.nunique().max()} aggregation levels (probetype) each and up to {pdf.groupby("probetype").genename.nunique().max()} genes each.')
        pdf.to_csv(f"{self.output_dir}/probe_aggregation.csv")

        if pdf[pdf["probetype"] == ""].shape[0] > 0:
            logerr(
                f"Failure decoding the name of one or more probe types: \n {pdf[pdf['probetype'] == '']} \n Please check your probe naming conventions are compatible with Castanet")

        return pdf

    def add_depth(self, probelengths):
        ''' Calculate read depth per position. '''
        loginfo('Calculating read depth information.')
        depth = None
        metrics = {}
        odir = f'{self.output_dir}/Depth_output'
        raw_readcount = get_read_num(self.a, self.bam_fname)
        if os.path.isdir(odir):
            '''Clear dir if already exists'''
            shell(f"rm -r {odir}")
        try:
            os.mkdir(odir)
        except OSError:
            odir = os.getcwd()
            logerr(
                f'Cannot create output directory {odir} for saving depth information. Proceeding with current working directory.')

        loginfo(
            'INFO: Calculating read depth statistics for all probes, for all samples.')
        for (sampleid, probetype, organism), g in self.df.groupby(['sampleid', 'probetype', 'organism']):
            orig_probetype = probetype
            gene_list = g.genename.unique()
            n_genes = len(gene_list)
            target_list = g.target_id.unique()
            n_targets = len(target_list)
            Dd, D1d = {}, {}
            for genename, gg in g.groupby('genename'):
                loginfo(f'Processing {sampleid} - {genename}')
                '''Generate two arrays for each probetype-genename group (D = number of occurrences, D1 = unique start/end positions)'''
                D = np.zeros(int(gg.target_len.max()), dtype=np.uint32)
                D1 = np.zeros(D.shape, dtype=np.uint32)
                for target_id, gt in gg.groupby('target_id'):
                    try:
                        '''Don't give the user the hashed header name, it will only upset them'''
                        sneaky_name = f'{target_id.split("_")[0]}_{self.lut[self.lut["key"] == target_id.split("_")[-1]]["description"].item()[0:100]}'
                    except TypeError:
                        '''If user has somehow broken fasta header'''
                        sneaky_name = f'{target_id.split("_")[0]}'
                    loginfo(f'..... target: {sneaky_name}.')
                    for _, row in gt.iterrows():
                        D[row.startpos-1:row.startpos-1+row.maplen] += row.n
                        D1[row.startpos-1:row.startpos-1+row.maplen] += 1
                Dd[genename] = D
                D1d[genename] = D1

            '''Max possible number of targets for this probetype, aggregating all genes'''
            nmax_targets = probelengths[probelengths.probetype ==
                                        probetype].target_id.nunique()
            '''Max possible number of genes for this probetype, aggregating all genes'''
            nmax_genes = probelengths[probelengths.probetype ==
                                      probetype].genename.nunique()
            '''Max possible positions for this probetype, aggregating all genes'''
            nmax_probetype = probelengths[probelengths.probetype == probetype].groupby(
                'genename').target_len.max().sum()

            '''Collapse to a single array for the entire genename group for this probetype in this sample'''
            D = np.hstack(list(Dd.values()))
            D1 = np.hstack(list(D1d.values()))
            valid_mask = (D1 != 0) & ~np.isnan(D) & ~np.isnan(D1)
            amprate = (D[valid_mask] / D1[valid_mask])
            '''Max possible positions for the genes that were actually in this BAM (accounts for some genes not being captured)'''
            npos = len(D)
            '''Now pad out with zeros to the total number of mappable positions for this probetype (nmax_probetype above)'''
            D = np.pad(D, (0, nmax_probetype - npos),
                       'constant', constant_values=0)
            D1 = np.pad(D1, (0, nmax_probetype - npos),
                        'constant', constant_values=0)

            loginfo(
                f'Mean depth (all reads) for {orig_probetype}: {D.mean()}')
            loginfo(
                f'Mean depth (deduplicated) for {orig_probetype}: {D1.mean()}')
            '''Amplification rate calculations use the unpadded (mapped) number of sites as the denominator'''
            loginfo(
                f'Mean amplification ratio for {orig_probetype}: {amprate.mean()}')

            if self.a["DebugMode"]:
                '''Save arrays as CSV'''
                with open(f'{odir}/{orig_probetype}-{sampleid}_depth_by_pos.csv', 'a') as o:
                    np.savetxt(o, D, fmt='%d', newline=',')
                    o.write('\n')
                    np.savetxt(o, D1, fmt='%d', newline=',')
                    o.write('\n')

            '''Save array plots as pdf if significant'''
            if D1.mean() >= 0.01:
                plot_df = pd.DataFrame()
                plot_df["position"], plot_df["All Reads"], plot_df["Deduplicated Reads"] = np.arange(
                    0, D.shape[0]), D, D1
                fig = px.line(plot_df, x="position", y=[
                    "All Reads", "Deduplicated Reads"], title=f'{sampleid}\n{orig_probetype} ({n_targets}/{nmax_targets} targets in {n_genes}/{nmax_genes} genes)',
                    labels={"position": "Position", "value": "Num Reads"})
                fig.update_layout(legend={"title_text": "", "orientation": "h", "entrywidth": 100,
                                          "yanchor": "bottom", "y": 1.02, "xanchor": "right", "x": 1})
                fig.write_image(f'{odir}/{orig_probetype}-{sampleid}.png')

            '''Build up dictionary of depth metrics for this sample and probetype'''
            metrics[sampleid, probetype, organism] = (g.n.sum(), g.n.count(), n_targets, n_genes, nmax_targets, nmax_genes, nmax_probetype, npos,
                                                      amprate.mean(), amprate.std(), np.median(amprate),
                                                      D.mean(), D.std(), np.percentile(D, 25), np.median(D), np.percentile(D, 75),
                                                      raw_readcount,
                                                      (D > 0).sum(),
                                                      (D >= 2).sum(),
                                                      (D >= 5).sum(),
                                                      (D >= 10).sum(),
                                                      (D >= 100).sum(),
                                                      (D >= 1000).sum(),
                                                      D1.mean(), D1.std(), np.percentile(D1, 25), np.median(D), np.percentile(D1, 75),
                                                      (D1 > 0).sum(),
                                                      (D1 >= 2).sum(),
                                                      (D1 >= 5).sum(),
                                                      (D1 >= 10).sum(),
                                                      (D1 >= 100).sum(),
                                                      (D1 >= 1000).sum())

        '''Data frame of all depth metrics'''
        depth = pd.DataFrame(metrics, index=['reads_for_mapping', 'n_reads_dedup', 'n_targets', 'n_genes', 'nmax_targets', 'nmax_genes', 'npos_max_probetype', 'npos_cov_probetype',
                                             'amprate_mean', 'amprate_std', 'amprate_median',
                                             'depth_mean', 'depth_std', 'depth_25pc', 'depth_median', 'depth_75pc',
                                             'raw_readcount',
                                             'npos_cov_mindepth1',
                                             'npos_cov_mindepth2',
                                             'npos_cov_mindepth5',
                                             'npos_cov_mindepth10',
                                             'npos_cov_mindepth100',
                                             'npos_cov_mindepth1000',
                                             'udepth_mean', 'udepth_std', 'udepth_25pc', 'udepth_median', 'udepth_75pc',
                                             'npos_dedup_cov_mindepth1',
                                             'npos_dedup_cov_mindepth2',
                                             'npos_dedup_cov_mindepth5',
                                             'npos_dedup_cov_mindepth10',
                                             'npos_dedup_cov_mindepth100',
                                             'npos_dedup_cov_mindepth1000']).T.reset_index()
        depth.rename(columns={'level_0': 'sampleid',
                              'level_1': 'AGGREGATE',
                              'level_2': "organism"}, inplace=True)

        '''Add reads on target (rot)'''
        try:
            rot = depth.groupby(
                'sampleid').reads_for_mapping.sum().reset_index()
        except:
            stoperr(
                f'ERROR: Depth dataframe seems to be empty. This usually happens when the naming conventions in your mapping reference file can\'t be resolved by Castanet: please see readme for details of correct formatting.')

        rot.rename(
            columns={'reads_for_mapping': 'reads_on_target'}, inplace=True)
        depth = depth.merge(rot, on='sampleid', how='left')
        rot_dedup = depth.groupby(
            'sampleid').n_reads_dedup.sum().reset_index()
        rot_dedup.rename(
            columns={'n_reads_dedup': 'reads_on_target_dedup'}, inplace=True)
        depth = depth.merge(rot_dedup, on='sampleid', how='left')
        depth['prop_of_reads_on_target'] = depth.reads_for_mapping / \
            depth.reads_on_target
        depth['prop_npos_cov1'] = depth.npos_cov_mindepth1 / \
            depth.npos_max_probetype
        depth['prop_npos_cov2'] = depth.npos_cov_mindepth2 / \
            depth.npos_max_probetype
        depth['prop_npos_cov5'] = depth.npos_cov_mindepth5 / \
            depth.npos_max_probetype
        depth['prop_npos_cov10'] = depth.npos_cov_mindepth10 / \
            depth.npos_max_probetype
        depth['prop_npos_cov100'] = depth.npos_cov_mindepth100 / \
            depth.npos_max_probetype
        depth['prop_npos_cov1000'] = depth.npos_cov_mindepth1000 / \
            depth.npos_max_probetype
        depth['prop_ntargets'] = depth.n_targets/depth.nmax_targets
        depth['prop_ngenes'] = depth.n_genes/depth.nmax_genes
        '''Add log transforms'''
        depth['log10_depthmean'] = depth.depth_mean.apply(
            lambda x: np.log10(x+1))
        depth['log10_udepthmean'] = depth.udepth_mean.apply(
            lambda x: np.log10(x+1))
        '''Duplicated reads only'''
        depth['clean_reads_for_mapping'] = depth.reads_for_mapping - \
            depth.n_reads_dedup
        depth['clean_prop_of_reads_on_target'] = (
            depth.reads_for_mapping-depth.n_reads_dedup)/depth.reads_on_target
        '''Save and log'''
        loginfo(f'Saving {self.output_dir}/{self.a["ExpName"]}_depth.csv.')
        depth.to_csv(
            f'{self.output_dir}/{self.a["ExpName"]}_depth.csv', index=False)
        loginfo(
            f'Mean read depth per sample: \n{depth.groupby("sampleid").depth_mean.mean().to_string()}')
        return depth

    def add_read_d_and_clin(self, depth):
        ''' Add raw read numbers and any external categorical/clinical data.
        If specified, samples file must supply at least the following columns: {}.
        If not specified, infer raw read num from input bam (assumes no prior filtering!!)'''
        loginfo('Adding sample information and clinical data.')
        read_num = get_read_num(self.a, self.bam_fname)
        samples = pd.DataFrame(
            [{"sampleid": str(self.a["ExpName"]), "pt": "", "rawreadnum": read_num}])

        depth['sampleid'] = depth['sampleid'].astype(str)
        samples['sampleid'] = samples['sampleid'].astype(str)
        '''Merge read n (and clin data if supplied) to depth counts, return'''
        cdf = depth.merge(samples, on='sampleid', how='left')
        cdf['readprop'] = cdf.reads_for_mapping/cdf.rawreadnum
        loginfo(
            f'Added the following columns to depth csv: {list(samples.columns)}')
        cdf.to_csv(
            f'{self.output_dir}/{self.a["ExpName"]}_depth.csv', index=False)
        return cdf

    def save_tophits(self, depth):
        ''' Simple output of likely best hit for each sample.'''
        loginfo('Target with highest proportion of captured reads:')
        tophits = depth.sort_values(['sampleid', 'clean_prop_of_reads_on_target'],
                                    ascending=True).drop_duplicates('sampleid', keep='last')
        loginfo(tophits[['sampleid', 'probetype',
                'clean_prop_of_reads_on_target']])
        tophits.to_csv(
            f'{self.output_dir}{self.a["ExpName"]}_tophit.csv', index=False)
        loginfo(
            f'Saved top hits to {self.output_dir}{self.a["ExpName"]}_tophits.csv')

    def read_dist_piechart(self):
        df = pd.read_csv(
            f"{self.output_dir}{self.a['ExpName']}_depth.csv")
        fig = px.pie(df, values=df["reads_for_mapping"], names=df["probetype"],
                     title=f"Read distribution, {self.a['ExpName']}")
        fig.update_traces(textposition='inside',
                          textinfo='percent+label+value')
        fig.write_image(
            f"{self.output_dir}/{self.a['ExpName']}_read_distributions.png")

    def get_cov(self, row, all_cov):
        cov = row["npos_cov_mindepth2"] / row["npos_max_probetype"]
        if not row["sampleid"] in all_cov.keys():
            all_cov[row["sampleid"]] = {}
        all_cov[row["sampleid"]][row["probetype"]] = round(cov, 2)

    def read_coverage_chart(self):
        df = pd.read_csv(
            f"{self.output_dir}{self.a['ExpName']}_depth.csv")

        df = df.rename(columns={"AGGREGATE": "probetype"})
        all_cov = {}
        df.apply(lambda x: self.get_cov(x, all_cov), axis=1)
        df_cov = pd.DataFrame(all_cov).T
        df_cov = df_cov.reindex(sorted(df_cov.columns), axis=1)
        df_cov.to_csv(f"{self.output_dir}{self.a['ExpName']}_coverage.csv")
        df.to_csv(f"{self.output_dir}{self.a['ExpName']}_depth.csv")

    def main(self):
        '''Entrypoint. Extract & merge probe lengths, reassign dupes if specified, then call anlysis & save'''
        end_sec_print("INFO: Analysis started.")
        probelengths = self.add_probelength()
        '''Depth calculation'''
        depth = self.add_depth(probelengths)
        '''Merge in sample info  (including total raw reads) and participant data if specified'''
        depth = self.add_read_d_and_clin(depth)
        if self.a["DebugMode"]:
            self.df["sampleid"] = self.df["sampleid"].astype(str)
            self.df = self.df.merge(
                depth, on=['sampleid', 'AGGREGATE'], how='left')
            self.df = self.df.rename(columns={"AGGREGATE": "probetype"})
            self.df.to_csv(
                f'{self.output_dir}/{self.a["ExpName"]}_fullself.df.csv.gz', index=False, compression='gzip')
        self.read_coverage_chart()
        self.read_dist_piechart()
        loginfo(
            f'Finished. Saved final data frame as {self.output_dir}/{self.a["ExpName"]}_fullself.df.csv.gz')
        end_sec_print("INFO: Analysis complete.")
