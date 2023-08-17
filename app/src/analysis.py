from __future__ import division
import os
import re
import numpy as np
import pandas as pd
import plotly.express as px

from app.utils.system_messages import end_sec_print
from app.utils.argparsers import parse_args_analysis
from app.utils.shell_cmds import loginfo, stoperr, logerr
from app.utils.error_handlers import error_handler_analysis
from app.utils.basic_cli_calls import samtools_read_num


class Analysis:
    def __init__(self, argies, api_entry=True) -> None:
        self.a = argies
        self.output_dir = f"experiments/{self.a['ExpName']}/"
        if api_entry:
            self.a["input_file"] = f"{self.output_dir}/{self.a['SeqName']}_PosCounts.csv"
        self.df = error_handler_analysis(self.a)

    def add_probelength(self):
        '''Add length of target_id to each row of master df after splitting probelength data.'''
        loginfo(f'Adding probe length information from {self.a["Probes"]}')
        try:
            probelengths = pd.read_csv(
                self.a["Probes"], dtype={'target_len': int})
        except:
            stoperr(
                f'Failed to read probe information. Is {self.a["Probes"]} a valid CSV file?')
        probelengths_mod = self.add_probetype(probelengths)
        self.df = self.df.merge(
            probelengths_mod, left_on='target_id', right_on='target_id', how='left')
        return probelengths_mod

    def add_probetype(self, pdf):
        ''' For probe aggregation, determine organism and gene ID.
        Splits genename and probetype into separate cols, then does manual adjustments
        '''
        loginfo('Aggregating by organism and gene name.')
        pdf['genename'] = pdf.target_id.apply(
            lambda x: x.replace('_', '-').split('-')[0])
        pdf.loc[pdf.target_id == 'roseolovirus_allrecords_cluster_1',
                'genename'] = 'HHV7_roseolovirus_allrecords_cluster_1'
        pdf.loc[pdf.target_id == 'roseolovirus_allrecords_cluster_2',
                'genename'] = 'HHV6_roseolovirus_allrecords_cluster_2'

        '''More precise definition for the different virus types'''
        pdf.loc[pdf.genename == 'enterovirus', 'genename'] = pdf.loc[pdf.genename ==
                                                                     'enterovirus'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))
        pdf.loc[pdf.genename == 'Coronaviridae', 'genename'] = pdf.loc[pdf.genename ==
                                                                       'Coronaviridae'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:3]))
        pdf.loc[pdf.genename == 'Adenoviridae', 'genename'] = pdf.loc[pdf.genename ==
                                                                      'Adenoviridae'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))
        pdf.loc[pdf.genename == 'Flaviviridae', 'genename'] = pdf.loc[pdf.genename ==
                                                                      'Flaviviridae'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))
        pdf.loc[pdf.genename == 'Influenza', 'genename'] = pdf.loc[pdf.genename ==
                                                                   'Influenza'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))
        pdf.loc[pdf.genename == 'Paramyxoviridae', 'genename'] = pdf.loc[pdf.genename ==
                                                                         'Paramyxoviridae'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))
        pdf.loc[pdf.genename == 'Parvoviridae', 'genename'] = pdf.loc[pdf.genename ==
                                                                      'Parvoviridae'].target_id.apply(lambda x: '_'.join(x.replace('_', '-').split('-')[:2]))
        pdf['probetype'] = pdf.genename
        pat = re.compile('BACT[0-9]+_([A-Za-z]+)-[0-9]+[|_]([A-Za-z]+)')
        backup_pat = re.compile('BACT[0-9]+_[0-9]+_[A-Z]+_([A-Za-z_]+)')

        def _pat_search(s, pat=pat, backup_pat=backup_pat):
            '''Private function to return empty string instead of error when pattern is not matched.'''
            res = pat.findall(s)
            if not res:
                res = (backup_pat.findall(s),)
            if not res:
                return ''
            return '_'.join(res[0])

        pdf.loc[pdf.genename.str.startswith('BACT'), 'probetype'] = pdf.loc[pdf.genename.str.startswith(
            'BACT')].target_id.apply(_pat_search)
        '''Streptococcus mitis group (pneumo, mitis, oralis) are cross-mapping, so classify as S. pneumoniae all targets that are found in S.pneumo at least once in the database'''
        pdf.loc[(pdf.target_id.apply(lambda x: 'pneumoniae' in x)) & (pdf.probetype.isin(('Streptococcus_pneumoniae',
                                                                                          'Streptococcus_pseudopneumoniae', 'Streptococcus_mitis', 'Streptococcus_oralis'))), 'probetype'] = 'Streptococcus_pneumoniae'
        '''Streptococcus_agalactiae and Streptococcus_pyogenes cross-map as well, aggregate them'''
        pdf.loc[(pdf.target_id.apply(lambda x: 'pyogenes' in x or 'agalactiae' in x)) & (pdf.probetype.isin(
            ('Streptococcus_pyogenes', 'Streptococcus_agalactiae'))), 'probetype'] = 'Streptococcus_agalactiae_pyogenes'
        '''Enterobacteriacae are not distinguishable at this level so group them all'''
        pdf.loc[pdf.probetype.apply(lambda x: x.startswith('Escherichia') or x.startswith(
            'Klebsiella') or x.startswith('Enterobacter')), 'probetype'] = 'Enterobacteriaceae'
        loginfo(
            f'Organism and gene summary: {pdf.probetype.nunique()} organisms, up to {pdf.groupby("probetype").genename.nunique().max()} genes each.')
        pdf.to_csv(f"{self.output_dir}/probe_aggregation.csv")
        return pdf

    def reassign_dups(self):
        ''' Reassign positional duplicates (reads with same start/end) to the sample with the highest count.
        This should NOT be used with unique dual indexed samples, as relies on heavily oversequenced data.
        Regardless may be useful to clean contaminated data.'''
        loginfo('Reassaigning duplicates.')
        ttl_prededup_mapped_reads = self.df.groupby('sampleid').n.sum()
        grouped_prededup_mapped_reads = self.df.groupby(
            ['sampleid', 'probetype']).n.sum().unstack().fillna(0)
        self.df.sort_values(
            ['target_id', 'startpos', 'maplen', 'n'], ascending=True, inplace=True)
        self.df['n'] = self.df.groupby(
            ['target_id', 'startpos', 'maplen']).n.cumsum()
        nrows0 = len(self.df)
        duprows = self.df.duplicated(
            ['target_id', 'startpos', 'maplen'], keep='last')
        loginfo(f'Saving {self.a["ExpName"]}_reads_to_drop.csv')
        self.df[duprows][['sampleid', 'target_id', 'startpos', 'maplen']].to_csv(
            f'{self.a["ExpDir"]}/{self.a["ExpName"]}_reads_to_drop.csv', index=False)
        self.df = self.df[~duprows]
        loginfo(
            f'Kept {len(self.df)} of {nrows0} rows ({float(len(self.df))/nrows0}).')
        ttl_dup_mapped_reads = self.df.groupby('sampleid').n.sum()
        grouped_dup_mapped_reads = self.df.groupby(
            ['sampleid', 'probetype']).n.sum().unstack().fillna(0)
        try:
            assert ttl_dup_mapped_reads.sum() == ttl_prededup_mapped_reads.sum()
        except AssertionError:
            stoperr(
                'Total reads after duplicate reassignment does not match starting total.')
        duprate_ttl = ttl_dup_mapped_reads/ttl_prededup_mapped_reads
        duprate_by_probetype = grouped_dup_mapped_reads/grouped_prededup_mapped_reads
        loginfo(
            f'Saving {self.a["ExpName"]}_duprate.csv and {self.a["ExpName"]}_duprate_by_probetype.csv.')
        duprate_ttl.to_csv(
            f'{self.a["ExpDir"]}/{self.a["ExpName"]}_duprate.csv', header=True)
        duprate_by_probetype.to_csv(
            f'{self.a["ExpDir"]}/{self.a["ExpName"]}_duprate_by_probetype.csv')
        loginfo(f'Duplication rate: \n{duprate_ttl.describe().to_string()}')
        return

    def add_depth(self, probelengths):
        ''' Calculate read depth per position. '''
        if self.a["DepthInf"]:
            loginfo(
                f'Reading read depth information from {self.a["DepthInf"]}.')
            try:
                depth = pd.read_csv(self.a["DepthInf"])
                assert len(depth) > 0
            except (IOError, AssertionError, TypeError):
                logerr(
                    'Failed to read read depth information from file. Will calculate from data.')
                depth = None
        else:
            loginfo('Calculating read depth information.')
            depth = None

        if not depth:
            metrics = {}
            odir = f'{self.output_dir}/Depth_output'
            if not os.path.isdir(odir):
                try:
                    os.mkdir(odir)
                except OSError:
                    odir = os.getcwd()
                    logerr(
                        f'Cannot create output directory {odir} for saving depth information. Proceeding with current working directory.')

            loginfo(
                'INFO: Calculating read depth statistics for all probes, for all samples. This is a slow step.')
            for (sampleid, probetype), g in self.df.groupby(['sampleid', 'probetype']):
                gene_list = g.genename.unique()
                n_genes = len(gene_list)
                target_list = g.target_id.unique()
                n_targets = len(target_list)
                Dd = {}
                D1d = {}
                for genename, gg in g.groupby('genename'):
                    loginfo(f'Processing {sampleid} - {genename}')
                    '''Generate two arrays for each probetype-genename group (D = number of occurrences, D1 = unique start/end positions)'''
                    D = np.zeros(int(gg.target_len.max()), dtype=np.uint32)
                    D1 = np.zeros(D.shape, dtype=np.uint32)
                    D.fill(0)
                    D1.fill(0)  # Clear target arrays
                    for target_id, gt in gg.groupby('target_id'):
                        loginfo(f'..... target: {target_id[:100]}')
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
                amprate = (D/D1)[~np.isnan(D/D1)]
                '''Max possible positions for the genes that were actually in this BAM (accounts for some genes not being captured)'''
                npos = len(D)
                '''Now pad out with zeros to the total number of mappable positions for this probetype (nmax_probetype above)'''
                D = np.pad(D, (0, nmax_probetype - npos),
                           'constant', constant_values=0)
                D1 = np.pad(D1, (0, nmax_probetype - npos),
                            'constant', constant_values=0)
                loginfo(f'Mean depth (all reads) for {probetype}: {D.mean()}')
                loginfo(
                    f'Mean depth (deduplicated) for {probetype}: {D1.mean()}')
                '''Amplification rate calculations use the unpadded (mapped) number of sites as the denominator'''
                loginfo(
                    f'Mean amplification ratio for {probetype}: {amprate.mean()}')

                '''Save arrays as CSV'''
                with open(f'{odir}/{probetype}-{sampleid}_depth_by_pos.csv', 'a') as o:
                    np.savetxt(o, D, fmt='%d', newline=',')
                    o.write('\n')
                    np.savetxt(o, D1, fmt='%d', newline=',')
                    o.write('\n')

                '''Save array plots as pdf if significant'''
                if D1.mean() >= 0.1:
                    plot_df = pd.DataFrame()
                    plot_df["position"], plot_df["All Reads"], plot_df["Deduplicated Reads"] = np.arange(
                        0, D.shape[0]), D, D1
                    fig = px.line(plot_df, x="position", y=[
                                  "All Reads", "Deduplicated Reads"], title=f'{sampleid}\n{probetype} ({n_targets}/{nmax_targets} targets in {n_genes}/{nmax_genes} genes)',
                                  labels={"position": "Position", "value": "Num Reads"})
                    fig.update_layout(legend={"title_text": "", "orientation": "h", "entrywidth": 100,
                                      "yanchor": "bottom", "y": 1.02, "xanchor": "right", "x": 1})
                    fig.write_image(f'{odir}/{probetype}-{sampleid}.png')

                '''Build up dictionary of depth metrics for this sample and probetype'''
                metrics[sampleid, probetype] = (g.n.sum(), g.n.count(), n_targets, n_genes, nmax_targets, nmax_genes, nmax_probetype, npos,
                                                amprate.mean(), amprate.std(), np.median(amprate),
                                                D.mean(), D.std(), np.percentile(D, 25), np.median(D), np.percentile(D, 75),
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
            depth = pd.DataFrame(metrics, index=['n_reads_all', 'n_reads_dedup', 'n_targets', 'n_genes', 'nmax_targets', 'nmax_genes', 'npos_max_probetype', 'npos_cov_probetype',
                                                 'amprate_mean', 'amprate_std', 'amprate_median',
                                                 'depth_mean', 'depth_std', 'depth_25pc', 'depth_median', 'depth_75pc',
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
                         'level_1': 'probetype'}, inplace=True)

            '''Add reads on target (rot)'''
            try:
                rot = depth.groupby('sampleid').n_reads_all.sum().reset_index()
            except:
                stoperr(
                    f'ERROR: Depth dataframe seems to be empty. Please check trimming and mapping steps completed successfully.')

            rot.rename(
                columns={'n_reads_all': 'reads_on_target'}, inplace=True)
            depth = depth.merge(rot, on='sampleid', how='left')
            rot_dedup = depth.groupby(
                'sampleid').n_reads_dedup.sum().reset_index()
            rot_dedup.rename(
                columns={'n_reads_dedup': 'reads_on_target_dedup'}, inplace=True)
            depth = depth.merge(rot_dedup, on='sampleid', how='left')
            depth['prop_of_reads_on_target'] = depth.n_reads_all / \
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
            depth['clean_n_reads_all'] = depth.n_reads_all-depth.n_reads_dedup
            depth['clean_prop_of_reads_on_target'] = (
                depth.n_reads_all-depth.n_reads_dedup)/depth.reads_on_target
            '''Save and log'''
            loginfo(f'Saving {self.output_dir}/{self.a["ExpName"]}_depth.csv.')
            depth.to_csv(
                f'{self.output_dir}/{self.a["ExpName"]}_depth.csv', index=False)
            loginfo(
                f'Mean read depth per sample: \n{depth.groupby("sampleid").depth_mean.mean().to_string()}')
        return depth

    def add_clin(self, req_cols_clin=['pt', 'clin_int']):
        '''Merge to create simpler metadata for each sample, including patient ID and clinical category.
        Clin file may additionally supply clinical/demographic data,
        with at least the following columns: {}'''.format(req_cols_clin)
        clin = pd.read_csv(self.a["Clin"], dtype={'pt': str})
        if 'pt' not in clin.columns:
            stoperr(
                f'Clin file {self.a["Clin"]} must contain at least the following columns: {req_cols_clin}')
        if 'sampleid' in clin.columns:
            logerr(
                f'Ignoring "sampleid" column from clinical info file {self.a["Clin"]}.')
            clin.drop('sampleid', axis=1, inplace=True)
            '''Merge sample information with participant data'''
            samples = samples.merge(clin, on='pt', how='left')
        return samples

    def add_read_d_and_clin(self, depth, req_cols_samples=['sampleid', 'pt', 'rawreadnum']):
        ''' Add raw read numbers and any external categorical/clinical data.
        Samples file must supply at least the following columns: {}.'''.format(req_cols_samples)
        loginfo('Adding sample information and clinical data.')
        read_num = samtools_read_num(self.output_dir, self.a["SeqName"])
        samples = pd.DataFrame(
            [{"sampleid": self.a["SeqName"], "pt": "", "rawreadnum": read_num}])

        if self.a["Clin"] != "":
            '''If supplied, merge clinical data'''
            samples = self.add_clin(samples)

        '''Merge read n (and clin data if supplied) to depth counts, return'''
        cdf = depth.merge(samples, on='sampleid', how='left')
        cdf['readprop'] = cdf.n_reads_all/cdf.rawreadnum
        loginfo(f'Added the following columns: {list(samples.columns)}')
        loginfo(
            f'Saving {self.output_dir}{self.a["SeqName"]}_depth_with_clin.csv.')
        cdf.to_csv(
            f'{self.output_dir}{self.a["SeqName"]}_depth_with_clin.csv', index=False)
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

    def main(self):
        '''Entrypoint. Extract & merge probe lengths, reassign dupes if specified, then call anlysis & save'''
        end_sec_print("INFO: Analysis started.")
        probelengths = self.add_probelength()
        if not self.a["KeepDups"]:
            self.reassign_dups()
        '''Depth calculation'''
        depth = self.add_depth(probelengths)
        '''Merge in sample info  (including total raw reads) and participant data if specified'''
        depth = self.add_read_d_and_clin(depth)
        self.save_tophits(depth)
        self.df = self.df.merge(
            depth, on=['sampleid', 'probetype'], how='left')
        self.df.to_csv(
            f'{self.output_dir}/{self.a["ExpName"]}_fullself.df.csv.gz', index=False, compression='gzip')
        loginfo(
            f'Finished. Saved final data frame as {self.output_dir}/{self.a["ExpName"]}_fullself.df.csv.gz')
        end_sec_print("INFO: Analysis complete.")


if __name__ == '__main__':
    '''CLI entry'''
    cls = Analysis(parse_args_analysis(), api_entry=False)
    cls.main()
    end_sec_print("INFO: Analysis complete.")
