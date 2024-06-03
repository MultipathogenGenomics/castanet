import argparse
import os
import sys
from app.utils.argparsemappings import map_args_filter_keep_reads, map_args_analysis


def parse_args_filter_keep_reads():
    '''Parse filter keep reads arguments, if script called from command line. Not used for API entry.'''
    parser = argparse.ArgumentParser(
        description=__doc__, epilog='Example: {progname} -i /path/to/my/fastq.gz -x 9606'.format(progname=os.path.basename(sys.argv[0])))
    parser.add_argument('-i', required=True, nargs='*',
                        help='Fastq[.gz] files to filter. Files must NOT contain blank lines. If gzipped, filenames must end in .gz. Each file will produce output with the suffix _filt.fastq in the current working directory.')
    parser.add_argument('-k', default=None, required=True,
                        help='Path to corresponding Kraken[.gz] output file, or other whitespace-separated text file where the second column is the read name and the third is the TaxID.')
    parser.add_argument(
        '-r', default='', help='NCBI TaxID(s) to retain. Use comma-separated list without spaces for multiple TaxIDs). Default is none.')
    parser.add_argument(
        '-x', default='', help='NCBI TaxID(s) to exclude. Use comma-separated list without spaces for multiple TaxIDs). Default is none. To remove all reads marked as human set to 9606.')
    parser.add_argument('--rT', default='', help='Text names of top-level taxa, to retain this and all taxa below it. Requires "--lineagefile". Eg. "--rT Bacteria,Viruses" to retain all bacterial and viral sequences.')
    parser.add_argument('--xT', default='', help='Text names of top-level taxa, to exclude this and all taxa below it. Requires "--lineagefile". Eg. "--xT \"Homo sapiens,Fungi\"" to exclude all human and fungal sequences.')
    parser.add_argument('--lineagefile', default='data/ncbi_lineages_2023-06-15.csv.gz',
                        help='Path to CSV file containing lineages of all NCBI taxa. Default is "lineages-2018-03-12.csv.gz".')
    return map_args_filter_keep_reads(parser.parse_args())


def parse_args_bam_parse():
    '''Parse args for filter and parse BAM functions. Programmatic only - no user interaction should be necessary.'''
    parser = argparse.ArgumentParser(
        description=__doc__, epilog='Do you want to parse or filter?')
    parser.add_argument('-ExpDir', required=False,
                        default="parse", help="Pass API arg in via shell.")
    parser.add_argument('-SingleEnded', required=True,
                        default="parse", help="Pass API arg in via shell.")
    parser.add_argument('-SeqName', required=True,
                        default="parse", help="Pass API arg in via shell.")
    parser.add_argument('-ExpName', required=False,
                        default="parse", help="Pass API arg in via shell.")
    parser.add_argument('-FilterFile', required=False,
                        default="parse", help="Pass API arg in via shell.")
    parser.add_argument('-Mode', required=True,
                        default="parse", help="parse or filter.")
    parser.add_argument('-ExpRoot', required=True,
                        default="./experiments", help="Pass API arg in via shell.")
    parser.add_argument('-MatchLength', required=True,
                        default="./experiments", help="Pass API arg in via shell.")
    return parser.parse_args()


def parse_args_analysis():
    '''Parse analysis command-line arguments. , if script called from command line. Not used for API entry.'''
    parser = argparse.ArgumentParser(
        description=__doc__, epilog='Example:-i /path/to/my/dataframe.csv.gz --samples /path/to/sampleinfo.csv')
    parser.add_argument(
        '-b', required=True, help='Batch name for these samples. Must be alphanumeric.')
    parser.add_argument(
        '-i', required=True, help='Data frame (csv[.gz]) file to process. If gzipped, filename must end in .gz.')
    parser.add_argument('-o', default=os.getcwd(),
                        help='Output directory. Default: current working directory.')
    parser.add_argument('-p', default='{}/probelengths_rmlst_virus_extra_ercc.csv'.format(os.path.dirname(sys.argv[0])),
                        help='Path to file containing probe lengths. Default: {}/probelengths_rmlst_virus_extra_ercc.csv'.format(os.path.dirname(sys.argv[0])))
    parser.add_argument('-d', '--keepdups', action='store_true', default=True,
                        help='If true, do not reassign duplicates to the sample with the majority in each duplicate cluster (Default: True).')
    parser.add_argument('--samples', required=True,
                        help='Path to CSV file containing information about raw reads (must have at least following fields: sampleid, pt, rawreadnum). Field "pt" must match clinical data.')
    parser.add_argument('--clin', required=False,
                        help='Path to CSV file containing clinical data (must have at least following fields: pt, clin_int; the field "sampleid" if present will be ignored). Other fields will be ignored.')
    parser.add_argument('--depth_inf', required=False,
                        help='(For regenerating full CSV with new clinical info): Path to previously generated CSV file of read depth per position for each probe, for all samples in this batch. Must contain the following fields: sampleid, target_id, depth_mean, depth_std, depth_25pc, depth_median, depth_75pc, prop_target_covered, prop_target_covered_mindepth2, prop_target_covered_mindepth5, prop_target_covered_mindepth10, udepth_mean, udepth_std, udepth_25pc, udepth_median, udepth_75pc, uprop_target_covered, uprop_target_covered_mindepth2, uprop_target_covered_mindepth5, uprop_target_covered_mindepth10')
    return map_args_analysis(parser.parse_args())
