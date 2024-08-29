import argparse


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
    parser.add_argument('-SaveDir', required=True,
                        default="./experiments", help="Pass API arg in via shell.")
    parser.add_argument('-MatchLength', required=True,
                        default="./experiments", help="Pass API arg in via shell.")
    return parser.parse_args()


def parse_arguments_lite():
    parser = argparse.ArgumentParser(
        description="Castanet Lite (Beta)"
    )
    '''N.b. Argparse is SO UNBELIEVABLY FUCKING SHIT that it can't evaluate booleans on optional arguments with a default, so we need to eval() strings passed to it later'''
    parser.add_argument('-Batch', required=False, default=False, type=str,
                        help="If True, conduct a batch run analysing multiple datasets; expects your ExpDir folder to contain sub-folders, each containing two (paired) read files.")
    parser.add_argument('-BAM', required=False, default=False, type=str,
                        help="If True, launch Castanet in BAM process mode, in which the software will skip initial processing and look in your input folder/s for BAM files rather than .fastq.gz")
    parser.add_argument('-ExpDir', required=True, type=str,
                        help="Folder containing your two paired read files, OR if batch = True, folder containing your experiment sub folders.")
    parser.add_argument('-ExpName', required=True, type=str,
                        help="Name for your experiment data.")
    parser.add_argument('-SaveDir', required=True, type=str,
                        help="Folder to save your experiment data.")
    parser.add_argument('-RefStem', required=True, type=str,
                        help="File location for mapping reference (fasta).")
    parser.add_argument('-DoKrakenPrefilter', required=False, default=True, type=str,
                        help="If True, do an initial Kraken pre-filter, using database specifid in -KrakenDbDir and list of NCBI TaxID(s) to exclude from -ExcludeIds.")
    parser.add_argument('-KrakenDbDir', required=False, default="kraken2_human_db", type=str,
                        help="If -DoKrakenPrefilter = True, path to Kraken2 database to do filtering.")
    parser.add_argument('-ExcludeIds', required=False, default="9606", type=str,
                        help="If -DoKrakenPrefilter = True, filter this/these NBCI TaxIDs (list of integers, separated by commas with no spaces)")
    parser.add_argument('-DoTrimming', required=False, default=True, type=str,
                        help="If True, use Trimmomatic to remove adapters and low quality sequences.")
    parser.add_argument('-DoConsensus', required=False, default=True, type=str,
                        help="If True, run the Castanet consensus sequence pipeline stage.")
    bool_fields = ["Batch", "BAM", "DoKrakenPrefilter",
                   "DoTrimming", "DoConsensus"]
    return parser, bool_fields


def parse_arguments_deptest():
    parser = argparse.ArgumentParser(
        description="Castanet Lite (Beta)"
    )
    parser.add_argument('-KrakenDbDir', required=False, default="kraken2_human_db/",
                        help="Directory path to your KrakenDbDir. If left blank, default will be Castanet's default install location for the human db.")
    return parser
