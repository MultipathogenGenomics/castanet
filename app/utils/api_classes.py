from fastapi import Query
from typing import Union, Literal, Optional
from pydantic import BaseModel, Field, FilePath, DirectoryPath, validator

'''Primitives'''


class Data_BatchName(BaseModel):
    BatchName: DirectoryPath = Query('./data/my_experiments/',
                                     description="Path to recursively read for individual datasets.")


class Data_ExpDir(BaseModel):
    ExpDir: DirectoryPath = Query('./data/',
                                  description="Path to retrieve and store data.")


class Data_AdaptP(BaseModel):
    AdaptP: FilePath = Query('data/all_adapters.fa',
                             description='Location of your Trimmomatic adapter sequences - may be in your Trimmomatic path, but a backup is in the data dir.')


class Data_RefStem(BaseModel):
    RefStem: FilePath = Query("data/my_probes.fasta",
                              description="Path to mapping file, in fasta format.")


class Data_PostFilt(BaseModel):
    PostFilt: bool = Query(False,
                           description="Post hoc filter BAM file to remove reads marked as contaminations")


class Data_GenerateCounts(BaseModel):
    SingleEndedReads: bool = Query(False,
                                   description="Set to true if using single-ended reads, e.g. if sequencing run ended half-way through.")


class Data_NThreads(BaseModel):
    NThreads: Union[int, str] = Query('auto',
                                      description="Specify the number of threads for multi-core processing. Options: integer == this many threads; 'auto' == let Castanet choose number of threads; 'hpc' == select when running on a compute cluster (hard codes to 1).")


class Data_ExpName(BaseModel):
    ExpName: str = Query('myexperiment',
                         description="Name your experiment/batch")


class Data_KrakenDir(BaseModel):
    KrakenDbDir: DirectoryPath = Query('kraken2_human_db/',
                                       description="Path to Kraken2 database for filtering human/other unwanted species reads. Only used if DoKrakenPrefilter = true")


class Data_FilterFilters(BaseModel):
    DoKrakenPrefilter: bool = Query(True,
                                    description="If true, run an initial pre-filtering step to label reads with Kraken2 and exclude those belonging to taxonomies indicated in Exclude/Retain IDs/Names arguments. If false, these other fields are ignored.")

    LineageFile: Union[None, str] = Query('data/ncbi_lineages_2023-06-15.csv.gz',
                                          description="(OPTIONAL) Path to CSV file containing lineages of all NCBI taxa. Only used if DoKrakenPrefilter = true.")

    ExcludeIds: Union[None, str] = Query("9606",
                                         description="(OPTIONAL) Exclude these NCBI TaxID/s from filter keep reads step. Comma separate without spaces. Set to 9606 to exclude Human. Only used if DoKrakenPrefilter = true.")

    RetainIds: Union[None, str] = Query("",
                                        description="(OPTIONAL) Exclude these NCBI TaxID/s from filter keep reads step. Comma separate without spaces. Only used if DoKrakenPrefilter = true.")

    RetainNames: Union[None, str] = Query("",
                                          description="(OPTIONAL) Retain these species names from filter keep reads step. Comma separate without spaces. Will be ignored if no Linneage file specified. Only used if DoKrakenPrefilter = true.")

    ExcludeNames: Union[None, str] = Query("Homo",
                                           description="(OPTIONAL) Exclude these species names from filter keep reads step. Comma separate without spaces. Will be ignored if no Linneage file specified. Only used if DoKrakenPrefilter = true.")


class Data_TrimmomaticParams(BaseModel):
    DoTrimming: bool = Query(True,
                             description="If true, use Trimmomatic to remove low quality reads and adapters. Minimum trim length can be set with TrimMinLen argument; if false, this extra argument is ignored.")
    TrimMinLen: int = Query(36,
                            description="Values < than min trim length will be removed by Trimmomatic tool. Only used if DoTrimming = true")


class Data_AnalysisExtras(BaseModel):
    KeepDups: bool = Query(True,
                           description='(OPTIONAL) If true, do not reassign duplicates to the sample with the majority in each duplicate cluster (Default: True).')
    Clin: Optional[str] = Query("",
                                description='(OPTIONAL) Path to CSV file containing clinical data (must have at least following fields: pt, clin_int; the field "sampleid" if present will be ignored). Other fields will be ignored.')
    DepthInf: str = Query("",
                          description='(OPTIONAL, For regenerating full CSV with new clinical info): Path to previously generated CSV file of read depth per position for each probe, for all samples in this batch. Must contain the following fields: sampleid, target_id, depth_mean, depth_std, depth_25pc, depth_median, depth_75pc, prop_target_covered, prop_target_covered_mindepth2, prop_target_covered_mindepth5, prop_target_covered_mindepth10, udepth_mean, udepth_std, udepth_25pc, udepth_median, udepth_75pc, uprop_target_covered, uprop_target_covered_mindepth2, uprop_target_covered_mindepth5, uprop_target_covered_mindepth10')
    SamplesFile: Optional[str] = Query("",
                                       description="(OPTIONAL) If specified, read raw read numbers from this CSV (needs cols 'sampleid', 'pt', 'rawreadnum'). If not specified, CASTANET will read the raw read numbers from the input bam file, i.e. it will assume you haven't pre-filtered the file.")


class Data_ConsensusParameters(BaseModel):
    ConsensusMinD: int = Query(10,
                               description="Minimum base depth required to make a call, for consensus calling functions")
    ConsensusCoverage: float = Query(30.0,
                                     description="Do not generate consensus if coverage < n. Applies to both target consensuses and final, remapped consensus.")
    ConsensusMapQ: float = Query(1.0,
                                 description="Minimum quality value for a target consensus to be included in the remapped consensus.")
    ConsensusCleanFiles: bool = Query(True,
                                      description="If True, consensus generator will delete BAM files for reads aggregated to each target organism. Disable to retain files for use in downstream analysis.")
    GtFile: Optional[str] = Query('',
                                  description="(OPTIONAL - EVAL MODE ONLY) CSV file containing at least columns: `Primary_accession` and `GenBank_accession` for evaluating consensus seqs vs ground truth")
    GtOrg: Optional[str] = Query('',
                                 description="(OPTIONAL - EVAL MODE ONLY) Name of target organism to measure ground truth sequence against.")


'''Endpoint objects'''


class E2e_data(Data_NThreads, Data_AdaptP,
               Data_PostFilt, Data_AnalysisExtras, Data_KrakenDir, Data_FilterFilters,
               Data_ConsensusParameters, Data_TrimmomaticParams, Data_GenerateCounts,
               Data_RefStem, Data_ExpName, Data_ExpDir):
    pass


class E2e_eval_data(Data_ExpDir, Data_ExpName, Data_NThreads, Data_AdaptP, Data_RefStem,
                    Data_PostFilt, Data_AnalysisExtras, Data_KrakenDir, Data_FilterFilters,
                    Data_ConsensusParameters, Data_TrimmomaticParams, Data_GenerateCounts):
    pass


class Bam_workflow_data(Data_NThreads, Data_PostFilt, Data_AnalysisExtras, Data_ConsensusParameters,
                        Data_GenerateCounts, Data_RefStem, Data_ExpName, Data_ExpDir):
    pass


class Preprocess_data(Data_ExpDir, Data_ExpName, Data_NThreads, Data_KrakenDir, Data_TrimmomaticParams):
    pass


class Filter_keep_reads_data(Data_ExpDir, Data_ExpName, Data_NThreads, Data_FilterFilters, Data_TrimmomaticParams):
    pass


class Trim_data(Data_ExpDir, Data_ExpName, Data_NThreads, Data_AdaptP, Data_TrimmomaticParams):
    pass


class Mapping_data(Data_ExpDir, Data_ExpName, Data_NThreads, Data_RefStem):
    pass


class Count_map_data(Data_ExpDir, Data_ExpName, Data_NThreads, Data_GenerateCounts):
    pass


class Analysis_data(Data_ExpDir, Data_ExpName, Data_NThreads, Data_AnalysisExtras, Data_RefStem):
    pass


class Post_filter_data(Data_ExpDir, Data_NThreads, Data_ExpName):
    pass


class Batch_eval_data(Data_BatchName, Data_ExpName, Data_NThreads, Data_AdaptP, Data_RefStem,
                      Data_PostFilt, Data_AnalysisExtras, Data_KrakenDir,
                      Data_FilterFilters, Data_ConsensusParameters, Data_TrimmomaticParams, Data_GenerateCounts):
    pass


class Consensus_data(Data_ExpName, Data_NThreads, Data_RefStem, Data_ConsensusParameters):
    pass


class Eval_data(Data_ExpName, Data_NThreads, Data_RefStem,
                Data_ConsensusParameters):
    pass


class Dep_check_data(Data_KrakenDir, Data_AdaptP):
    pass


class Convert_probe_data(BaseModel):
    InputFolder: str = Query("",
                             description="Folder path, wherein all .fasta files will be imported, collated and converted to a Castanet-compatible probe set.")
    OutFolder: str = Query("",
                           description="Output path for collated probes fasta and probe length csv")
    OutFileName: str = Query("",
                             description="Output file name for collated probes fasta and probe length csv")
