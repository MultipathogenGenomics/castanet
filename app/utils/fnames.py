'''File names for persistent file storage'''


def get_consensus_fnames(args):
    return {
        "master_bam": f"{args['folder_stem']}{args['SeqName']}.bam",
        "flat_cons_refs": f"{args['folder_stem']}consensus_data/temp_refs.fasta",
        "flat_cons_seqs": f"{args['folder_stem']}consensus_data/temp_seqs.fasta",
        "grouped_reads_dir": f"{args['folder_stem']}grouped_reads/",
        "collated_reads_bam": f"{args['folder_stem']}consensus_data/collated_reads.bam",
        "temp_folder": f"{args['folder_stem']}/tempfolder/",
        "temp_ref_seq": f"{args['folder_stem']}/tempfolder/refseq.fasta",
        "collated_reads_fastq": f"{args['folder_stem']}consensus_data/collated_reads.fastq"
    }
