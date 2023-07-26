from Bio import Entrez, SeqIO
import os


def DownloadGenBankFile(GenomeSeqFile, SeqIDLists, email):
    '''Hit GenBank to get data. Needs user to provide email authentication'''
    if not os.path.isfile(GenomeSeqFile):
        if not os.path.exists("/".join(GenomeSeqFile.split("/")[:-1])):
            os.makedirs("/".join(GenomeSeqFile.split("/")[:-1]))

        Entrez.email = email
        try:
            handle = Entrez.efetch(
                db="nucleotide", id=SeqIDLists, rettype="gb", retmode="text")
        except Exception as e:
            raise SystemExit(f"Failed to pull Genbank Data from NCBI Entrez with exception: {e}"
                             f"This is usually a temporary problem due to NCBI server down time, try again in a few minutes!")
        with open(GenomeSeqFile, "w") as GenomeSeqFile_handle:
            GenomeSeqFile_handle.write(handle.read())

        GenBankDict = SeqIO.index(GenomeSeqFile, "genbank")
        handle.close()
    else:
        GenBankDict = SeqIO.index(GenomeSeqFile, "genbank")
    return {k.split(".")[0]: v for k, v in GenBankDict.items()}
