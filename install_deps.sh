# Install kraken2 via conda, for removing human reads
conda install -c bioconda kraken2

# Install Kallisto for quantifying target seqs
conda install -c bioconda kallisto

# Install samtools
# RM < TODO Install from binaries as conda install is rubbish
conda install -c bioconda samtools=1.9

# # Download pre-built kraken2 database with human genome only
mkdir kraken2_human_db
curl -L -o kraken2_human_db/kraken2_human_db.tar.gz https://ndownloader.figshare.com/files/23567780
tar -xzvf kraken2_human_db/kraken2_human_db.tar.gz

# Download trimmomatic and extract
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm Trimmomatic-0.39.zip
alias trimmomatic="java -jar Trimmomatic-0.39/trimmomatic-0.39.jar"

# Download BWA
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
alias bwa="bwa-mem2-2.2.1_x64-linux/bwa-mem2"


# Install viral consensus tool
apt install libbz2-dev libcurl4-openssl-dev liblzma-dev g++
git clone https://github.com/niemasd/ViralConsensus.git
cd ViralConsensus
make
sudo mv viral_consensus /usr/local/bin/
cd ..
