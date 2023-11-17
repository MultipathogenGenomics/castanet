# Adjust channel priority on conda env
conda config --add channels bioconda
conda config --add channels conda-forge

# Install kraken2 via conda, for removing human reads
conda install -y -c bioconda kraken2

# Install samtools
conda install -y "samtools>=1.10"

# Install mafft
conda install -y -c bioconda mafft

# Install MASH
conda install -y -c bioconda mash

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
git clone --recurse-submodules https://github.com/niemasd/ViralConsensus.git
cd ViralConsensus
make
sudo mv viral_consensus /usr/local/bin/
cd ..
alias viral_consensus="./ViralConsensus/viral_consensus"
