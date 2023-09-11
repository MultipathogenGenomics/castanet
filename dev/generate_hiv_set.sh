#!/bin/bash

################################################################################
# Download test set with ground truth labels for Castanet validation functions
# HUMAN IMMUNODEFICIENCY VIRUS
# Labels/metadata derived from "Supplementary Table S1_formatted":
# https://doi.org/10.1093/ve/vey007
################################################################################

DIR=data/hiv_set/
mkdir -p $DIR
STEM="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR732/"
TAIL=".fastq.gz"

declare -a DataIds=(
    "ERR732065/ERR732065"
    "ERR732066/ERR732066"
    "ERR732067/ERR732067"
    "ERR732068/ERR732068"
    "ERR732069/ERR732069"
    "ERR732070/ERR732070"
    "ERR732071/ERR732071"
    "ERR732072/ERR732072"
    "ERR732073/ERR732073"
    "ERR732074/ERR732074"
    "ERR732076/ERR732076"
    "ERR732077/ERR732077"
    "ERR732078/ERR732078"
    "ERR732079/ERR732079"
    "ERR732080/ERR732080"
    "ERR732081/ERR732081"
    "ERR732082/ERR732082"
    "ERR732083/ERR732083"
    "ERR732085/ERR732085"
    "ERR732086/ERR732086"
)

# Pull paired reads
for val in ${DataIds[@]}; do
    read acc < <(echo $val | cut -d "/" -f 2)
    echo acc is: $acc
    echo building set: $DIR$acc
    mkdir -p $DIR$acc
    wget -nc $STEM$val"_1"$TAIL -P $DIR$acc
    wget -nc $STEM$val"_2"$TAIL -P $DIR$acc
done
