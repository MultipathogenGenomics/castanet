#!/bin/bash

################################################################################
# Download test set with ground truth labels for Castanet validation functions
# Labels/metadata file found here ("jiaa448_suppl_Supplementary_Data.xlsx"):
# https://doi.org/10.1093/infdis/jiaa448
################################################################################

DIR=data/rsv/
mkdir -p $DIR
STEM="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/"
TAIL=".fastq.gz"

declare -a DataIds=(
    "074/ERR10812874/ERR10812874"
    "075/ERR10812875/ERR10812875"
    "076/ERR10812876/ERR10812876"
    "077/ERR10812877/ERR10812877"
    "078/ERR10812878/ERR10812878"
    "079/ERR10812879/ERR10812879"
    "080/ERR10812880/ERR10812880"
    "081/ERR10812881/ERR10812881"
    "082/ERR10812882/ERR10812882"
    "083/ERR10812883/ERR10812883"
    "084/ERR10812884/ERR10812884"
    "085/ERR10812885/ERR10812885"
    "086/ERR10812886/ERR10812886"
    "087/ERR10812887/ERR10812887"
    "088/ERR10812888/ERR10812888"
    "089/ERR10812889/ERR10812889"
    "090/ERR10812890/ERR10812890"
)

# Pull paired reads
for val in ${DataIds[@]}; do
    read acc < <(echo $val | cut -d "/" -f 3)
    echo building set: $DIR$acc
    mkdir $DIR$acc
    wget -nc $STEM$val"_1"$TAIL -P $DIR$acc
    wget -nc $STEM$val"_2"$TAIL -P $DIR$acc
done
