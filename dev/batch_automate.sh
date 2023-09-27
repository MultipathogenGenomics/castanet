DIR=data/hiv_set/

for FILE in $DIR*; do
    read acc < <(echo $FILE | cut -d "/" -f 3)
    echo Hitting API for dataset: $acc
    curl -X 'POST' \
        'http://127.0.0.1:8001/end_to_end/' \
        -H 'accept: application/json' \
        -H 'Content-Type: application/json' \
        -d '{
        "ConsensusMinD": 10,
        "ConsensusCoverage": 30,
        "ConsensusMapQ": 1,
        "LineageFile": "data/ncbi_lineages_2023-06-15.csv.gz",
        "ExcludeIds": "9606",
        "RetainIds": "",
        "RetainNames": "",
        "ExcludeNames": "Homo,Alteromonas,Achromobacter",
        "KrakenDbDir": "kraken2_human_db/",
        "Probes": "data/2023_panel/all_seqs.csv",
        "KeepDups": true,
        "Clin": "",
        "DepthInf": "",
        "SamplesFile": "",
        "PostFilt": false,
        "RefStem": "data/2023_panel/all_seqs.fasta",
        "AdaptP": "data/all_adapters.fa",
        "ExpName": "test",
        "SeqName": "'$val'",
        "ExpDir": "./data/rsv_set/'$val'"
        }'
done
