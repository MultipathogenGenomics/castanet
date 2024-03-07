#!/bin/bash
curl -X 'POST' \
'http://127.0.0.1:8001/end_to_end/' \
-H 'accept: application/json' \
-H 'Content-Type: application/json' \
-d '{
"ExpDir": "./data/eval/",
"ExpName": "CastanetTest",
"RefStem": "data/eval/ref.fa",
"SingleEndedReads": false,
"DoTrimming": true,
"TrimMinLen": 36,
"ConsensusMinD": 10,
"ConsensusCoverage": 30,
"ConsensusMapQ": 1,
"ConsensusCleanFiles": true,
"GtFile": "",
"GtOrg": "",
"DoKrakenPrefilter": true,
"LineageFile": "data/ncbi_lineages_2023-06-15.csv.gz",
"ExcludeIds": "9606",
"RetainIds": "",
"RetainNames": "",
"ExcludeNames": "Homo",
"KrakenDbDir": "kraken2_human_db/",
"KeepDups": true,
"Clin": "",
"DepthInf": "",
"SamplesFile": "",
"PostFilt": false,
"AdaptP": "data/all_adapters.fa",
"NThreads": "auto"
}'
