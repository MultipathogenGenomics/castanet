ADD=http://127.0.0.1:8001/check_dependencies/

resp=$(curl -X 'POST' \
    ''$ADD'' \
    -H 'accept: application/json' \
    -H 'Content-Type: application/json' \
    -d '{
        "AdaptP": "data/all_adapters.fa",
        "KrakenDbDir": "kraken2_human_db/"
    }')
echo Completed processing: $resp
