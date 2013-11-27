#!/bin/bash

VARIANTDB=/data/apagan/variantdb
REDISDIR=/home/apagan/code/redis-2.6.16/src

for file in "$VARIANTDB"/*.vd; do
    echo "Loading $file into redis"
    python redis_inserts_for_vardb.py $file | "$REDISDIR/"redis-cli --pipe
done
