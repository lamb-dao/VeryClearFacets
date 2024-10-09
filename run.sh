#!/usr/bin/env bash
#strict mode
set -euo pipefail
IFS=$'\n\t'

# ./run.sh NULL single input interactive mode
# ./run.sh full, single input noninteractive
# ./run.sh batch, multiple inputs in inputBatch dir, noninteractive

flag="${1:-default}"
if [ "${flag}" == "full" ]; then
    ./filter.sh batch 2>&1 | tee -a log.txt
elif [ "${flag}" == "batch" ]; then
    cd inputBatch
    for v in *.vcf; do
        cp -rav ${v} ../input
	cd ..
	rn=$(echo ${v} | sed s/.vcf// )
        ./filter.sh batch 2>&1 | tee -a log_${rn}.txt
	cd inputBatch
        rm ${v}
    done
else
    ./filter.sh 2>&1 | tee -a log.txt
fi
