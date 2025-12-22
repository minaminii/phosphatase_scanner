#!/usr/bin/env bash

input_protein_fasta=$1
prefix=$2
outdir=$3

BASE_DIR=$(dirname "$0")
LOG_FILE="${outdir}/log"
THREDS=4

echo "The script you are running has:"
echo "basename: [$(basename "$0")]"
echo "dirname : [${BASE_DIR}]"
echo "pwd     : [$(pwd)]"

if [[ ! -f $input_protein_fasta ]]; then
    echo "File not found: ${input_protein_fasta}"
    exit 1
fi

if [[ ! -d $outdir ]]; then
    mkdir -p $outdir
fi

for hmmdb in general PD; do
    hmmscan \
        -o "${outdir}/${prefix}.${hmmdb}.hmmscan.out" \
        --tblout "${outdir}/${prefix}.${hmmdb}.hmmscan.tab" \
        --domtblout "${outdir}/${prefix}.${hmmdb}.hmmscan.dom.tab" \
        --cpu $THREDS \
        "${BASE_DIR}/data/${hmmdb}.hmm3" \
        $input_protein_fasta | tee -a ${LOG_FILE}
    hmmscan_exit=$?
    if [[ hmmscan_exit -eq 0 ]]; then
        echo "[CHECKPOINT]: hmmscan with ${hmmdb} [DONE]" | tee -a ${LOG_FILE}
    else
        echo "[CHECKPOINT]: hmmscan with ${hmmdb} [FAILED - exit: ${hmmscan_exit}]" | tee -a ${LOG_FILE}
    fi
done

if [[ $run_blastp ]]; then
    blastp
fi
