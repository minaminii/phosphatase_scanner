#!/usr/bin/env bash

# Script: hmmscan_pipeline.sh
# Version: 1.0.0
#
# Author: Sara Wattabasombat
# Contact: sara_watt@cmu.ac.th
# Affiliation: Department of Microbiology, Faculty of Medicine, Chiang Mai University
#
# Description:
#   Run HMMER hmmscan against one or more HMM databases and summarize results.
#
# Last updated: 2026-01-21
#
# REQUIREMENTS:
#   - hmmscan (HMMER3)
#   - python3
#   - hmmscan_summary.py (located in the same directory as this script)
#
# USAGE:
#   hmmscan_pipeline.sh -i <protein_fasta> -p <prefix> -o <outdir> [options]
#

SCRIPT_NAME=$(basename "$0")
VERSION="1.0.0"
AUTHOR="W. Sara, Department of Microbiology, Faculty of Medicine, Chiang Mai University"
CONTACT="sara_watt@cmu.ac.th"

THREDS=4
prefix="untitled"

usage() {
    echo "${SCRIPT_NAME} v${VERSION}"
    echo "Author: ${AUTHOR}"
    echo "Contact: ${CONTACT}"
    echo
    echo "Usage:"
    echo "  ${SCRIPT_NAME} -i <protein_fasta> -p <prefix> -o <outdir> [options]"
    echo
    echo "Required arguments:"
    echo "  -i FILE   Input protein FASTA file"
    echo "  -p STR    Prefix for output files"
    echo "  -o DIR    Output directory"
    echo
    echo "Optional arguments:"
    echo "  -t INT    Number of CPU threads (default: ${THREDS})"
    echo "  -h        Show this help message and exit"
    echo
    echo "Example:"
    echo "  ${SCRIPT_NAME} \\"
    echo "      -i proteins.faa \\"
    echo "      -p sample1 \\"
    echo "      -o results/sample1 \\"
    echo "      -t 8"
    exit 1
}
OPTIONS=$(getopt -o i:p:o:t:ch \
  --long input:,prefix:,outdir:,threads:,config,help \
  -n "$0" -- "$@")
eval set -- "$OPTIONS"
while true; do
  case "$1" in
    -i|--input)
      input_protein_fasta="$2"
      shift 2
      ;;
    -p|--prefix)
      prefix="$2"
      shift 2
      ;;
    -o|--outdir)
      outdir="$2"
      shift 2
      ;;
    -t|--threads)
      THREADS="$2"
      shift 2
      ;;
    -c|--config)
      config=true
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    *)
      echo "ERROR: Invalid option $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "${input_protein_fasta:-}" || -z "${prefix:-}" || -z "${outdir:-}" ]]; then
    echo "ERROR: Missing required arguments" >&2
    usage
fi

if [[ ! -f "${input_protein_fasta}" ]]; then
    echo "ERROR: File not found: ${input_protein_fasta}" >&2
    exit 1
fi

BASE_DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
LOG_FILE="${outdir}/log"

if [[ ! -f $input_protein_fasta ]]; then
    echo "File not found: ${input_protein_fasta}"
    exit 1
fi

if [[ ! -d $outdir ]]; then
    mkdir -p $outdir
fi

for hmmdb in general PD; do
    hmmscan_out="${outdir}/${prefix}.${hmmdb}.hmmscan.out"
    hmmscan_tblout="${outdir}/${prefix}.${hmmdb}.hmmscan.tab"
    hmmscan_dom="${outdir}/${prefix}.${hmmdb}.hmmscan.dom.tab"
    echo "[RUNNING]: Running hmmscan with ${hmmdb}" | tee -a ${LOG_FILE}
    hmmscan \
        -o $hmmscan_out \
        --tblout $hmmscan_tblout \
        --domtblout $hmmscan_dom \
        --cpu $THREDS \
        "${BASE_DIR}/data/${hmmdb}.hmm3" \
        $input_protein_fasta | tee -a ${LOG_FILE}
    hmmscan_exit=$?
    if [[ hmmscan_exit -eq 0 ]]; then
        echo "[DONE]: hmmscan with ${hmmdb}.hmm profile" | tee -a ${LOG_FILE}
    else
        echo "[FAILED]: hmmscan with ${hmmdb}.hmm [exit: ${hmmscan_exit}]" | tee -a ${LOG_FILE}
    fi
done
echo "[RUNNING]: Generate HTML report" | tee -a ${LOG_FILE}
python3 $BASE_DIR/hmmscan_summary.py -t $prefix -tbl ""${outdir}/${prefix}.general.hmmscan.tab"" -fa $input_protein_fasta -o $outdir
