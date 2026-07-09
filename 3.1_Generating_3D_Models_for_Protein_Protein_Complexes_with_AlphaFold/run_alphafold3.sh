#!/bin/bash
set -euo pipefail

# CASE INPUTS
# Case ID
export CASE="${CASE:-4POU}"
#export CASE="${CASE:-2FJG}"

# Protein sequences separated by ":".
# Each ":"-separated block will be treated as a different chain:
#   first sequence  -> chain A
#   second sequence -> chain B
#   third sequence  -> chain C
#   etc.
export SEQ="${SEQ:-QVQLVESGGGLVQAGGSLRLSCAASGYPHPYLHMGWFRQAPGKEREGVAAMDSGGGGTLYADSVKGRFTISRDKGKNTVYLQMDSLKPEDTATYYCAAGGYQLRDRTYGHWGQGTQVTVSS:KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQSYSTMSITDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV}"
#export SEQ="${SEQ:-EVQLVESGGGLVQPGGSLRLSCAASGFTISDYWIHWVRQAPGKGLEWVAGITPAGGYTYYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARFVFFLPYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSC:DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNR:EVVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIGEMSFLQHNKCECRPK:EVVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIGEMSFLQHNKCECRPK}"

# AlphaFold 3 model seeds.
# Use one seed or several comma-separated seeds.
export MODEL_SEEDS="${MODEL_SEEDS:-13441}"

# SOFTWARE AND PATHS
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APPDIR="${APPDIR:-/home/user/Programs}"
ALPHAFOLD3DIR="${ALPHAFOLD3DIR:-${APPDIR}/alphafold3}"
CIF_TO_PDB="${CIF_TO_PDB:-${SCRIPT_DIR}/cif_to_pdb.py}"
# HMMER binaries.
# Default assumes HMMER was installed in the active conda environment.
# To use system binaries, run for example:
#   HMMER3_BINDIR=/usr/bin ./run_af3_case.sh
HMMER3_BINDIR="${HMMER3_BINDIR:-${CONDA_PREFIX:-}/bin}"

DB_DIR="${DB_DIR:-${ALPHAFOLD3DIR}/public_databases}"
MODEL_DIR="${MODEL_DIR:-${ALPHAFOLD3DIR}/models}"

# AlphaFold 3 runner
RUN_AF3="${RUN_AF3:-run_alphafold.py}"

# AlphaFold 3 buckets
BUCKETS="${BUCKETS:-256,512,768,1024,1280,1536,2048,2560,3072,3584,4096,4608,5120}"

# WORKING DIRECTORIES
WORK_DIR="$(pwd)"
CASE_DIR="${WORK_DIR}/${CASE}/AF3"
JSON_FILE="${CASE_DIR}/${CASE}.json"
OUTPUT_DIR="${CASE_DIR}/output"
LOG_FILE="${OUTPUT_DIR}/af3_run.log"

mkdir -p "${CASE_DIR}"
mkdir -p "${OUTPUT_DIR}"

# BASIC CHECKS
if [[ -z "${SEQ}" ]]; then
    echo "ERROR: SEQ is empty." >&2
    exit 1
fi

if [[ -z "${HMMER3_BINDIR}" || ! -d "${HMMER3_BINDIR}" ]]; then
    echo "ERROR: HMMER3_BINDIR does not exist or is empty: ${HMMER3_BINDIR}" >&2
    echo "Activate the conda environment or define HMMER3_BINDIR manually." >&2
    exit 1
fi

for bin in jackhmmer nhmmer hmmalign hmmsearch hmmbuild; do
    if [[ ! -x "${HMMER3_BINDIR}/${bin}" ]]; then
        echo "ERROR: could not find executable ${bin} in ${HMMER3_BINDIR}" >&2
        exit 1
    fi
done

if [[ ! -d "${DB_DIR}" ]]; then
    echo "ERROR: DB_DIR does not exist: ${DB_DIR}" >&2
    exit 1
fi

if [[ ! -d "${MODEL_DIR}" ]]; then
    echo "ERROR: MODEL_DIR does not exist: ${MODEL_DIR}" >&2
    exit 1
fi

# GENERATE ALPHAFOLD 3 JSON INPUT
python - <<PY
import json
import re
from pathlib import Path

case = "${CASE}"
seq_raw = """${SEQ}"""
seeds_raw = """${MODEL_SEEDS}"""
json_file = Path("${JSON_FILE}")

# Split input sequences by ":"
chains = [
    seq.strip().replace(" ", "").replace("\\n", "")
    for seq in seq_raw.split(":")
    if seq.strip()
]

if not chains:
    raise SystemExit("ERROR: no sequences were found in SEQ.")

if len(chains) > 26:
    raise SystemExit("ERROR: this automatic generator only supports up to 26 chains: A-Z.")

# Basic protein sequence validation.
# Extend this block if you later need ligands, DNA, RNA, glycans, or custom modifications.
aa_pattern = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXBZUOJ]+$", re.IGNORECASE)

sequences = []

for idx, seq in enumerate(chains):
    chain_id = chr(ord("A") + idx)
    seq = seq.upper()

    if not aa_pattern.match(seq):
        raise SystemExit(
            f"ERROR: chain {chain_id} contains characters not supported by this protein-only generator: {seq}"
        )

    sequences.append({
        "protein": {
            "id": chain_id,
            "sequence": seq
        }
    })

# Parse model seeds
try:
    model_seeds = [
        int(seed.strip())
        for seed in seeds_raw.split(",")
        if seed.strip()
    ]
except ValueError:
    raise SystemExit("ERROR: MODEL_SEEDS must be a comma-separated list of integers.")

if not model_seeds:
    raise SystemExit("ERROR: MODEL_SEEDS is empty.")

data = {
    "name": case,
    "dialect": "alphafold3",
    "version": 1,
    "modelSeeds": model_seeds,
    "sequences": sequences
}

json_file.parent.mkdir(parents=True, exist_ok=True)

with json_file.open("w") as handle:
    json.dump(data, handle, indent=2)

print(f"AlphaFold 3 JSON written to: {json_file}")
print(f"Number of chains: {len(sequences)}")
print(f"Chain IDs: {', '.join(chr(ord('A') + i) for i in range(len(sequences)))}")
print(f"Number of model seeds: {len(model_seeds)}")
PY

# RUN ALPHAFOLD 3
echo "Running AlphaFold 3..."
echo "CASE:       ${CASE}"
echo "JSON_FILE:  ${JSON_FILE}"
echo "OUTPUT_DIR: ${OUTPUT_DIR}"
echo "LOG_FILE:   ${LOG_FILE}"

"${RUN_AF3}" \
    --jackhmmer_binary_path="${HMMER3_BINDIR}/jackhmmer" \
    --nhmmer_binary_path="${HMMER3_BINDIR}/nhmmer" \
    --hmmalign_binary_path="${HMMER3_BINDIR}/hmmalign" \
    --hmmsearch_binary_path="${HMMER3_BINDIR}/hmmsearch" \
    --hmmbuild_binary_path="${HMMER3_BINDIR}/hmmbuild" \
    --db_dir="${DB_DIR}" \
    --model_dir="${MODEL_DIR}" \
    --json_path="${JSON_FILE}" \
    --output_dir="${OUTPUT_DIR}" \
    --buckets="${BUCKETS}" \
    2>&1 | tee -a "${LOG_FILE}"

echo "AlphaFold 3 finished."
echo "Results directory: ${OUTPUT_DIR}"

# CONVERT CIF TO PDB
echo "Converting CIF to PDB files"
for cif in "${OUTPUT_DIR}/${CASE}"/seed-*/*cif; do
    [[ -e "${cif}" ]] || continue
    "${CIF_TO_PDB}" "${cif}" "${cif/.cif/.pdb}"
done
