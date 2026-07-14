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
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_STUDIES_DIR="${REPO_ROOT}/Case_Studies_Included"
APPDIR="${APPDIR:-/home/user/Programs}"
ALPHAFOLD3DIR="${ALPHAFOLD3DIR:-${APPDIR}/alphafold3}"
CIF_TO_PDB="${CIF_TO_PDB:-${SCRIPT_DIR}/cif_to_pdb.py}"
AF3_CONDA_ENV="${AF3_CONDA_ENV:-Alphafold3}"
# HMMER binaries.
# If HMMER3_BINDIR is not set, use the active conda environment or PATH.
# To use system binaries, run for example:
#   HMMER3_BINDIR=/usr/bin ./run_alphafold3.sh
HMMER3_BINDIR="${HMMER3_BINDIR:-}"

DB_DIR="${DB_DIR:-${ALPHAFOLD3DIR}/public_databases}"
MODEL_DIR="${MODEL_DIR:-${ALPHAFOLD3DIR}/models}"

# AlphaFold 3 runner
RUN_AF3="${RUN_AF3:-run_alphafold.py}"

# AlphaFold 3 buckets
BUCKETS="${BUCKETS:-256,512,768,1024,1280,1536,2048,2560,3072,3584,4096,4608,5120}"

# WORKING DIRECTORIES
WORK_DIR="${CASE_STUDIES_DIR}"
CASE_DIR="${WORK_DIR}/${CASE}/AF3"
JSON_FILE="${CASE_DIR}/${CASE}.json"
OUTPUT_DIR="${CASE_DIR}/output"
LOG_FILE="${OUTPUT_DIR}/af3_run.log"

# BASIC CHECKS
if [[ -z "${SEQ}" ]]; then
    echo "ERROR: SEQ is empty." >&2
    exit 1
fi

ACTIVE_CONDA_ENV="${CONDA_DEFAULT_ENV:-}"

if RUN_AF3_PATH="$(command -v "${RUN_AF3}" 2>/dev/null)"; then
    if [[ "${ACTIVE_CONDA_ENV}" == "${AF3_CONDA_ENV}" ]]; then
        echo "INFO: expected Conda environment is active: ${AF3_CONDA_ENV}"
    else
        echo "WARNING: Conda environment '${AF3_CONDA_ENV}' is not active (current: ${ACTIVE_CONDA_ENV:-none})." >&2
        echo "WARNING: continuing because ${RUN_AF3} was found at: ${RUN_AF3_PATH}" >&2
    fi
else
    if [[ "${ACTIVE_CONDA_ENV}" == "${AF3_CONDA_ENV}" ]]; then
        echo "ERROR: Conda environment '${AF3_CONDA_ENV}' is active, but ${RUN_AF3} could not be resolved as an executable." >&2
        echo "The AlphaFold 3 installation in this environment may be incomplete." >&2
    else
        echo "ERROR: Conda environment '${AF3_CONDA_ENV}' is not active and ${RUN_AF3} could not be resolved as an executable." >&2
        echo "Activate it with: conda activate ${AF3_CONDA_ENV}" >&2
        echo "Alternatively, install AlphaFold 3 and expose ${RUN_AF3} through PATH or set RUN_AF3 explicitly." >&2
    fi
    exit 1
fi

declare -A HMMER_BIN_PATHS=()
HMMER_BINS=(jackhmmer nhmmer hmmalign hmmsearch hmmbuild)

if [[ -n "${HMMER3_BINDIR}" ]]; then
    if [[ ! -d "${HMMER3_BINDIR}" ]]; then
        echo "ERROR: HMMER3_BINDIR does not exist: ${HMMER3_BINDIR}" >&2
        exit 1
    fi

    for bin in "${HMMER_BINS[@]}"; do
        HMMER_BIN_PATHS["${bin}"]="${HMMER3_BINDIR}/${bin}"
    done
elif [[ -n "${CONDA_PREFIX:-}" ]]; then
    HMMER3_BINDIR="${CONDA_PREFIX}/bin"

    for bin in "${HMMER_BINS[@]}"; do
        HMMER_BIN_PATHS["${bin}"]="${HMMER3_BINDIR}/${bin}"
    done
else
    for bin in "${HMMER_BINS[@]}"; do
        resolved_bin="$(command -v "${bin}" 2>/dev/null || true)"

        if [[ -z "${resolved_bin}" ]]; then
            echo "ERROR: Conda environment '${AF3_CONDA_ENV}' is not active and ${bin} was not found in PATH." >&2
            echo "Activate the environment, define HMMER3_BINDIR, or expose the HMMER executables through PATH." >&2
            exit 1
        fi

        HMMER_BIN_PATHS["${bin}"]="${resolved_bin}"
    done
fi

for bin in "${HMMER_BINS[@]}"; do
    if [[ ! -x "${HMMER_BIN_PATHS[${bin}]}" ]]; then
        echo "ERROR: could not find executable ${bin} at ${HMMER_BIN_PATHS[${bin}]}" >&2
        if [[ -n "${CONDA_PREFIX:-}" && "${HMMER3_BINDIR}" == "${CONDA_PREFIX}/bin" ]]; then
            echo "The active Conda environment may have an incomplete HMMER installation." >&2
        else
            echo "Define HMMER3_BINDIR or activate an environment containing all required HMMER executables." >&2
        fi
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

mkdir -p "${CASE_DIR}"
mkdir -p "${OUTPUT_DIR}"

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

"${RUN_AF3_PATH}" \
    --jackhmmer_binary_path="${HMMER_BIN_PATHS[jackhmmer]}" \
    --nhmmer_binary_path="${HMMER_BIN_PATHS[nhmmer]}" \
    --hmmalign_binary_path="${HMMER_BIN_PATHS[hmmalign]}" \
    --hmmsearch_binary_path="${HMMER_BIN_PATHS[hmmsearch]}" \
    --hmmbuild_binary_path="${HMMER_BIN_PATHS[hmmbuild]}" \
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
