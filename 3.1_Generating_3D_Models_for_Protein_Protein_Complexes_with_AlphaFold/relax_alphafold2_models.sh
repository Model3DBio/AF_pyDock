#!/bin/bash
set -euo pipefail

# REPOSITORY PATHS
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_STUDIES_DIR="${REPO_ROOT}/Case_Studies_Included"

# Optional case selector. If empty, all Case_Studies_Included/*/AF2 folders are used.
CASE="${CASE:-}"

# Override this to relax models from any custom AF2 output directory.
INPUT_DIR="${INPUT_DIR:-}"

# Base directory used when INPUT_DIR is not set.
INPUT_ROOT="${INPUT_ROOT:-${CASE_STUDIES_DIR}}"

if [[ -n "${INPUT_DIR}" && "${INPUT_DIR}" != /* ]]; then
    INPUT_DIR="${REPO_ROOT}/${INPUT_DIR}"
fi

if [[ "${INPUT_ROOT}" != /* ]]; then
    INPUT_ROOT="${REPO_ROOT}/${INPUT_ROOT}"
fi

# By default, relax unrelaxed AF2 PDB files and skip outputs that already exist.
PDB_PATTERN="${PDB_PATTERN:-*_unrelaxed_*.pdb}"
OVERWRITE="${OVERWRITE:-0}"

# Relaxation command configuration.
GREASY="${GREASY:-/path/to/software/greasy/}"
GREASY_NWORKERS="${GREASY_NWORKERS:-4}"
RUN_GREASY="${RUN_GREASY:-1}"
JOB_FILE="${JOB_FILE:-${SCRIPT_DIR}/${CASE:-case_studies}.af2_relax.greasy}"

if [[ "${JOB_FILE}" != /* ]]; then
    JOB_FILE="${SCRIPT_DIR}/${JOB_FILE}"
fi

input_dirs=()

if [[ -n "${INPUT_DIR}" ]]; then
    input_dirs=("${INPUT_DIR}")
elif [[ -n "${CASE}" ]]; then
    input_dirs=("${INPUT_ROOT}/${CASE}/AF2")
else
    shopt -s nullglob
    input_dirs=("${INPUT_ROOT}"/*/AF2)
    shopt -u nullglob
fi

if [[ ${#input_dirs[@]} -eq 0 ]]; then
    echo "ERROR: no AF2 input directories found." >&2
    echo "Set CASE, INPUT_DIR, or INPUT_ROOT." >&2
    exit 1
fi

: > "${JOB_FILE}"
job_count=0
skipped_count=0
valid_dir_count=0

for dir in "${input_dirs[@]}"; do
    if [[ ! -d "${dir}" ]]; then
        echo "WARNING: input directory does not exist, skipping: ${dir}" >&2
        continue
    fi

    valid_dir_count=$((valid_dir_count + 1))

    while IFS= read -r -d '' pdb_file; do
        relaxed_file="${pdb_file/_unrelaxed_/_relaxed_}"

        if [[ "${relaxed_file}" == "${pdb_file}" ]]; then
            skipped_count=$((skipped_count + 1))
            continue
        fi

        if [[ "${OVERWRITE}" != "1" && -e "${relaxed_file}" ]]; then
            skipped_count=$((skipped_count + 1))
            continue
        fi

        printf 'colabfold_relax --max-iterations 2000 --tolerance 2.39 --stiffness 10.0 --max-outer-iterations 3 --use-gpu "%s" "%s"\n' \
            "${pdb_file}" "${relaxed_file}" >> "${JOB_FILE}"
        job_count=$((job_count + 1))
    done < <(find "${dir}" -type f -name "${PDB_PATTERN}" -print0 | sort -z)
done

if [[ "${valid_dir_count}" -eq 0 ]]; then
    echo "ERROR: none of the selected AF2 input directories exists." >&2
    exit 1
fi

echo "Relaxation jobs written to: ${JOB_FILE}"
echo "Jobs: ${job_count}"
echo "Skipped existing or unmatched files: ${skipped_count}"

if [[ "${job_count}" -eq 0 ]]; then
    echo "No relaxation jobs to run."
    exit 0
fi

if [[ "${RUN_GREASY}" == "1" ]]; then
    export GREASY_NWORKERS
    "${GREASY}/greasy" "${JOB_FILE}"
else
    echo "RUN_GREASY=0, not launching Greasy."
fi
