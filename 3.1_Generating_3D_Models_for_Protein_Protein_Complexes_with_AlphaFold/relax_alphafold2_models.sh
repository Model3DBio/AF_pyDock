#!/bin/bash
set -euo pipefail

# REPOSITORY PATHS
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_STUDIES_DIR="${REPO_ROOT}/4_Case_Studies"

# Optional case selector. If empty, all 4_Case_Studies/*/AF2 folders are used.
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
AF2_RELAX_VERSIONS="${AF2_RELAX_VERSIONS-v2 v3}"

AF2_RELAX_VERSION_LIST=()
read -r -a AF2_RELAX_VERSION_LIST <<< "${AF2_RELAX_VERSIONS}" || true

if [[ ${#AF2_RELAX_VERSION_LIST[@]} -eq 0 ]]; then
    echo "ERROR: AF2_RELAX_VERSIONS must contain at least one version." >&2
    exit 1
fi

SEEN_AF2_RELAX_VERSIONS=""
for version in "${AF2_RELAX_VERSION_LIST[@]}"; do
    case "${version}" in
        v1|v2|v3)
            ;;
        *)
            echo "ERROR: unsupported AlphaFold2-Multimer relaxation version: ${version}" >&2
            echo 'Set AF2_RELAX_VERSIONS to a space-separated selection of: v1 v2 v3' >&2
            exit 1
            ;;
    esac

    if [[ " ${SEEN_AF2_RELAX_VERSIONS} " == *" ${version} "* ]]; then
        echo "ERROR: duplicate AlphaFold2-Multimer relaxation version: ${version}" >&2
        exit 1
    fi
    SEEN_AF2_RELAX_VERSIONS+=" ${version}"
done

# Relaxation command configuration.
AF2_CONDA_ENV="${AF2_CONDA_ENV:-Alphafold2}"
COLABFOLD_RELAX="${COLABFOLD_RELAX:-colabfold_relax}"
GREASY_HOME="${GREASY_HOME:-/path/to/software/GREASY_2.2}"
GREASY="${GREASY:-${GREASY_HOME%/}/bin}"
GREASY_NWORKERS="${GREASY_NWORKERS:-4}"
RELAX_TIMEOUT="${RELAX_TIMEOUT:-8m}"
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
version_skipped_count=0
valid_dir_count=0

for dir in "${input_dirs[@]}"; do
    if [[ ! -d "${dir}" ]]; then
        echo "WARNING: input directory does not exist, skipping: ${dir}" >&2
        continue
    fi

    valid_dir_count=$((valid_dir_count + 1))

    while IFS= read -r -d '' pdb_file; do
        model_version=""
        case "${pdb_file}" in
            *.colabfold.v1/*) model_version="v1" ;;
            *.colabfold.v2/*) model_version="v2" ;;
            *.colabfold.v3/*) model_version="v3" ;;
        esac

        if [[ -n "${model_version}" && " ${AF2_RELAX_VERSION_LIST[*]} " != *" ${model_version} "* ]]; then
            version_skipped_count=$((version_skipped_count + 1))
            continue
        fi

        relaxed_file="${pdb_file/_unrelaxed_/_relaxed_}"

        if [[ "${relaxed_file}" == "${pdb_file}" ]]; then
            skipped_count=$((skipped_count + 1))
            continue
        fi

        if [[ "${OVERWRITE}" != "1" && -e "${relaxed_file}" ]]; then
            skipped_count=$((skipped_count + 1))
            continue
        fi

        printf 'timeout %q %q --max-iterations 2000 --tolerance 2.39 --stiffness 10.0 --max-outer-iterations 3 --use-gpu "%s" "%s"\n' \
            "${RELAX_TIMEOUT}" "${COLABFOLD_RELAX}" "${pdb_file}" "${relaxed_file}" >> "${JOB_FILE}"
        job_count=$((job_count + 1))
    done < <(find "${dir}" -type f -name "${PDB_PATTERN}" -print0 | sort -z)
done

if [[ "${valid_dir_count}" -eq 0 ]]; then
    echo "ERROR: none of the selected AF2 input directories exists." >&2
    exit 1
fi

echo "Relaxation jobs written to: ${JOB_FILE}"
echo "Jobs: ${job_count}"
echo "Skipped models from unselected AF2 versions: ${version_skipped_count}"
echo "Skipped existing or unmatched files: ${skipped_count}"

if [[ "${job_count}" -eq 0 ]]; then
    echo "No relaxation jobs to run."
    exit 0
fi

ACTIVE_CONDA_ENV="${CONDA_DEFAULT_ENV:-}"

if COLABFOLD_RELAX_PATH="$(command -v "${COLABFOLD_RELAX}" 2>/dev/null)"; then
    if [[ "${ACTIVE_CONDA_ENV}" == "${AF2_CONDA_ENV}" ]]; then
        echo "INFO: expected Conda environment is active: ${AF2_CONDA_ENV}"
    else
        echo "WARNING: Conda environment '${AF2_CONDA_ENV}' is not active (current: ${ACTIVE_CONDA_ENV:-none})." >&2
        echo "WARNING: ${COLABFOLD_RELAX} is available at: ${COLABFOLD_RELAX_PATH}" >&2
    fi
else
    if [[ "${RUN_GREASY}" == "1" ]]; then
        if [[ "${ACTIVE_CONDA_ENV}" == "${AF2_CONDA_ENV}" ]]; then
            echo "ERROR: Conda environment '${AF2_CONDA_ENV}' is active, but ${COLABFOLD_RELAX} could not be resolved as an executable." >&2
            echo "The ColabFold installation in this environment may be incomplete." >&2
        else
            echo "ERROR: Conda environment '${AF2_CONDA_ENV}' is not active and ${COLABFOLD_RELAX} could not be resolved as an executable." >&2
            echo "Activate it with: conda activate ${AF2_CONDA_ENV}" >&2
            echo "Alternatively, install ColabFold and expose ${COLABFOLD_RELAX} through PATH or set COLABFOLD_RELAX explicitly." >&2
        fi
        exit 1
    fi

    echo "WARNING: ${COLABFOLD_RELAX} could not be resolved as an executable." >&2
    echo "WARNING: the task file was created because RUN_GREASY=0, but it cannot be executed in the current environment." >&2
fi

if [[ "${RUN_GREASY}" == "1" ]]; then
    GREASY_EXECUTABLE="${GREASY%/}/greasy"

    if [[ ! -x "${GREASY_EXECUTABLE}" ]]; then
        echo "ERROR: Greasy executable not found or not executable: ${GREASY_EXECUTABLE}" >&2
        echo "Set GREASY to the directory containing the greasy executable." >&2
        exit 1
    fi

    export GREASY_NWORKERS
    "${GREASY_EXECUTABLE}" "${JOB_FILE}"
else
    echo "RUN_GREASY=0, not launching Greasy."
fi
