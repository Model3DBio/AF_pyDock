#!/bin/bash
# Repository paths are derived from this script, not from the caller's cwd.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_STUDIES_DIR="${REPO_ROOT}/Case_Studies_Included"

# Case ID
export CASE="${CASE:-4POU}"
export PYDOCK="${PYDOCK:-/usr/local/software/pyDock3/}"
export GREASY="${GREASY:-/path/to/software/greasy/}"
export GREASY_NWORKERS=${GREASY_NWORKERS:-8}
export CHAINS_REC_LIG_VALUES=${CHAINS_REC_LIG_VALUES:-"A B"}

CASE_DIR="${CASE_STUDIES_DIR}/${CASE}"
AF2_DIR="${CASE_DIR}/AF2"
AF3_DIR="${CASE_DIR}/AF3"
AF3_LOCAL_DIR="${AF3_DIR}/output/${CASE}"
JOB_FILE="${CASE_DIR}/Greasy_BindE_mul_ligs.txt"

if [[ ! -d "${AF2_DIR}" && ! -d "${AF3_DIR}" ]]; then
    echo "ERROR: neither AF2 nor AF3 directory exists for case ${CASE}." >&2
    echo "Expected at least one of:" >&2
    echo "  ${AF2_DIR}" >&2
    echo "  ${AF3_DIR}" >&2
    exit 1
fi

shopt -s nullglob
model_dirs=(
    "${AF2_DIR}"/*cola*v*/
    "${AF3_DIR}"/fold*/
    "${AF3_LOCAL_DIR}"/seed-*/
)
shopt -u nullglob

if [[ ${#model_dirs[@]} -eq 0 ]]; then
    echo "ERROR: no AF2 or AF3 model directories found for case ${CASE}." >&2
    exit 1
fi

: > "${JOB_FILE}"

for h in "${model_dirs[@]}"; do
    echo "${h}"

    shopt -s nullglob
    pdb_files=(
        "${h}"/*[0-9].pdb
        "${h}"/*_sample-[0-9]*_model.pdb
    )
    shopt -u nullglob

    for pdb_file in "${pdb_files[@]}"; do
        j="$(basename "${pdb_file}")"
        echo "$j"

        for value in "${CHAINS_REC_LIG_VALUES[@]}"; do
            REC=$(echo "$value" | cut -d" " -f1)
            LIG=$(echo "$value" | cut -d" " -f2)
            CH=${LIG/,}

            ini_file="${j/.pdb}_LIG_${CH}.ini"
            ini_path="${h}/${ini_file}"

            cat <<EOF > "${ini_path}"
[receptor]
pdb = $j
mol = $REC
newmol = $REC

[ligand]
pdb = $j
mol = $LIG
newmol = $LIG
EOF

            printf 'cd "%s"; timeout 5m "%s/pyDock3" "%s" bindEy\n' \
                "${h%/}" "${PYDOCK%/}" "${ini_file/.ini}" >> "${JOB_FILE}"
        done
    done
done
"${GREASY%/}/greasy" "${JOB_FILE}"
