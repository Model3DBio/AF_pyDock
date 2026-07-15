#!/bin/bash
# Repository paths are derived from this script, not from the caller's cwd.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_STUDIES_DIR="${REPO_ROOT}/Case_Studies_Included"

# Case ID
export CASE="${CASE:-4POU}"
export PYDOCK="${PYDOCK:-/usr/local/software/pyDock3/}"
export GREASY_HOME="${GREASY_HOME:-/path/to/software/GREASY_2.2}"
export GREASY="${GREASY:-${GREASY_HOME%/}/bin}"
export GREASY_NWORKERS="${GREASY_NWORKERS:-8}"
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

RECEPTOR_VALUES=()
LIGAND_VALUES=()
RECEPTOR_TAGS=()
LIGAND_TAGS=()
declare -A SEEN_CONFIGURATIONS=()

for value in "${CHAINS_REC_LIG_VALUES[@]}"; do
    read -r REC LIG extra <<< "${value}"

    if [[ -z "${REC}" || -z "${LIG}" || -n "${extra}" ]]; then
        echo "ERROR: invalid receptor/ligand definition: ${value}" >&2
        echo 'Expected: "RECEPTOR_CHAINS LIGAND_CHAINS"' >&2
        exit 1
    fi

    REC_TAG="${REC//,/-}"
    LIG_TAG="${LIG//,/-}"
    configuration_key="${REC_TAG}|${LIG_TAG}"

    if [[ ${SEEN_CONFIGURATIONS["${configuration_key}"]+_} ]]; then
        echo "ERROR: duplicate receptor/ligand configuration after filename normalization:" >&2
        echo "  ${SEEN_CONFIGURATIONS["${configuration_key}"]}" >&2
        echo "  ${value}" >&2
        exit 1
    fi

    SEEN_CONFIGURATIONS["${configuration_key}"]="${value}"
    RECEPTOR_VALUES+=("${REC}")
    LIGAND_VALUES+=("${LIG}")
    RECEPTOR_TAGS+=("${REC_TAG}")
    LIGAND_TAGS+=("${LIG_TAG}")
done

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

        for config_index in "${!RECEPTOR_VALUES[@]}"; do
            REC="${RECEPTOR_VALUES[${config_index}]}"
            LIG="${LIGAND_VALUES[${config_index}]}"
            REC_TAG="${RECEPTOR_TAGS[${config_index}]}"
            LIG_TAG="${LIGAND_TAGS[${config_index}]}"

            ini_file="${j%.pdb}_LIG_${LIG_TAG}_REC_${REC_TAG}.ini"
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
