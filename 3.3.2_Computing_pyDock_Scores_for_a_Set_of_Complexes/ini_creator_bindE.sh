#!/bin/bash
# Repository paths are derived from this script, not from the caller's cwd.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_STUDIES_DIR="${REPO_ROOT}/4_Case_Studies"

# Case ID
export CASE="${CASE:-4POU}"
export PYDOCK="${PYDOCK:-/usr/local/software/pyDock3/}"
export PYDOCK_BINARY="${PYDOCK_BINARY:-pyDock3}"
export GREASY_HOME="${GREASY_HOME:-/path/to/software/GREASY_2.2}"
export GREASY="${GREASY:-${GREASY_HOME%/}/bin}"
export GREASY_NWORKERS="${GREASY_NWORKERS:-8}"
export CHAINS_REC_LIG_VALUES=${CHAINS_REC_LIG_VALUES:-"A B"}
RUN_CAPRI_RMSD="${RUN_CAPRI_RMSD:-0}"
REFERENCE_PDB="${REFERENCE_PDB:-}"
REFERENCE_CHAINS_DECLARATION="$(declare -p REFERENCE_CHAINS_REC_LIG_VALUES 2>/dev/null || true)"
if [[ -z "${REFERENCE_CHAINS_DECLARATION}" ]]; then
    REFERENCE_CHAINS_REC_LIG_VALUES=()
fi
CIF_TO_PDB="${CIF_TO_PDB:-${REPO_ROOT}/3.1_Generating_3D_Models_for_Protein_Protein_Complexes_with_AlphaFold/cif_to_pdb.py}"

CASE_DIR="${CASE_STUDIES_DIR}/${CASE}"
AF2_DIR="${CASE_DIR}/AF2"
AF3_DIR="${CASE_DIR}/AF3"
AF3_LOCAL_DIR="${AF3_DIR}/output/${CASE}"
JOB_FILE="${CASE_DIR}/Greasy_BindE_mul_ligs.txt"

case "${RUN_CAPRI_RMSD}" in
    0|1) ;;
    *)
        echo "ERROR: RUN_CAPRI_RMSD must be 0 or 1." >&2
        exit 1
        ;;
esac

if [[ "$(declare -p CHAINS_REC_LIG_VALUES)" == declare\ -A* ||
      "${REFERENCE_CHAINS_DECLARATION}" == declare\ -A* ]]; then
    echo "ERROR: chain definitions must be a scalar or an indexed Bash array." >&2
    exit 1
fi

if [[ ! -d "${AF2_DIR}" && ! -d "${AF3_DIR}" ]]; then
    echo "ERROR: neither AF2 nor AF3 directory exists for case ${CASE}." >&2
    echo "Expected at least one of:" >&2
    echo "  ${AF2_DIR}" >&2
    echo "  ${AF3_DIR}" >&2
    exit 1
fi

if [[ "${RUN_CAPRI_RMSD}" == "0" ]]; then
    if [[ -n "${REFERENCE_PDB}" || ${#REFERENCE_CHAINS_REC_LIG_VALUES[@]} -gt 0 ]]; then
        echo "ERROR: reference variables require RUN_CAPRI_RMSD=1." >&2
        exit 1
    fi
else
    if [[ -z "${REFERENCE_PDB}" ]]; then
        echo "ERROR: REFERENCE_PDB is required when RUN_CAPRI_RMSD=1." >&2
        exit 1
    fi

    if [[ "${REFERENCE_PDB}" = /* ]]; then
        REFERENCE_PDB_PATH="${REFERENCE_PDB}"
    else
        REFERENCE_PDB_PATH="${CASE_DIR}/${REFERENCE_PDB}"
    fi

    if [[ ! -s "${REFERENCE_PDB_PATH}" ]]; then
        echo "ERROR: reference PDB does not exist or is empty:" >&2
        echo "  ${REFERENCE_PDB_PATH}" >&2
        exit 1
    fi

    REFERENCE_PDB_PATH="$(cd "$(dirname "${REFERENCE_PDB_PATH}")" && pwd)/$(basename "${REFERENCE_PDB_PATH}")"
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

REFERENCE_RECEPTOR_VALUES=()
REFERENCE_LIGAND_VALUES=()

if [[ "${RUN_CAPRI_RMSD}" == "1" ]]; then
    if [[ ${#REFERENCE_CHAINS_REC_LIG_VALUES[@]} -ne ${#RECEPTOR_VALUES[@]} ]]; then
        echo "ERROR: REFERENCE_CHAINS_REC_LIG_VALUES must contain one entry for each" >&2
        echo "CHAINS_REC_LIG_VALUES entry when RUN_CAPRI_RMSD=1." >&2
        exit 1
    fi

    for value in "${REFERENCE_CHAINS_REC_LIG_VALUES[@]}"; do
        read -r REF_REC REF_LIG extra <<< "${value}"

        if [[ -z "${REF_REC}" || -z "${REF_LIG}" || -n "${extra}" ]]; then
            echo "ERROR: invalid reference receptor/ligand definition: ${value}" >&2
            echo 'Expected: "REFERENCE_RECEPTOR_CHAINS REFERENCE_LIGAND_CHAINS"' >&2
            exit 1
        fi

        REFERENCE_RECEPTOR_VALUES+=("${REF_REC}")
        REFERENCE_LIGAND_VALUES+=("${REF_LIG}")
    done

    for config_index in "${!RECEPTOR_VALUES[@]}"; do
        model_rec_commas="${RECEPTOR_VALUES[${config_index}]//[^,]/}"
        model_lig_commas="${LIGAND_VALUES[${config_index}]//[^,]/}"
        reference_rec_commas="${REFERENCE_RECEPTOR_VALUES[${config_index}]//[^,]/}"
        reference_lig_commas="${REFERENCE_LIGAND_VALUES[${config_index}]//[^,]/}"

        if [[ ${#model_rec_commas} -ne ${#reference_rec_commas} ||
              ${#model_lig_commas} -ne ${#reference_lig_commas} ]]; then
            echo "ERROR: model and reference chain counts do not match for configuration" >&2
            echo "  model:     ${RECEPTOR_VALUES[${config_index}]} ${LIGAND_VALUES[${config_index}]}" >&2
            echo "  reference: ${REFERENCE_RECEPTOR_VALUES[${config_index}]} ${REFERENCE_LIGAND_VALUES[${config_index}]}" >&2
            exit 1
        fi
    done
fi

# AlphaFold Server models are CIF files directly below AF3/fold*/.
for cif_file in "${AF3_DIR}"/fold*/*_model_*.cif; do
    [[ -f "${cif_file}" && "${cif_file}" =~ _model_[0-9]+\.cif$ ]] || continue
    pdb_file="${cif_file%.cif}.pdb"
    [[ -s "${pdb_file}" ]] && continue

    echo "Converting AlphaFold Server model: ${cif_file}"
    if ! "${CIF_TO_PDB}" "${cif_file}" "${pdb_file}" || [[ ! -s "${pdb_file}" ]]; then
        echo "ERROR: CIF-to-PDB conversion failed: ${cif_file}" >&2
        rm -f -- "${pdb_file}"
        exit 1
    fi
done

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
        if [[ "${RUN_CAPRI_RMSD}" == "1" && "${pdb_file}" == "${REFERENCE_PDB_PATH}" ]]; then
            continue
        fi

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

            if [[ "${RUN_CAPRI_RMSD}" == "1" ]]; then
                REF_REC="${REFERENCE_RECEPTOR_VALUES[${config_index}]}"
                REF_LIG="${REFERENCE_LIGAND_VALUES[${config_index}]}"

                cat <<EOF >> "${ini_path}"

[reference]
pdb = ${REFERENCE_PDB_PATH}
recmol = $REF_REC
ligmol = $REF_LIG
newrecmol = $REC
newligmol = $LIG
EOF

                printf 'cd "%s" && timeout 5m "%s/%s" "%s" bindEy && timeout 5m "%s/%s" "%s" capriRMSD\n' \
                    "${h%/}" "${PYDOCK%/}" "${PYDOCK_BINARY}" \
                    "${ini_file%.ini}" "${PYDOCK%/}" "${PYDOCK_BINARY}" \
                    "${ini_file%.ini}" >> "${JOB_FILE}"
            else
                printf 'cd "%s" && timeout 5m "%s/%s" "%s" bindEy\n' \
                    "${h%/}" "${PYDOCK%/}" "${PYDOCK_BINARY}" \
                    "${ini_file%.ini}" >> "${JOB_FILE}"
            fi
        done
    done
done
"${GREASY%/}/greasy" "${JOB_FILE}"
