#!/bin/bash

# REPOSITORY PATHS
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_STUDIES_DIR="${REPO_ROOT}/Case_Studies_Included"

# SOFTWARE CHECKS
AF2_CONDA_ENV="${AF2_CONDA_ENV:-Alphafold2}"
COLABFOLD_BATCH="${COLABFOLD_BATCH:-colabfold_batch}"

ACTIVE_CONDA_ENV="${CONDA_DEFAULT_ENV:-}"

if COLABFOLD_BATCH_PATH="$(command -v "${COLABFOLD_BATCH}" 2>/dev/null)"; then
    if [[ "${ACTIVE_CONDA_ENV}" == "${AF2_CONDA_ENV}" ]]; then
        echo "INFO: expected Conda environment is active: ${AF2_CONDA_ENV}"
    else
        echo "WARNING: Conda environment '${AF2_CONDA_ENV}' is not active (current: ${ACTIVE_CONDA_ENV:-none})." >&2
        echo "WARNING: continuing because ${COLABFOLD_BATCH} was found at: ${COLABFOLD_BATCH_PATH}" >&2
    fi
else
    if [[ "${ACTIVE_CONDA_ENV}" == "${AF2_CONDA_ENV}" ]]; then
        echo "ERROR: Conda environment '${AF2_CONDA_ENV}' is active, but ${COLABFOLD_BATCH} could not be resolved as an executable." >&2
        echo "The ColabFold installation in this environment may be incomplete." >&2
    else
        echo "ERROR: Conda environment '${AF2_CONDA_ENV}' is not active and ${COLABFOLD_BATCH} could not be resolved as an executable." >&2
        echo "Activate it with: conda activate ${AF2_CONDA_ENV}" >&2
        echo "Alternatively, install ColabFold and expose ${COLABFOLD_BATCH} through PATH or set COLABFOLD_BATCH explicitly." >&2
    fi
    exit 1
fi

# CASE INPUTS
# Case ID
export CASE="${CASE:-4POU}"
#export CASE="${CASE:-2FJG}"

# Protein sequences separated by ":".
# Each ":"-separated block will be treated as a different chain.
export SEQ="${SEQ:-QVQLVESGGGLVQAGGSLRLSCAASGYPHPYLHMGWFRQAPGKEREGVAAMDSGGGGTLYADSVKGRFTISRDKGKNTVYLQMDSLKPEDTATYYCAAGGYQLRDRTYGHWGQGTQVTVSS:KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQSYSTMSITDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV}"
#export SEQ="${SEQ:-EVQLVESGGGLVQPGGSLRLSCAASGFTISDYWIHWVRQAPGKGLEWVAGITPAGGYTYYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARFVFFLPYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSC:DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNR:EVVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIGEMSFLQHNKCECRPK:EVVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIGEMSFLQHNKCECRPK}"

CASE_DIR="${CASE_STUDIES_DIR}/${CASE}/AF2"
FASTA_FILE="${CASE_DIR}/${CASE}.fasta"

mkdir -p "${CASE_DIR}"

# Create FASTA file (60 chars per line)
{
    echo ">${CASE}"
    echo "${SEQ}" | fold -w 60
} > "${FASTA_FILE}"

# Run AF2-Multimer
for h in v1 v2 v3; do
    "${COLABFOLD_BATCH_PATH}" "${FASTA_FILE}" "${FASTA_FILE%.fasta}.colabfold.${h}" \
        --save-recycles \
        --model-type "alphafold2_multimer_${h}" \
        --recycle-early-stop-tolerance 0 \
        --rank multimer \
        --use-dropout \
        --num-recycle 20
done
