#!/usr/bin/env bash
set -euo pipefail

# Adapted from Yoshitaka Moriwaki's LocalColabFold Linux installer:
# https://github.com/YoshitakaMo/localcolabfold
# Original work copyright (c) 2021 Yoshitaka Moriwaki and licensed under MIT.
# The upstream license notice is preserved in LOCALCOLABFOLD_LICENSE.

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ENV_NAME="Alphafold2"
CURRENT_PATH="$(pwd)"
COLABFOLD_DIR="${CURRENT_PATH}/localcolabfold"
LOCAL_CONDA_DIR="${COLABFOLD_DIR}/conda"

# ---------------------------------------------------------------------------
# Basic requirements
# ---------------------------------------------------------------------------

command -v wget >/dev/null 2>&1 || {
    echo "ERROR: wget is not installed."
    echo "Install it using apt, dnf, yum, or the package manager of your system."
    exit 1
}

mkdir -p "${COLABFOLD_DIR}"

# ---------------------------------------------------------------------------
# Initialize an existing Conda installation or install Miniforge locally
# ---------------------------------------------------------------------------

initialize_conda() {
    local conda_base

    # Conda may be available as a shell function.
    if type conda >/dev/null 2>&1; then
        conda_base="$(conda info --base 2>/dev/null || true)"

        if [[ -n "${conda_base}" && -f "${conda_base}/etc/profile.d/conda.sh" ]]; then
            # shellcheck disable=SC1091
            source "${conda_base}/etc/profile.d/conda.sh"
            return 0
        fi
    fi

    # Conda may exist as an executable without having been initialized
    # in the current shell.
    if command -v conda >/dev/null 2>&1; then
        conda_base="$(conda info --base 2>/dev/null || true)"

        if [[ -n "${conda_base}" && -f "${conda_base}/etc/profile.d/conda.sh" ]]; then
            # shellcheck disable=SC1091
            source "${conda_base}/etc/profile.d/conda.sh"
            return 0
        fi
    fi

    # Search common Conda installation locations.
    local candidate

    for candidate in \
        "${HOME}/miniforge3" \
        "${HOME}/mambaforge" \
        "${HOME}/miniconda3" \
        "${HOME}/anaconda3" \
        "/opt/conda" \
        "/usr/local/miniforge3" \
        "/usr/local/miniconda3"
    do
        if [[ -f "${candidate}/etc/profile.d/conda.sh" ]]; then
            # shellcheck disable=SC1091
            source "${candidate}/etc/profile.d/conda.sh"
            return 0
        fi
    done

    return 1
}

if initialize_conda; then
    CONDA_BASE="$(conda info --base)"
    echo "Existing Conda installation detected:"
    echo "  ${CONDA_BASE}"
else
    echo "No functional Conda installation was detected."
    echo "Installing Miniforge locally in:"
    echo "  ${LOCAL_CONDA_DIR}"

    MINIFORGE_INSTALLER="${COLABFOLD_DIR}/Miniforge3-Linux-x86_64.sh"

    wget -q \
        -O "${MINIFORGE_INSTALLER}" \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"

    bash "${MINIFORGE_INSTALLER}" \
        -b \
        -p "${LOCAL_CONDA_DIR}"

    rm -f "${MINIFORGE_INSTALLER}"

    # shellcheck disable=SC1091
    source "${LOCAL_CONDA_DIR}/etc/profile.d/conda.sh"

    CONDA_BASE="${LOCAL_CONDA_DIR}"

    echo "Miniforge installation completed."
fi

# Ensure that the selected Conda installation is available.
export PATH="${CONDA_BASE}/condabin:${PATH}"

echo "Conda executable:"
echo "  $(command -v conda)"
echo "Conda base:"
echo "  $(conda info --base)"

# ---------------------------------------------------------------------------
# Create or reuse the Alphafold2 environment
# ---------------------------------------------------------------------------

if conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
    echo "Conda environment '${ENV_NAME}' already exists."
    echo "Reusing the existing environment."
else
    echo "Creating Conda environment '${ENV_NAME}'..."

    conda create \
        --name "${ENV_NAME}" \
        --channel conda-forge \
        --channel bioconda \
        --strict-channel-priority \
        git \
        python=3.10 \
        openmm=8.2.0 \
        pdbfixer \
        kalign2=2.04 \
        hhsuite=3.3.0 \
        mmseqs2 \
        pip \
        -y
fi

conda activate "${ENV_NAME}"

ENV_PREFIX="${CONDA_PREFIX}"
ENV_PYTHON="${ENV_PREFIX}/bin/python"
ENV_COLABFOLD_BATCH="${ENV_PREFIX}/bin/colabfold_batch"

echo "Active environment:"
echo "  ${CONDA_DEFAULT_ENV}"
echo "Environment location:"
echo "  ${ENV_PREFIX}"

# ---------------------------------------------------------------------------
# Isolate the environment from packages installed under ~/.local
# ---------------------------------------------------------------------------

conda env config vars set \
    -n "${ENV_NAME}" \
    PYTHONUSERBASE=intentionally-disabled

# conda env config vars are normally applied on the next activation.
# Export the variable now so that it also affects this installation process.
export PYTHONUSERBASE="intentionally-disabled"

# Also prevent Python from loading the user site-packages directory.
conda env config vars set \
    -n "${ENV_NAME}" \
    PYTHONNOUSERSITE=1

export PYTHONNOUSERSITE=1

# Reactivate to ensure that the stored environment variables are loaded.
conda deactivate
conda activate "${ENV_NAME}"

ENV_PREFIX="${CONDA_PREFIX}"
ENV_PYTHON="${ENV_PREFIX}/bin/python"
ENV_COLABFOLD_BATCH="${ENV_PREFIX}/bin/colabfold_batch"

# ---------------------------------------------------------------------------
# Install ColabFold, AlphaFold dependencies, JAX and TensorFlow
# ---------------------------------------------------------------------------

"${ENV_PYTHON}" -m pip install --upgrade pip setuptools wheel

"${ENV_PYTHON}" -m pip install \
    --no-warn-conflicts \
    "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"

"${ENV_PYTHON}" -m pip install \
    "colabfold[alphafold]"

"${ENV_PYTHON}" -m pip install \
    --upgrade \
    "jax[cuda12]==0.5.3"

"${ENV_PYTHON}" -m pip install \
    --upgrade \
    tensorflow \
    silence_tensorflow

# ---------------------------------------------------------------------------
# Download the LocalColabFold updater
# ---------------------------------------------------------------------------

wget -qnc \
    -O "${COLABFOLD_DIR}/update_linux.sh" \
    "https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/update_linux.sh"

chmod +x "${COLABFOLD_DIR}/update_linux.sh"

# ---------------------------------------------------------------------------
# Locate and modify the installed ColabFold package
# ---------------------------------------------------------------------------

COLABFOLD_PACKAGE_DIR="$(
    "${ENV_PYTHON}" -c \
        'import pathlib, colabfold; print(pathlib.Path(colabfold.__file__).resolve().parent)'
)"

if [[ ! -d "${COLABFOLD_PACKAGE_DIR}" ]]; then
    echo "ERROR: ColabFold package directory was not found."
    exit 1
fi

echo "ColabFold package directory:"
echo "  ${COLABFOLD_PACKAGE_DIR}"

pushd "${COLABFOLD_PACKAGE_DIR}" >/dev/null

# Store downloaded parameters and related data in localcolabfold/colabfold.
if grep -qF 'appdirs.user_cache_dir(__package__ or "colabfold")' download.py; then
    sed -i \
        "s#appdirs.user_cache_dir(__package__ or \"colabfold\")#\"${COLABFOLD_DIR}/colabfold\"#g" \
        download.py
else
    echo "WARNING: The expected cache-directory expression was not found in download.py."
    echo "The ColabFold source may have changed, so download.py was not modified."
fi

# Suppress TensorFlow informational warnings.
if ! grep -qF "from silence_tensorflow import silence_tensorflow" batch.py; then
    if grep -qF "from io import StringIO" batch.py; then
        sed -i \
            '/from io import StringIO/a from silence_tensorflow import silence_tensorflow\nsilence_tensorflow()' \
            batch.py
    else
        echo "WARNING: The expected import was not found in batch.py."
        echo "TensorFlow warning suppression was not inserted."
    fi
fi

rm -rf __pycache__

popd >/dev/null

# ---------------------------------------------------------------------------
# Download AlphaFold2 weights
# ---------------------------------------------------------------------------

"${ENV_PYTHON}" -m colabfold.download

# ---------------------------------------------------------------------------
# Final information
# ---------------------------------------------------------------------------

echo
echo "Download of AlphaFold2 weights finished."
echo "------------------------------------------------------------"
echo "Installation of ColabFold finished."
echo
echo "Conda environment:"
echo "  ${ENV_NAME}"
echo
echo "Environment location:"
echo "  ${ENV_PREFIX}"
echo
echo "To use ColabFold:"
echo
echo "  conda activate ${ENV_NAME}"
echo "  colabfold_batch --help"
echo
echo "The environment has the following isolation variables:"
echo "  PYTHONUSERBASE=intentionally-disabled"
echo "  PYTHONNOUSERSITE=1"
echo
echo "ColabFold executable:"
echo "  ${ENV_COLABFOLD_BATCH}"
