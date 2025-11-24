#!/bin/bash

# Case ID
export CASE="${CASE:-4POU}"

# Sequence
export SEQ="QVQLVESGGGLVQAGGSLRLSCAASGYPHPYLHMGWFRQAPGKEREGVAAMDSGGGGTLYADSVKGRFTISRDKGKNTVYLQMDSLKPEDTATYYCAAGGYQLRDRTYGHWGQGTQVTVSS:KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQSYSTMSITDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV"

# Greasy executable path
export GREASY="${GREASY:-/path/to/software/greasy}"

mkdir -p "${CASE}/AF2/"

# Create FASTA file (60 chars per line)
{
    echo ">${CASE}"
    echo "${SEQ}" | fold -w 60
} > "${CASE}/AF2/${CASE}.fasta"

# Run AF2-Multimer
for h in v1 v2 v3; do
    for i in */AF2/*fasta; do
        colabfold_batch "$(pwd)/$i" "$(pwd)/${i/.fasta/.colabfold.$h}" \
            --save-recycles \
            --model-type "alphafold2_multimer_${h}" \
            --re-cycle-early-stop-tolerance 0 \
            --rank multimer \
            --use-dropout \
            --num-recycle 20 \
            --amber \
            --use-gpu-relax
    done
done

# Generate relaxation job list
for f in $(find -name "*_*\.r*[0-9].pdb"); do
    echo "colabfold_relax --max-iterations 2000 --tolerance 2.39 \
--stiffness 10.0 --max-outer-iterations 3 --use-gpu $f ${f/unrelaxed/relaxed}"
done > "${CASE}.colabfold.relax.greasy"

# Number of simultaneous GPU relaxation jobs
export GREASY_NWORKERS=4

# Run Greasy
"$GREASY/greasy" "${CASE}.colabfold.relax.greasy"

