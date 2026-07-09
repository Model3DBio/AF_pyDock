#!/bin/bash

# Case ID
export CASE="${CASE:-4POU}"

# Sequence
export SEQ="${SEQ:-QVQLVESGGGLVQAGGSLRLSCAASGYPHPYLHMGWFRQAPGKEREGVAAMDSGGGGTLYADSVKGRFTISRDKGKNTVYLQMDSLKPEDTATYYCAAGGYQLRDRTYGHWGQGTQVTVSS:KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQSYSTMSITDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV}"

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
            --num-recycle 20
    done
done
