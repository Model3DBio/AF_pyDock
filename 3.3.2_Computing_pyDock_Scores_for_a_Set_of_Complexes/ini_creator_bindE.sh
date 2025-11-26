#!/bin/bash
# Case ID
export CASE="${CASE:-4POU}"
export PYDOCK="${PYDOCK:-/usr/local/software/pyDock3/}"
export GREASY="${GREASY:-/path/to/software/greasy/}"
export GREASY_NWORKERS=${GREASY_NWORKERS:-8}
export CHAINS_REC_LIG_VALUES=${CHAINS_REC_LIG_VALUES:-"A B"}

cd ${CASE}/AF2/ || exit

rm -f Greasy_BindE_mul_ligs.txt

for h in *cola*v*/ fold*/; do
    echo "$h"
    cd "$h" || exit

    for j in *[0-9].pdb; do
        echo "$j"

        for value in "${CHAINS_REC_LIG_VALUES[@]}"; do
            REC=$(echo "$value" | cut -d" " -f1)
            LIG=$(echo "$value" | cut -d" " -f2)
            CH=${LIG/,}

            ini_file="${j/.pdb}_LIG_${CH}.ini"

            cat <<EOF > "$ini_file"
[receptor]
pdb = $j
mol = $REC
newmol = $REC

[ligand]
pdb = $j
mol = $LIG
newmol = $LIG
EOF

            echo "cd ${h}; timeout 5m $PYDOCK/pyDock3 $ini_file bindEy; cd -" >> ../Greasy_BindE_mul_ligs.txt
        done
    done

    cd - || exit
done
${GREASY}/greasy Greasy_BindE_mul_ligs.txt
