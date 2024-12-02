#!/bin/bash
#SBATCH --partition=local.q
#SBATCH --nodelist=galaxy3
#SBATCH -c 16
#SBATCH -o dynamic_contacts.log
#SBATCH -e dynamic_contacts.log


entries=("2ycw_A" "2ycw_B" "3qak_A" "4eiy_A" "5wf5_A" "6gdg" "6gdg_A")

for entry in "${entries[@]}"; do
    for i in {0..2}; do
        python ./get_dynamic_contacts.py --cores 8  --topology ../../data/water_md/md_nolig/${entry}_${i}/amber/step5_input.parm7 \
                                --trajectory ../../data/water_md/md_nolig/${entry}_${i}/amber/production_2fs_gpu/prod_50ns.nc \
                                --itypes all \
                                --ligand "resname 'Na+'" \
                                --output ../../data/water_md/md_nolig_interaction/${entry}_${i}_contacts.tsv
    done
done
