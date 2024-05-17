#!/bin/bash
#SBATCH --array 0-3
#SBATCH --job-name clump_LDmatrix
#SBATCH --mem-per-cpu 10G -c 4 -t 1-00:00:00
#SBATCH --mail-type ALL
#SBATCH --output LD_clump_matrix.out

SS_DIR="/home/cl2375/project/sumstats"
CLUMP_DIR="/home/cl2375/project/finemap/LD_clumping"
MATRIX_DIR="/home/cl2375/project/finemap/LD_calculation"
sample=("MS" "ALS" "AD" "PD")

for i in $SLURM_ARRAY_TASK_ID; do

    ### LD clumping
    module load PLINK/1.9b_6.21-x86_64
    plink --bfile /home/cl2375/project/GRCh37/plink1_file/original.EUR/1000G.EUR.plink \
    --clump ${SS_DIR}/${sample[${i}]}/*_chr6rm.* \
    --clump-p1 0.00000005 --clump-p2 0.01 --clump-r2 0.2 --clump-kb 250 \
    --out ${sample[${i}]}

    ### clumping info extraction
    awk 'FNR >= 2 {print}' ${sample[${i}]}.clumped | awk '{print $3, $12}' > clumped_list_${sample[${i}]}_all.txt
    # remove extra information
    sed -i -e "s/ /\n/g" -e "s/,//g" -e "s/(1)/\n/g" clumped_list_${sample[${i}]}_all.txt
    # clumping result: each locus are separated by either a blank line or NONE
    # remove blank line at the end of the file
    tac clumped_list_${sample[${i}]}_all.txt | awk 'NF {p=1} p' | tac > clumped_list_${sample[${i}]}_all_rmextra.txt
    # calculate the amount of NONE in the file, meaning single SNP locus
    #grep -w NONE clumped_list_ALS_all_rmextra.txt | wc -l 
    # Remove all NONE in the file and replace with empty line
    sed -i 's/NONE//g' clumped_list_${sample[${i}]}_all_rmextra.txt
    # calculate the number of blank line
    #grep -n ^$ clumped_list_ALS_all_rmextra.txt | wc -l
    # move into clumping directory
    cd ${CLUMP_DIR}/${sample[${i}]}_clump
    # separate file with empty line
    awk -v sample_var="${sample[${i}]}" 'BEGIN{file=sample_var "_locus_"++i".txt"} !NF{file=sample_var "_locus_"++i".txt";next} {print > file}' ../../clumped_list_${sample[${i}]}_all_rmextra.txt

    ### LD calculation
    cd ${MATRIX_DIR}
    for y in ${CLUMP_DIR}/${sample[${i}]}_clump/*; do
        name=$(echo ${y##*/} | sed s/${sample[${i}]}_locus_// | sed s/.txt// )
        plink --bfile /home/cl2375/project/GRCh37/plink1_file/original.EUR/1000G.EUR.plink \
        --extract $y \
        --r --matrix \
        --out ${MATRIX_DIR}/${sample[${i}]}_LD/LDmatrix_${sample[${i}]}_$name
    done
done



