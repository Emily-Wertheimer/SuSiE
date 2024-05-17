: '
Oct 18 2023
Chia-Yi Joy Lee (cl2375)
Version 2 LD-based clumping / LD calculation
'
### data cleaning for SumStats
## make sure all data have exact header as LD climping page state (mostly SNP and P for GWAS SumStats)
# change the heading from p_value to P, rsid to SNP
sed 's/rsid/SNP/;s/p_value/P/' GCST90027163_buildGRCh37.tsv > ALS_SumStats_title.tsv


### LD clumping
# load PLINK
module load PLINK/1.9b_6.21-x86_64
# remove HMC region
awk -F'\t' '!($3=="6" && $4 >= 28000000 && $4 <= 34000000)' AD_SumStats_title.tsv > AD_SumStats_chr6rm.txt #GRCh38
awk -F'\t' '!($2=="6" && $3 >= 28000000 && $3 <= 34000000)' ALS_SumStats_title.tsv > ALS_SumStats_chr6rm.txt #GRCh37
awk -F'\t' '!($1=="6" && $2 >= 28000000 && $2 <= 34000000)' PD_SumStats_sorted.tab > PD_SumStats_chr6rm.txt #GRCh37

# clumping -v1
plink --bfile /home/cl2375/project/GRCh37/plink1_file/clean_up/GRCh37.rmDups.snps \
--clump /home/cl2375/project/sumstats/AD/AD_SumStats_chr6rm.txt \
--clump-p1 0.00000005 --clump-p2 0.05 --clump-r2 0.2 --clump-kb 250

--> MS: 180 clumps formed from 3675 top variants (479 more top variant IDs missing: rsID are not part of file but chr:bp is)
--> AD: 362 clumps formed from 4732 top variants (429 top variant IDs missing)
        147 clumps when I use r2 > 0.01
--> ALS: 34 clumps formed from 489 top variants
# clumping -v2
plink --bfile /home/cl2375/project/GRCh37/plink1_file/original.EUR/1000G.EUR.plink \
--clump /home/cl2375/project/sumstats/AD/AD_SumStats_chr6rm.txt \
--clump-p1 0.00000005 --clump-p2 0.01 --clump-r2 0.2 --clump-kb 250
--> AD: 223 clumps
--> ALS: 23 clumps
--> PD: 42 clumps
--> MS: 129 clumps


### clean up clumping data for fine-mapping
# get SNPs and SP2
awk 'FNR >= 2 {print}' plink.clumped.ALS | awk '{print $3, $12}' > clumped_list_ALS_all.txt
# remove extra information
sed -i -e "s/ /\n/g" -e "s/,//g" -e "s/(1)/\n/g" clumped_list_ALS_all.txt
# clumping result: each locus are separated by either a blank line or NONE
# remove blank line at the end of the file
tac clumped_list_ALS_all.txt | awk 'NF {p=1} p' | tac > clumped_list_ALS_all_rmextra.txt
# calculate the amount of NONE in the file, meaning single SNP locus
#grep -w NONE clumped_list_ALS_all_rmextra.txt | wc -l 
# Remove all NONE in the file and replace with empty line
sed -i 's/NONE//g' clumped_list_ALS_all_rmextra.txt
# calculate the number of blank line
#grep -n ^$ clumped_list_ALS_all_rmextra.txt | wc -l
# separate file with empty line
awk 'BEGIN{file="ALS_locus_"++i".txt"} !NF{file="ALS_locus_"++i".txt";next} {print > file}' ../clumped_list_ALS_all_rmextra.txt

## LD calculation

#!/bin/bash
#SBATCH -J plink_LDmatrix
#SBATCH --mem-per-cpu 10G -c 4 -t 12:00:00
#SBATCH --mail-type ALL
#SBATCH --output LD_calculation_ALS_all.out

module load PLINK/1.9b_6.21-x86_64

for sample in /gpfs/gibbs/project/huckins/cl2375/finemap/LDbased_clumping/ALS_locus/*; do
    name=$(echo ${sample##*/} | sed s/ALS_locus_// | sed s/.txt// )
    plink --bfile /gpfs/gibbs/project/huckins/cl2375/GRCh37/plink1_file/clean_up/GRCh37.rmDups.snps \
    --extract $sample \
    --r --matrix \
    --out /gpfs/gibbs/project/huckins/cl2375/finemap/LD_calculation/LD_ALS/LDmatrix_ALS_$name
done
