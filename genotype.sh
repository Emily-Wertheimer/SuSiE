## convert plink2 file into plink1
plink2 --pgen all_phase3.pgen --pvar all_phase3.pvar.zst --psam phase3_corrected.psam  --make-bed --out GRCh37

### plink 1.9 load module 
module load PLINK/1.90-beta5.3

## remove duplicates in genotype bfile - batch job
cut -f 2 GRCh37.bim | sort | uniq -d > GRCh37.bim.dups
plink --bfile GRCh37  --exclude GRCh37.bim.dups --make-bed --out GRCh37.rmDups --allow-extra-chr

## remove indel and structural variants - batch job
plink --bfile GRCh37.rmDups --snps-only --make-bed --out GRCh37.rmDups.snps --allow-extra-chr

## find variant ID for all SNPs in summary statistic data -batch job
# make files executable
chmod u+x variant_ID_conversion_v3.2.sh
# input: location of summary statistic, location of genotype file, which chromosome, bp column in SS, variant name column in SS, bp column in bim, variant name column in bim
sbatch variant_ID_conversion_v3.2.sh /home/cl2375/scratch60/summary_statistics/MS_SS/MS_SS_1.csv /home/cl2375/scratch60/GRCh37/plink1_file/clean_up/GRCh37.rmDups.snps.bim 1 2 3 4 2
--> chr22 take 9h30, chr1/2 take 10D15h


### use EUR only 1000G data
## merging plink file
'''
#!/bin/bash
#SBATCH -J merge_plink
#SBATCH --mem-per-cpu 10G -c 4 -t 12:00:00
#SBATCH --mail-type ALL

rm mergelist.txt

for i in {1..22}
do
    echo 1000G.EUR.QC.$i >> mergelist.txt
done

module load PLINK/1.9b_6.21-x86_64
plink --merge-list mergelist.txt --make-bed --out 1000G.EUR.plink
'''
## check if SNPs are mostly included in the file
# v1 - only check the SNP and get the ratio
# v2 - try running in parallel
'''
#!/bin/bash
#SBATCH -J check_SNP
#SBATCH --mem-per-cpu 10G -c 4 -t 12:00:00
#SBATCH --mail-type ALL

# Path to your files
MS_LOCUS_DIR="/home/cl2375/project/finemap/LDbased_clumping/MS_locus"
PLINK_BIM_FILE="/home/cl2375/project/GRCh37/plink1_file/original.EUR/1000G.EUR.plink.bim"

# Initialize counters
total_snps=0
snps_in_reference=0

# Check each MS_locus_* file
for file in ${MS_LOCUS_DIR}/MS_locus_*; do
    filename=$(basename "$file")

    # Check SNPs in MS_locus file against 1000G.EUR.plink.bim
    while IFS= read -r line; do
        snp=$(echo "$line" | cut -f 2)
        ((total_snps++))
        if grep -q -w "$snp" "$PLINK_BIM_FILE"; then
            ((snps_in_reference++))
        fi
    done < "$file"
done

snps_not_in_reference=$((total_snps - snps_in_reference))
echo "Total SNPs: $total_snps"
echo "SNPs in reference file ($PLINK_BIM_FILE): $snps_in_reference"
echo "SNPs not in reference file: $snps_not_in_reference"
echo "Ratio of SNPs in reference file: $((snps_in_reference * 100 / total_snps))%"
echo "Ratio of SNPs not in reference file: $((snps_not_in_reference * 100 / total_snps))%"
'''
'''
Total SNPs: 8943
SNPs in reference file (/home/cl2375/project/GRCh37/plink1_file/original.EUR/1000G.EUR.plink.bim): 8863
SNPs not in reference file: 80
Ratio of SNPs in reference file: 99%
Ratio of SNPs not in reference file: 0%
'''
### was going to clean up with the same step reviously used, but turn out this is already cleaned






------ temp ------


### v2 is so far not working -- need to talk to HPC about it
#!/bin/bash
#SBATCH -J check_SNP
#SBATCH --mem-per-cpu 10G -c 4 -t 12:00:00
#SBATCH --mail-type ALL

# Path to your files
MS_LOCUS_DIR="/home/cl2375/project/finemap/LDbased_clumping/MS_locus"
PLINK_BIM_FILE="/home/cl2375/project/GRCh37/plink1_file/original.EUR/1000G.EUR.plink.bim"

# Function to check SNPs in one file
check_snps() {
    file=$1
    local total_snps=0
    local snps_in_reference=0

    echo "Checking SNPs in $file"

    while IFS= read -r snp; do
        if [ -n "$snp" ]; then
            echo "SNP found: $snp"
            ((total_snps++))
            if grep -q -w "$snp" "$PLINK_BIM_FILE"; then
                ((snps_in_reference++))
            fi
        else
            echo "No SNP found in the line"
        fi
    done < "$file"

    snps_not_in_reference=$((total_snps - snps_in_reference))
    echo "Total SNPs in $file: $total_snps"
    echo "SNPs in reference file ($PLINK_BIM_FILE): $snps_in_reference"
    echo "SNPs not in reference file: $snps_not_in_reference"
    echo "Ratio of SNPs in reference file: $((snps_in_reference * 100 / total_snps))%"
    echo "Ratio of SNPs not in reference file: $((snps_not_in_reference * 100 / total_snps))%"
}

export -f check_snps

module load parallel/20210322-GCCcore-10.2.0
ls ${MS_LOCUS_DIR}/MS_locus_* | parallel -j0 check_snps {}





## 2023.10.28 plan
# merge those files in /home/cl2375/palmer_scratch/LSDC/ldsc_webfiles/1000G_EUR_Phase3_plink for bed bim fam
# the LD score you might want: 
# from LDSC: /gpfs/gibbs/pi/huckins/software/eur_w_ld_chr/ #106137
# from private: /gpfs/gibbs/pi/huckins/software/LDSC_ref_PGCED/hc1kgp3.b38.eur.l2.jz2023.chr/ #2050831