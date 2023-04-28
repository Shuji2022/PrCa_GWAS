``
# Sample QC

```bash


HLA_START=24000000
HLA_STOP=36000000

# HLA SNP 
gawk '$1==6 && $4 > '${HLA_START}' && $4 < '${HLA_STOP}'{print $2}' Auto.bim | sort > qc/hla.snps

# Hardy Weinberg SNP
plink --bfile Auto --autosome --missing --out Auto
plink --bfile Auto --hardy --out Auto

# call rate > 99 %
gawk 'NR>1 && $5 > 0.01{print $2}' Auto.lmiss > Auto.miss01.snps
# HWE > 1.0 * 10^-6
gawk '/ALL/ && $9<1e-6{print $2}' Auto.hwe > Auto.hwp1e06.snps
done

# SNPs for sample QC

cat Auto.miss01.snps Auto.hwp1e06.snps hla.snps| sort | uniq > Auto.miss01+hwp1e06+hla.snps



```

# Identical sample check

```bash

plink --bfile Auto --autosome --exclude Auto.miss01+hwp1e06+hla.snps --keep QC.keep.sexok_1st_muscle --maf 0.05 --indep 50 5 2 --out Auto

plink --bfile Auto --extract Auto.prune.in --genome gz --min 0.1 --out Auto

plink --bfile Auto --exclude Auto.miss01+hwp1e06+hla.snps --maf 0.01 --missing --out Auto.maf01


IBD_THRES=0.9

python
    --genome Auto.genome.gz \
    --imiss Auto.maf01.imiss \
    --thres ${IBD_THRES} \
    --out Auto.relates_PIHAT.out


sort -k1,1 Auto.relates_PIHAT.out | join -v1 <(cat Auto.fam | awk '{print $1}' | sort) - | awk '{print $1}' > QC_ExcludeClose.sample


```
# Principal compornent analysis

```bash

Excluding outliers based on Hapmap PCA data


```
## SNP level QC

```bash

plink --bfile Auto --keep QC.keep --autosome --missing --out auto.sampQC
plink --bfile Auto --keep QC.keep.pcaok --autosome --hardy  --out auto.sampQC
plink --bfile Auto --keep QC.keep.pcaok --autosome --freq --out auto.sampQC


# call rate > 99 %
gawk '$5 > 0.01{print $2}' auto.sampQC.lmiss > auto.miss01.sampQC.snps
# HWE > 1.0 * 10^-6
gawk '$9<1e-6{print $2}' auto.sampQC.hwe > auto.hwp1e06.sampQC.snps
gawk '$5 < 0.01{print $2}' auto.sampQC.frq > auto.maf0.01.sampQC.snps
done


```

## Call Rate


```bash
# Call rate < 98%
plink --bfile auto --missing --out auto2

gawk 'NR>1 && $6 < 0.02{print $2}' auto2.imiss | awk '{print $1}' | sort -k1,1 > CallRateError.sample


```
# Imputation


```bash

for CHR in {1..22}
do
plink --bfile auto2 \
	--chr $CHR \
	--make-bed \
	--out auto.$CHR
done


# Phasing

for CHR in {1..22}
do
eagle  \
--bfile=${BED_GWAS_NAME} \
--geneticMapFile=${HOMEDIR}/Applications/Eagle/Eagle/tables/genetic_map_hg19_withX.txt.gz \
--outPrefix=Eagle_${FILE_GWAS_NAME} \
--maxMissingPerSnp=1 \
--maxMissingPerIndiv=1 \
--numThreads=${THR} \
--chrom=${CHR}
done

# Convert to vcf
for CHR in {1..22}
do
zcat ${DIR_EAGLE}/Eagle_${FILE_GWAS_NAME}.haps.gz > Eagle_${FILE_GWAS_NAME}.haps
cp ${DIR_EAGLE}/Eagle_${FILE_GWAS_NAME}.sample haps
shapeit \
	-convert \
	--input-haps Eagle_${FILE_GWAS_NAME} \
	--output-vcf Eagle_${FILE_GWAS_NAME}.vcf \
	--thread 10
done


for CHR in {1..22}
do
gzip Eagle_${FILE_GWAS_NAME}_${CHR}.vcf
done


## Imputation
for CHR in {1..22}
do
minimac4 \
	--refHaps \
    --haps ${VCF} \
	--prefix ${OUT} \
	--format GT,DS,GP \
	--cpus 10
done

