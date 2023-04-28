``
## GWAS using SAIGE

```bash


# NullModel (Step 1)
sbatch execSAIGE_step1.sh BED PHENO pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 GMMAT_Auto


sbatch execSAIGE_step1.sh
-------------
#!/bin/sh

BED=$1
PHENO=$2
COV=$3
OUTPUT=$4
/opt/local/r/3.5.1/R-3.5.1/bin/Rscript step1_fitNULLGLMM.R \
--plinkFile=${BED} \
--phenoFile=${PHENO} \
--phenoCol=affected \
--covarColList=${COV} \
--sampleIDColinphenoFile=IID \
--traitType=binary \
--outputPrefix=${OUTPUT} \
--nThreads=20 \
--LOCO=FALSE
---------------


--------------------
#!/bin/sh

CHR=$1
START=$2
END=$3
VCF=$4
PHENOFILE=$5
RDA=$6
VARIANCE=$7
LD_LIBRARY_PATH=/opt/local/htslib-1.9/:$LD_LIBRARY_PATH
module load gcc/6.1.0

zcat ${VCF} | head -15 | tail -1 | cut -f10- | tr \\t \\n | awk '{print "Y" $1}' > vcf_shurink/sampleIDindosageChr${CHR}.${START}_${END}.Y

/opt/local/r/3.5.1/R-3.5.1/bin/Rscript step2_SPAtests.R \
--vcfFile=${VCF} \
--vcfFileIndex=${VCF}.tbi \
--vcfField=DS \
--chrom=${CHR} \
--minMAF=0.005 \
--minMAC=1 \
--sampleFile=vcf_shurink/sampleIDindosageChr${CHR}.${START}_${END}.Y \
--GMMATmodelFile=${RDA} \
--varianceRatioFile=${VARIANCE} \
--SAIGEOutputFile=ccs/SPAtests_Chr${CHR}.${START}_${END}.dose.txt \
--numLinesOutput=10000 \
--IsOutputAFinCaseCtrl=TRUE


----------------------


```

### Meta-analysis using METAL 

```bash

/metal
GENOMICCONTROL OFF
SCHEME STDERR
EFFECTLABEL EFFECT
AVERAGEFREQ ON
FREQLABEL FREQ
MARKERLABEL MARKER
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N
PROCESS sakamoto.input
PROCESS suetsugu.input
VERBOSE ON
ANALYZE
