
``
# Polygenic risk score

```bash

# clump
# r 
for r in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
mkdir -p r${r}
done
# input file
echo SNP CHR BP A1 A2 P BETA > input

vi clump.sh 
#!/bin/sh
r=$1
plink \
--bfile Auto \
--clump input \
--clump-p1 1 \
--clump-p2 1 \
--clump-r2 ${r} \
--clump-kb 250 \
--out r${r}/250

chmod 755 clump.sh
module load hpcs

for r in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
sbatch -o r${r}/250.log ./clump.sh ${r}
done

```
# extract OR
```bash


vi extract_OR.R
-----------------
#Rscript
rsq <- commandArgs(trailingOnly=TRUE)[1]
gwassumstate <- commandArgs(trailingOnly=TRUE)[2]
data <- read.table(paste0("r",rsq,"/250.clumped"), header = T, as.is = T)
snp.list <- subset(data, select=c(SNP))
discovery <- read.table(gwassumstate, header=T)
discovery.clump <- merge(snp.list,discovery, by="SNP")
discovery.clump$Beta <- discovery.clump$BETA
write.table(discovery.clump[,c(1:6,8)],file=paste0("r",rsq,"/discovery.clumped"), sep="\t",row.names=F, quote=F )
#flip the effect allele of the summary results
data <- discovery.clump
head(data)
data$A1 <- as.character(data$A1)
data$A2 <- as.character(data$A2)
data$A1.flip <- ifelse(data$Beta >=0, data$A1,data$A2)
data$A2.flip <- ifelse(data$Beta >=0, data$A2,data$A1)
data$Beta.flip <- ifelse(data$Beta >=0, data$Beta,data$Beta* -1)
data.2 <- data[,c("SNP","CHR","BP","A1.flip","A2.flip","P","Beta","Beta.flip","A1","A2")]
names(data.2) <- c("SNP","CHR","BP","A1","A2","P","Beta.noflip","Beta","A1.noflip","A2.noflip")
write.table(data.2, file=paste0("r",rsq,"/discovery.flipBeta.clumped"), row.names=FALSE, col.names=TRUE,quote=FALSE, sep="\t")
#make the pvalue threshold file
S1 <- c(0,5*10^-8); S2 <- c(0,5*10^-7); S3 <- c(0,1*10^-6);S4 <- c(0,5*10^-6); S5 <- c(0,10^-5)
S6 <- c(0,5*10^-5); S7 <- c(0,5*10^-4); S8 <- c(0,5*10^-3); S9 <- c(0,0.01)
S10 <- c(0,0.05); S11 <- c(0,0.1); S12 <- c(0,0.2); S13 <- c(0,0.3)
S14 <- c(0,0.4); S15 <- c(0,0.5); S16 <- c(0,0.6); S17 <- c(0,0.7)
S18 <- c(0,0.8); S19 <- c(0,0.9);S20 <- c(0,1)
range <- as.data.frame(rbind(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,S19,S20))
write.table(range, file=paste0("r",rsq,"/range"), sep="\t", col.names= F,quote=F)
#check the number of SNPs within each threshold
range$snp.n <- NA
for(i in 1:nrow(range)){
        cat(range[i,2])
        cat("\n")
        range$snp.n[i] <- nrow(data[data$P < range[i,2],])
        print(nrow(data[data$P < range[i,2],]))
        cat("-------------------\n")
}
write.table(range, file=paste0("r",rsq,"/range.SNPnumber"), sep="\t",col.names= F,quote=F)
-----

chmod 755 extract_OR.R
module load hpcs

for r in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
sbatch -p s3 --mem=20G --qos snp-normal -o r${r}/extractOR.log ./extract_OR.R ${r} input
done


```
# Generate PRS
```bash

vi generate_PRS.sh
---------------------
#!/bin/sh
r=$1
plink \
--bfile Auto \
--keep sample \ 
--score r${r}/discovery.flipBeta.clumped 1 4 8 header sum double-dosage \
--q-score-range r${r}/range r${r}/discovery.flipBeta.clumped 1 6 header \
--out r${r}/prs
----------------------------------------

chmod 755 generate_PRS.sh
for r in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
sbatch -o r${r}/generatePRS.log ./generate_PRSv2.sh ${r}
done



