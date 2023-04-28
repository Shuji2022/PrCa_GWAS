``
### Meta-analysis using MANTRA

```bash


## Making Datasets of SNPs
{SNP CHR BP A1(ALT) A2(REF) A1_Freq A1_BETA SE P Case Ctrl}

MANTRA_BBJ
MANTRA_EUR
MANTRA_AFR
MANTRA_HIS

cat MANTRA_BBJ | sed 1,1d | sort -k1,1 | join - <(cat MANTRA_EUR | sed 1,1d | awk '{print $1,$6,$7,$8,$9,$10}' | sort -k1,1) | join - <(cat MANTRA_AFR | sed 1,1d | awk '{print $1,$6,$7,$8,$9,$10}' | sort -k1,1) | join - <(cat MANTRA_HIS | sed 1,1d | awk '{print $1,$6,$7,$8,$9,$10}' | sort -k1,1) > MANTRA_combine



# map file
cat MANTRA_combine | awk '{print $1,$2,$3,$4,$5}' | sed -e 's/ /\t/g' > map
# dat file
## BBJ
cat MANTRA_combine | awk '{print 1,$6,$7,$8,$9,$10}' | sed -e 's/ /\t/g' > MANTRA.J.dat
## EUR 
cat MANTRA_combine | awk '{print 1,$11,$12,$13,$14,$15}' | sed -e 's/ /\t/g' > MANTRA.E.dat
## AFR
cat MANTRA_combine | awk '{print 1,$16,$17,$18,$19,$20}' | sed -e 's/ /\t/g' > MANTRA.A.dat
## HIS
cat MANTRA_combine | awk '{print 1,$21,$22,$23,$24,$25}' | sed -e 's/ /\t/g' > MANTRA.H.dat


#mantra.in file
vi mantra.in
MANTRA.J.dat    J
MANTRA.E.dat    E
MANTRA.A.dat    A
MANTRA.H.dat    H

```

# MANTRA

mantra.v2
mantra.v2.f95
dmatcal.v2
dmatcal.v2.f95


```bash


vi exe.sh
-----------------------------------------
#!/bin/bash

./mantra.v2 < fname"$1".in

---------------------------------------

chmod 755 mantra.v2
chmod 755 dmatcal.v2
chmod 755 exe.sh

# creat file.in
numchunk=`wc -l map | awk '{print int($1/10000)}'`
for i in `seq 0 "$numchunk"`
do
num=`expr $i \* 10000`
echo -e "map\n"$num"\n10000\nout_"$i".txt\nout_"$i".bet.txt\n0" > fname"$i".in
done

./dmatcal.v2
map

module load hpcs
#test
for i in 0; do sbatch -c 1 --qos snp-normal -p s1 exe.sh ${i} ;usleep 10000; done
# run mantra
for i in $(seq 1 -); do sbatch -c 1 --qos snp-normal -p s1  exe.sh ${i} ;usleep 10000; done
