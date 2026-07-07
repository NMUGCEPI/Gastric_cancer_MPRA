
###################################################################################################################
# MID
###################################################################################################################

#---------------------------------------Step 1: Step1_Cutadaptor
mkdir Step1_Cutadaptor
fastp \
-i MID_R1.fq.gz \
-I MID_R2.fq.gz \
-o ./Step1_Cutadaptor/TarBa_R1.cutAdaptor.fq.gz \
-O ./Step1_Cutadaptor/TarBa_R2.cutAdaptor.fq.gz \
-w 30 -q 10 -l 50 \
--correction --detect_adapter_for_pe \
-j ./Step1_Cutadaptor/TarBa.json -h ./Step1_Cutadaptor/TarBa.html

#---------------------------------------Step 2: Step2_MergeR1R2
mkdir Step2_MergeR1R2
flash \
./Step1_Cutadaptor/TarBa_R1.cutAdaptor.fq.gz \
./Step1_Cutadaptor/TarBa_R2.cutAdaptor.fq.gz \
-p 33 -r 150 -f 206 -s 8 -m 80 \
-o ./Step2_MergeR1R2/TarBa.merge \
-z -t 20

#---------------------------------------Step 3: Step3_Trim  
mkdir Step3_Trim
fastp \
-i ./Step2_MergeR1R2/TarBa.merge.extendedFrags.fastq.gz \
-o ./Step3_Trim/TarBa.trim.fq.gz \
-w 16 -q 30 -l 180 \
-j ./Step3_Trim/TarBa.trim.json \
-h ./Step3_Trim/TarBa.trim.html

#---------------------------------------Step 4: Step4_Rmdup
mkdir Step4_Rmdup
seqkit rmdup --by-seq --ignore-case ./Step3_Trim/TarBa.trim.fq.gz > ./Step4_Rmdup/TarBa.trim.rmdup.fq

#---------------------------------------Step 5: Step5_BWA
mkdir Step5_BWA

bwa index TAsnp_indelAA_noEXsite_ref.fa

split -l 20000000 ./Step4_Rmdup/TarBa.trim.rmdup.fq -d -a 3 --additional-suffix=.fastq ./Step5_BWA/rmdup_splitforbwa_

mkdir -p Step6_BWAout
mkdir -p logs

for i in {000..055}
do
    echo "Processing ${i} ..."

    bwa mem -t 20 TAsnp_indelAA_noEXsite_ref.fa \
        ./Step5_BWA/rmdup_splitforbwa_${i}.fastq \
        -L 100 -k 8 -O 5 |
    awk 'NR >= 12631 {print $3, $6, $10}' |
    awk '
    {
        match_sum = 0;
        cigar = $2
        n = split(cigar, arr, /[MSIDHPN=X]/, seps)
        for (j = 1; j <= n; ++j) {
            if (seps[j] == "M") {
                match_sum += arr[j]
            }
        }
        print $0, match_sum
    }' > ./Step6_BWAout/Oligo_barcode_${i}.out

if [[ -s ./Step6_BWAout/Oligo_barcode_${i}.out ]]; then
    rm -f ./Step5_BWA/rmdup_splitforbwa_${i}.fastq
fi

    echo "Finished ${i}"
done

#---------------------------------------Step 6: Step6_matchOligo_Barcode
mkdir Step6_matchOligo_Barcode

for i in {000..055}
do
awk '$1 == "*" { print $3 }' ./Step6_BWAout/Oligo_barcode_${i}.out > ./Step6_matchOligo_Barcode/NoIDmatched_${i}.seq
awk '$1 != "*" { print $1, $3, $4 }' ./Step6_BWAout/Oligo_barcode_${i}.out > ./Step6_matchOligo_Barcode/IDmatched_${i}.out
done

# match ID_1
R
library(data.table)
library(dplyr)

for (i in sprintf("%03d", 0:55)) {

IDmatch = data.frame(fread(paste0('./Step6_matchOligo_Barcode/IDmatched_',i,'.out'),header=F))
length = data.frame(fread('TAsnp_indelAA_noEXsite.length',header=F))
colnames(length)[2] = 'V4'

IDmatch_length = left_join(IDmatch, length, by = "V1")
IDmatch_length$V5 = IDmatch_length$V3 / IDmatch_length$V4

Match_less_095 = IDmatch_length[which(IDmatch_length$V5 < 0.95),]
fwrite(Match_less_095,paste0('./Step6_matchOligo_Barcode/Match_less_095_',i,'.out'),quo=F,row=F,col=F,sep=' ')

Match_more095_less100 = IDmatch_length[which(IDmatch_length$V5 < 1 & IDmatch_length$V5 >= 0.95),]
fwrite(Match_more095_less100,paste0('./Step6_matchOligo_Barcode/Match_more095_less100_',i,'.out'),quo=F,row=F,col=F,sep=' ')

Match_equa_100 = IDmatch_length[which(IDmatch_length$V5 == 1),]
Match_equa_100 = Match_equa_100[,c('V1','V2')]

ID_1 = data.frame(fread('ID_1.seq',header=F))
colnames(ID_1) = c('V1','V3','V4')
Match_equa_100_ID1 = left_join(Match_equa_100, ID_1, by = "V1")

fwrite(Match_equa_100_ID1,paste0('./Step6_matchOligo_Barcode/Match_equa_100_ID1_',i,'.out'),quo=F,row=F,col=F,sep=' ')

}
q()
n

for i in {000..055}
do

awk '{ gsub($4, "#", $2); print }' \
  ./Step6_matchOligo_Barcode/Match_equa_100_ID1_${i}.out \
  > ./Step6_matchOligo_Barcode/Match_equa_100_ID1_gsub_${i}.out

grep '#' ./Step6_matchOligo_Barcode/Match_equa_100_ID1_gsub_${i}.out \
  | sed 's/#/ /' \
  | awk '{print $4, $3, $2}' \
  > ./Step6_matchOligo_Barcode/Oligo_MatchedID1_${i}.out

grep -v '#' ./Step6_matchOligo_Barcode/Match_equa_100_ID1_gsub_${i}.out \
  | awk '{print $1, $2}' \
  > ./Step6_matchOligo_Barcode/Oligo_NoMatchedID1_${i}.out

done

# match ID_2
R
library(data.table)
library(dplyr)

for (i in sprintf("%03d", 0:55)) {

ID1nomatch = data.frame(fread(paste0('./Step6_matchOligo_Barcode/Oligo_NoMatchedID1_',i,'.out'),header=F))

ID_2 = data.frame(fread('ID_2.seq',header=F))
colnames(ID_2) = c('V1','V3','V4')
ID1nomatch_joinID2 = left_join(ID1nomatch, ID_2, by = "V1")

fwrite(ID1nomatch_joinID2,paste0('./Step6_matchOligo_Barcode/Match_equa_100_ID2_',i,'.out'),quo=F,row=F,col=F,sep=' ')

}
q()
n

for i in {000..055}
do

awk '{ gsub($4, "#", $2); print }' \
  ./Step6_matchOligo_Barcode/Match_equa_100_ID2_${i}.out \
  > ./Step6_matchOligo_Barcode/Match_equa_100_ID2_gsub_${i}.out

grep '#' ./Step6_matchOligo_Barcode/Match_equa_100_ID2_gsub_${i}.out \
  | sed 's/#/ /' \
  | awk '{print $4, $3, $2}' \
  > ./Step6_matchOligo_Barcode/Oligo_MatchedID2_${i}.out

grep -v '#' ./Step6_matchOligo_Barcode/Match_equa_100_ID2_gsub_${i}.out \
  | awk '{print $2}' \
  > ./Step6_matchOligo_Barcode/Oligo_NoMatchedID1ID2_${i}.out

done

# merge files
for i in {000..055}
do

awk '{print $2}' \
  ./Step6_matchOligo_Barcode/Match_less_095_${i}.out \
  > ./Step6_matchOligo_Barcode/Match_less_095_${i}.seq

rm ./Step6_matchOligo_Barcode/Match_less_095_${i}.out

awk '{print $2}' \
  ./Step6_matchOligo_Barcode/Match_more095_less100_${i}.out \
  > ./Step6_matchOligo_Barcode/Match_more095_less100_${i}.seq

rm ./Step6_matchOligo_Barcode/Match_more095_less100_${i}.out

done

cd ./Step6_matchOligo_Barcode
mkdir IDmatched
mv IDmatched_*.out ./IDmatched
mkdir Match_equa_100
mv Match_equa_100_ID1_*.out ./Match_equa_100  
mkdir Match_equa_100_ID2    
mv Match_equa_100_ID2_*.out ./Match_equa_100_ID2  
mkdir Match_less_095 
mv Match_less_095_*.seq ./Match_less_095  
mkdir Match_more095_less100 
mv Match_more095_less100_*.seq ./Match_more095_less100  
mkdir NoIDmatched 
mv NoIDmatched_*.seq ./NoIDmatched
mkdir Oligo_MatchedID1  
mv Oligo_MatchedID1_*.out ./Oligo_MatchedID1 
mkdir Oligo_NoMatchedID1
mv Oligo_NoMatchedID1_*.out ./Oligo_NoMatchedID1 
mkdir Oligo_MatchedID2
mv Oligo_MatchedID2_*.out ./Oligo_MatchedID2  
mkdir Oligo_NoMatchedID1ID2 
mv Oligo_NoMatchedID1ID2_*.out ./Oligo_NoMatchedID1ID2  
cd ..

#---------------------------------------Step 7: Step7_HindIII_XbaI_Barcode20_Seq
mkdir Step7_HindIII_XbaI_Barcode20_Seq

for i in {000..055}
do

awk '$2 ~ /^GCTT[ATGC]*TCTAGA[ATGC]{20}CC/' ./Step6_matchOligo_Barcode/Oligo_MatchedID1/Oligo_MatchedID1_${i}.out > ./Step7_HindIII_XbaI_Barcode20_Seq/Oligo_MatchedID1_DEX_Barcode20_${i}.out
awk '$2 ~ /^GCTT[ATGC]*TCTAGA[ATGC]{20}CC/' ./Step6_matchOligo_Barcode/Oligo_MatchedID2/Oligo_MatchedID2_${i}.out > ./Step7_HindIII_XbaI_Barcode20_Seq/Oligo_MatchedID2_DEX_Barcode20_${i}.out

grep -E 'AAGCTT[ATCG]*TCTAGA[(?=A|T|G|C)]{20}CC' ./Step6_matchOligo_Barcode/Oligo_NoMatchedID1ID2/Oligo_NoMatchedID1ID2_${i}.out > ./Step7_HindIII_XbaI_Barcode20_Seq/Oligo_NoMatchedID1ID2_DEX_Barcode20_${i}.out
grep -E 'AAGCTT[ATCG]*TCTAGA[(?=A|T|G|C)]{20}CC' ./Step6_matchOligo_Barcode/Match_more095_less100/Match_more095_less100_${i}.seq > ./Step7_HindIII_XbaI_Barcode20_Seq/Match_more095_less100_DEX_Barcode20_${i}.seq
grep -E 'AAGCTT[ATCG]*TCTAGA[(?=A|T|G|C)]{20}CC' ./Step6_matchOligo_Barcode/NoIDmatched/NoIDmatched_${i}.seq > ./Step7_HindIII_XbaI_Barcode20_Seq/NoIDmatched_DEX_Barcode20_${i}.seq
grep -E 'AAGCTT[ATCG]*TCTAGA[(?=A|T|G|C)]{20}CC' ./Step6_matchOligo_Barcode/Match_less_095/Match_less_095_${i}.seq > ./Step7_HindIII_XbaI_Barcode20_Seq/Match_less_095_DEX_Barcode20_${i}.seq
done

#---------------------------------------Step 8: Step8_LV_Distance
mkdir Step8_LV_Distance

for i in {000..055}
do
awk '{
    match($2, /^GCTT[ATGC]*TCTAGA[ATGC]{20}CC/, arr)
    if (match($2, /^GCTT[ATGC]*TCTAGA[ATGC]{20}CC/)) {
        print $0, arr[0]
    }
}' ./Step7_HindIII_XbaI_Barcode20_Seq/Oligo_MatchedID1_DEX_Barcode20_${i}.out > ./Step8_LV_Distance/Oligo_MatchedID1_DEX_Barcode20_${i}.tmp1
done

for i in {000..055}
do
awk '{
    match($2, /^GCTT[ATGC]*TCTAGA[ATGC]{20}CC/, arr)
    if (match($2, /^GCTT[ATGC]*TCTAGA[ATGC]{20}CC/)) {
        print $0, arr[0]
    }
}' ./Step7_HindIII_XbaI_Barcode20_Seq/Oligo_MatchedID2_DEX_Barcode20_${i}.out > ./Step8_LV_Distance/Oligo_MatchedID2_DEX_Barcode20_${i}.tmp1
done

for i in {000..055}
do
# down_arm
awk '{ gsub($4, "", $2); print }' ./Step8_LV_Distance/Oligo_MatchedID1_DEX_Barcode20_${i}.tmp1 > ./Step8_LV_Distance/Oligo_MatchedID1_DEX_Barcode20_${i}.tmp2
awk '{ gsub($4, "", $2); print }' ./Step8_LV_Distance/Oligo_MatchedID2_DEX_Barcode20_${i}.tmp1 > ./Step8_LV_Distance/Oligo_MatchedID2_DEX_Barcode20_${i}.tmp2

# Barcode
awk '{print $0, substr($4, length($4)-21, 20)}' ./Step8_LV_Distance/Oligo_MatchedID1_DEX_Barcode20_${i}.tmp2 > ./Step8_LV_Distance/Oligo_MatchedID1_DEX_Barcode20_${i}.tmp3
awk '{print $0, substr($4, length($4)-21, 20)}' ./Step8_LV_Distance/Oligo_MatchedID2_DEX_Barcode20_${i}.tmp2 > ./Step8_LV_Distance/Oligo_MatchedID2_DEX_Barcode20_${i}.tmp3

# XbaI/HindIII
awk '{print $1, $3, $2, substr($4, 1, length($4)-22), $5}' ./Step8_LV_Distance/Oligo_MatchedID1_DEX_Barcode20_${i}.tmp3 > ./Step8_LV_Distance/Oligo_MatchedID1_DEX_Barcode20_${i}.tmp4
awk '{print $1, $3, $2, substr($4, 1, length($4)-22), $5}' ./Step8_LV_Distance/Oligo_MatchedID2_DEX_Barcode20_${i}.tmp3 > ./Step8_LV_Distance/Oligo_MatchedID2_DEX_Barcode20_${i}.tmp4

awk '{print $1, $2 "TA", "CC" $3, "AA" $4, $5}' ./Step8_LV_Distance/Oligo_MatchedID1_DEX_Barcode20_${i}.tmp4 > ./Step8_LV_Distance/Oligo_MatchedID1_DEX_Barcode20_${i}.forLVDcaculation
awk '{print $1, $2 "TA", "CC" $3, "AA" $4, $5}' ./Step8_LV_Distance/Oligo_MatchedID2_DEX_Barcode20_${i}.tmp4 > ./Step8_LV_Distance/Oligo_MatchedID2_DEX_Barcode20_${i}.forLVDcaculation
done

cd ./Step8_LV_Distance
rm *.tmp1
rm *.tmp2
rm *.tmp3
rm *.tmp4
cd ..

# LV_Distance
R
library(stringdist)
library(data.table)

for (i in sprintf("%03d", 0:55)) {

ID1 = data.frame(fread(paste0('./Step8_LV_Distance/Oligo_MatchedID1_DEX_Barcode20_',i,'.forLVDcaculation'),header=F))
ID2 = data.frame(fread(paste0('./Step8_LV_Distance/Oligo_MatchedID2_DEX_Barcode20_',i,'.forLVDcaculation'),header=F))

UpArm= "CGGTACCTGAGCTCGCTA"
DownArm = "CCGCTTCGAGCAGACATG" 
DEX = "AAGCTTGGGCTAGTCTAGA"

# calculate Levenshtein Distance
ID1$UpArmLV = stringdist(ID1$V2, UpArm, method = 'lv')
ID1$DownArmLV = stringdist(ID1$V3, DownArm, method = 'lv')
ID1$DEXLV = stringdist(ID1$V4, DEX, method = 'lv')
ID2$UpArmLV = stringdist(ID2$V2, UpArm, method = 'lv')
ID2$DownArmLV = stringdist(ID2$V3, DownArm, method = 'lv')
ID2$DEXLV = stringdist(ID2$V4, DEX, method = 'lv')

# Levenshtein Distance 3 or less
ID1_LVmatch = ID1[which(ID1$UpArmLV <=3 & ID1$DownArmLV <=3 & ID1$DEXLV <=3),]
ID1_LVnomatch = ID1[which(ID1$UpArmLV > 3 | ID1$DownArmLV > 3 | ID1$DEXLV > 3),]
ID2_LVmatch = ID2[which(ID2$UpArmLV <=3 & ID2$DownArmLV <=3 & ID2$DEXLV <=3),]
ID2_LVnomatch = ID2[which(ID2$UpArmLV > 3 | ID2$DownArmLV > 3 | ID2$DEXLV > 3),]
ID1_LVmatch = ID1_LVmatch[,c('V1','V5')]
ID2_LVmatch = ID2_LVmatch[,c('V1','V5')]
ID1_LVnomatch = ID1_LVnomatch[,c('V1','V5')]
ID2_LVnomatch = ID2_LVnomatch[,c('V1','V5')]

fwrite(ID1_LVmatch,paste0('./Step8_LV_Distance/Oligo_ID1_Barcode_',i,'.passLV'),quo=F,row=F,col=F,sep=' ')
fwrite(ID2_LVmatch,paste0('./Step8_LV_Distance/Oligo_ID2_Barcode_',i,'.passLV'),quo=F,row=F,col=F,sep=' ')
fwrite(ID1_LVnomatch,paste0('./Step8_LV_Distance/Oligo_ID1_Barcode_',i,'.failLV'),quo=F,row=F,col=F,sep=' ')
fwrite(ID2_LVnomatch,paste0('./Step8_LV_Distance/Oligo_ID2_Barcode_',i,'.failLV'),quo=F,row=F,col=F,sep=' ')

}
q()
n

cat ./Step8_LV_Distance/*.passLV > ./Step8_LV_Distance/OligoID_Barcode.passed
cat ./Step8_LV_Distance/*.failLV | awk '{print $2}' > ./Step8_LV_Distance/Barcode.failed.UpArm.DownArm.DEX.lv
rm ./Step8_LV_Distance/Oligo_ID*_Barcode_*.passLV
rm ./Step8_LV_Distance/Oligo_ID*_Barcode_*.failLV

sort ./Step8_LV_Distance/Barcode.failed.UpArm.DownArm.DEX.lv | uniq > ./Step8_LV_Distance/Barcode.failed.UpArm.DownArm.DEX.sort.uniq.lv
sort ./Step8_LV_Distance/OligoID_Barcode.passed | uniq > ./Step8_LV_Distance/OligoID_Barcode.sort.uniq.passed

#--------------------------------------- Step 9: Step9_BarcodeNeedExclude
mkdir Step9_BarcodeNeedExclude

# < 95%
cat ./Step6_matchOligo_Barcode/Match_less_095/Match_less_095_*.seq  > ./Step9_BarcodeNeedExclude/Match_less_095.seq.tmp
cat ./Step6_matchOligo_Barcode/NoIDmatched/NoIDmatched_*.seq  > ./Step9_BarcodeNeedExclude/NoIDmatched.seq.tmp
cat ./Step9_BarcodeNeedExclude/Match_less_095.seq.tmp ./Step9_BarcodeNeedExclude/NoIDmatched.seq.tmp > ./Step9_BarcodeNeedExclude/Match_less_095.seq
rm ./Step9_BarcodeNeedExclude/Match_less_095.seq.tmp
rm ./Step9_BarcodeNeedExclude/NoIDmatched.seq.tmp
find . -name "Match_less_095.seq" | xargs grep -P 'AAGCTT[(?=A|T|G|C)]*TCTAGA[(?=A|T|G|C)]{20}CC' -o  > ./Step9_BarcodeNeedExclude/Barcodes.Match_less_095.tmp
awk '{print substr($0, length($0)-21, 20)}' ./Step9_BarcodeNeedExclude/Barcodes.Match_less_095.tmp > ./Step9_BarcodeNeedExclude/Barcodes.Match_less_095

# 95% - 100%
cat ./Step6_matchOligo_Barcode/Match_more095_less100/Match_more095_less100_*.seq > ./Step9_BarcodeNeedExclude/Match_more095_less100.seq
find . -name "Match_more095_less100.seq" | xargs grep -P 'AAGCTT[(?=A|T|G|C)]*TCTAGA[(?=A|T|G|C)]{20}CC' -o  > ./Step9_BarcodeNeedExclude/Barcodes.Match_more095_less100.tmp
awk '{print substr($0, length($0)-21, 20)}' ./Step9_BarcodeNeedExclude/Barcodes.Match_more095_less100.tmp > ./Step9_BarcodeNeedExclude/Barcodes.Match_more095_less100

# = 100%
cat ./Step6_matchOligo_Barcode/Oligo_NoMatchedID1ID2/Oligo_NoMatchedID1ID2_*.out > ./Step9_BarcodeNeedExclude/Oligo_NoMatchedID1ID2.seq
find . -name "Oligo_NoMatchedID1ID2.seq" | xargs grep -P 'AAGCTT[(?=A|T|G|C)]*TCTAGA[(?=A|T|G|C)]{20}CC' -o  > ./Step9_BarcodeNeedExclude/Barcodes.Oligo_NoMatchedID1ID2.tmp
awk '{print substr($0, length($0)-21, 20)}' ./Step9_BarcodeNeedExclude/Barcodes.Oligo_NoMatchedID1ID2.tmp > ./Step9_BarcodeNeedExclude/Barcodes.Oligo_NoMatchedID1ID2

rm ./Step9_BarcodeNeedExclude/Match_more095_less100.seq
rm ./Step9_BarcodeNeedExclude/Match_less_095.seq
rm ./Step9_BarcodeNeedExclude/Oligo_NoMatchedID1ID2.seq
rm ./Step9_BarcodeNeedExclude/Barcodes.Match_more095_less100.tmp
rm ./Step9_BarcodeNeedExclude/Barcodes.Match_less_095.tmp
rm ./Step9_BarcodeNeedExclude/Barcodes.Oligo_NoMatchedID1ID2.tmp

mv ./Step8_LV_Distance/Barcode.failed.UpArm.DownArm.DEX.sort.uniq.lv ./Step9_BarcodeNeedExclude/Barcode.failedLV

cat \
  ./Step9_BarcodeNeedExclude/Barcode.failedLV \
  ./Step9_BarcodeNeedExclude/Barcodes.Match_more095_less100 \
  ./Step9_BarcodeNeedExclude/Barcodes.Oligo_NoMatchedID1ID2 \
  > ./Step9_BarcodeNeedExclude/Barcode.loose.needExclude

sort -T ./ ./Step9_BarcodeNeedExclude/Barcode.loose.needExclude | uniq > ./Step9_BarcodeNeedExclude/Barcode.loose.sort.uniq.needExclude

rm ./Step9_BarcodeNeedExclude/Barcode.failedLV
rm ./Step9_BarcodeNeedExclude/Barcodes.Match_more095_less100
rm ./Step9_BarcodeNeedExclude/Barcodes.Oligo_NoMatchedID1ID2

sort ./Step9_BarcodeNeedExclude/Barcodes.Match_less_095 | uniq > ./Step9_BarcodeNeedExclude/Barcode.Matchless095.sort.uniq.needExclude

#--------------------------------------- Step 10: Step10_Barcode_with_Exclude
mkdir Step10_Barcode_with_Exclude

awk 'FNR==NR{a[$1]=1; next} !($2 in a)' \
  ./Step9_BarcodeNeedExclude/Barcode.Matchless095.sort.uniq.needExclude \
  ./Step10_Barcode_with_Exclude/Oligo_Barcode.loose.exclude \
  > ./Step10_Barcode_with_Exclude/Oligo_Barcode.extra.exclude

#--------------------------------------- Step 11:  Step11_Barcode_Specific
mkdir Step11_Barcode_Specific

sort ./Step10_Barcode_with_Exclude/Oligo_Barcode.extra.exclude | uniq > ./Step11_Barcode_Specific/Oligo_Barcode.extra.sort.uniq.exclude

awk '{
    count[$2]++;
    lines[NR] = $0;
    keys[NR] = $2;
}
END {
    for (i = 1; i <= NR; i++)
        if (count[keys[i]] == 1)
            print lines[i];
}' ./Step11_Barcode_Specific/Oligo_Barcode.extra.sort.uniq.exclude > ./Step11_Barcode_Specific/Oligo_SpecificBarcode.strict.tsv

awk '{print $1}' \
  ./Step11_Barcode_Specific/Oligo_SpecificBarcode.strict.tsv \
  | uniq -c \
  > ./Step11_Barcode_Specific/Oligo.SpeciBarcCount.strict

awk '{print $2, $1}' \
  ./Step11_Barcode_Specific/Oligo_SpecificBarcode.strict.tsv \
  | sort \
  > ./Step11_Barcode_Specific/SpecificBarcode.sort_Oligo.strict.tsv

