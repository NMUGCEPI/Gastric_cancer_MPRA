
#--------------------------------------------------Step1_Cutadaptor
mkdir Step1_Cutadaptor  
for i in {1..16}
do
fastp -i ../Rawdata/cDNA${i}/cDNA${i}_R1.fq.gz -I ../Rawdata/cDNA${i}/cDNA${i}_R2.fq.gz -o ./Step1_Cutadaptor/RNA${i}_R1.cutAdaptor.fq.gz -O ./Step1_Cutadaptor/RNA${i}_R2.cutAdaptor.fq.gz -w 5 -q 10 -l 50 --correction --detect_adapter_for_pe -j ./Step1_Cutadaptor/RNA${i}.json -h ./Step1_Cutadaptor/RNA${i}.html
fastp -i ../Rawdata/DNA${i}/DNA${i}_R1.fq.gz -I ../Rawdata/DNA${i}/DNA${i}_R2.fq.gz -o ./Step1_Cutadaptor/DNA${i}_R1.cutAdaptor.fq.gz -O ./Step1_Cutadaptor/DNA${i}_R2.cutAdaptor.fq.gz -w 5 -q 10 -l 50 --correction --detect_adapter_for_pe -j ./Step1_Cutadaptor/DNA${i}.json -h ./Step1_Cutadaptor/DNA${i}.html
done

#--------------------------------------------------Step2_MergeR1R2
mkdir Step2_MergeR1R2   
for i in {1..16}
do
flash ./Step1_Cutadaptor/RNA${i}_R1.cutAdaptor.fq.gz ./Step1_Cutadaptor/RNA${i}_R2.cutAdaptor.fq.gz -p 33 -r 150 -f 129 -s 5 -m 110 -o ./Step2_MergeR1R2/RNA${i}.fq.gz -z -t 5
flash ./Step1_Cutadaptor/DNA${i}_R1.cutAdaptor.fq.gz ./Step1_Cutadaptor/DNA${i}_R2.cutAdaptor.fq.gz -p 33 -r 150 -f 129 -s 5 -m 110 -o ./Step2_MergeR1R2/DNA${i}.fq.gz -z -t 5
done

#--------------------------------------------------Step3_Trim
mkdir Step3_Trim        
for i in {1..16}
do
fastp -i ./Step2_MergeR1R2/DNA${i}.fq.gz.extendedFrags.fastq.gz -o ./Step3_Trim/DNA${i}.clean.fq.gz -w 5 -q 30 -l 100 -j ./Step3_Trim/DNA${i}.clean.json -h ./Step3_Trim/DNA${i}.clean.html
fastp -i ./Step2_MergeR1R2/RNA${i}.fq.gz.extendedFrags.fastq.gz -o ./Step3_Trim/RNA${i}.clean.fq.gz -w 5 -q 30 -l 100 -j ./Step3_Trim/RNA${i}.clean.json -h ./Step3_Trim/RNA${i}.clean.html
done

#--------------------------------------------------Step4_rmdup
mkdir Step4_rmdup     
for i in {1..16}
do
seqkit rmdup --by-seq --ignore-case ./Step3_Trim/DNA${i}.clean.fq.gz | gzip > ./Step4_rmdup/DNA${i}.rmdup.fastq.gz
seqkit rmdup --by-seq --ignore-case ./Step3_Trim/RNA${i}.clean.fq.gz | gzip > ./Step4_rmdup/RNA${i}.rmdup.fastq.gz
done

#--------------------------------------------------Step5_BWA
mkdir Step5_BWA       
bwa index constant_sequence_withBarcode.fa

sbatch ./run_bwa_loop.sbatch 

#--------------------------------------------------Step6_Uniq
mkdir Step6_Uniq

for i in {1..16}
do

echo "Processing unique ${i} ..."

awk '$4 >= 0.95 { print $3}' ./Step5_BWA/DNA${i}.bwa.out > ./Step6_Uniq/DNA${i}.seq
sort ./Step6_Uniq/DNA${i}.seq | uniq > ./Step6_Uniq/DNA${i}.sort.uniq.seq

awk '$4 >= 0.95 { print $3}' ./Step5_BWA/RNA${i}.bwa.out > ./Step6_Uniq/RNA${i}.seq
sort ./Step6_Uniq/RNA${i}.seq | uniq > ./Step6_Uniq/RNA${i}.sort.uniq.seq

done

#--------------------------------------------------Step7_Cut50
mkdir Step7_Cut50

for i in {1..16}
do
awk '{print substr($0,length($0)-49,length($0))}' ./Step6_Uniq/DNA${i}.sort.uniq.seq > ./Step7_Cut50/DNA${i}.cut.seq
done

for i in {1..16}
do
awk '{print substr($0,length($0)-49,length($0))}' ./Step6_Uniq/RNA${i}.sort.uniq.seq > ./Step7_Cut50/RNA${i}.cut.seq
done

#--------------------------------------------------Step8_getBarcode
mkdir Step8_getBarcode 
for i in {1..16}
do
awk '{
    match($1, /TCTAGA[ATGC]{20}CC/, arr)
    if (match($1, /TCTAGA[ATGC]{20}CC/)) {
        print $0, arr[0]
    }
}' ./Step7_Cut50/DNA${i}.cut.seq > ./Step8_getBarcode/DNA${i}.seq
done

for i in {1..16}
do
awk '{
    match($1, /TCTAGA[ATGC]{20}CC/, arr)
    if (match($1, /TCTAGA[ATGC]{20}CC/)) {
        print $0, arr[0]
    }
}' ./Step7_Cut50/RNA${i}.cut.seq > ./Step8_getBarcode/RNA${i}.seq
done

for i in {1..16}
do
awk '{print "N" $1 "N", $2}'  ./Step8_getBarcode/DNA${i}.seq > ./Step8_getBarcode/DNA${i}.NN.seq
awk '{ gsub($2, " ", $1); print }' ./Step8_getBarcode/DNA${i}.NN.seq | awk '{print $3, $2}' > ./Step8_getBarcode/DNA${i}.NN.tmp
awk '{print $0, substr($1, length($1)-21, 20)}' ./Step8_getBarcode/DNA${i}.NN.tmp > ./Step8_getBarcode/DNA${i}.NN.Barcode.forLVseq.tmp

awk '{print $3, "CC" $2}' ./Step8_getBarcode/DNA${i}.NN.Barcode.forLVseq.tmp > ./Step8_getBarcode/DNA${i}.NN.Barcode.forLVseq
done

for i in {1..16}
do
awk '{print "N" $1 "N", $2}'  ./Step8_getBarcode/RNA${i}.seq > ./Step8_getBarcode/RNA${i}.NN.seq
awk '{ gsub($2, " ", $1); print }' ./Step8_getBarcode/RNA${i}.NN.seq | awk '{print $3, $2}' > ./Step8_getBarcode/RNA${i}.NN.tmp
awk '{print $0, substr($1, length($1)-21, 20)}' ./Step8_getBarcode/RNA${i}.NN.tmp > ./Step8_getBarcode/RNA${i}.NN.Barcode.forLVseq.tmp

awk '{print $3, "CC" $2}' ./Step8_getBarcode/RNA${i}.NN.Barcode.forLVseq.tmp > ./Step8_getBarcode/RNA${i}.NN.Barcode.forLVseq
done

#--------------------------------------------------Step9_ForLVcaculation
mkdir Step9_ForLVcaculation 
cd ./Step9_ForLVcaculation
mkdir split  
for i in {1..16}
do
split -n l/20 --numeric-suffixes=1 --additional-suffix=.split ../Step8_getBarcode/DNA${i}.NN.Barcode.forLVseq ./split/DNA${i}.Barcode.forLV
split -n l/20 --numeric-suffixes=1 --additional-suffix=.split ../Step8_getBarcode/RNA${i}.NN.Barcode.forLVseq ./split/RNA${i}.Barcode.forLV
done

R
library(stringdist)
library(data.table)

for (i in 1:16) {
  for (j in 1:20) {
    
    file <- paste0("./split/DNA", i,".Barcode.forLV",sprintf("%02d", j),".split")
    
    DNA <- data.frame(fread(file, header = FALSE))

DownArm = "CCGCTTCGAGCAGACATGAN"

DNA$DownArmLV = stringdist(DNA$V2, DownArm, method = 'lv')

DNA_LVmatch = DNA[which(DNA$DownArmLV <= 4),]
DNA_LVnomatch = DNA[which(DNA$DownArmLV > 4),]

DNA_LVmatch = DNA_LVmatch[,c('V1','DownArmLV')]
DNA_LVnomatch = DNA_LVnomatch[,c('V1','DownArmLV')]
fwrite(DNA_LVmatch,paste0('DNA',i,'.Barcode.passLV_0',j,'.out'),quo=F,row=F,col=F,sep=' ')
fwrite(DNA_LVnomatch,paste0('DNA',i,'.Barcode.failLV_0',j,'.out'),quo=F,row=F,col=F,sep=' ')
}}
q()
n

R
library(stringdist)
library(data.table)

for (i in 1:16) {
  for (j in 1:20) {
    
    file <- paste0("./split/RNA", i,".Barcode.forLV",sprintf("%02d", j),".split")
    
    RNA <- data.frame(fread(file, header = FALSE))

DownArm = "CCGCTTCGAGCAGACATGAN"

RNA$DownArmLV = stringdist(RNA$V2, DownArm, method = 'lv')

RNA_LVmatch = RNA[which(RNA$DownArmLV <= 4),]
RNA_LVnomatch = RNA[which(RNA$DownArmLV > 4),]

RNA_LVmatch = RNA_LVmatch[,c('V1','DownArmLV')]
RNA_LVnomatch = RNA_LVnomatch[,c('V1','DownArmLV')]
fwrite(RNA_LVmatch,paste0('RNA',i,'.Barcode.passLV_0',j,'.out'),quo=F,row=F,col=F,sep=' ')
fwrite(RNA_LVnomatch,paste0('RNA',i,'.Barcode.failLV_0',j,'.out'),quo=F,row=F,col=F,sep=' ')
}}
q()
n

mkdir passLV
for i in {1..16}
do
cat RNA${i}*passLV* > ./passLV/RNA${i}.passLV
cat DNA${i}*passLV* > ./passLV/DNA${i}.passLV
done
cd ..

#--------------------------------------------------Step10_OligoCount
mkdir Step10_OligoCount  

for i in {1..16}
do

awk '{print $1}' ./Step9_ForLVcaculation/passLV/DNA${i}.passLV > ./Step10_OligoCount/DNA${i}.tmp
awk '{print $1}' ./Step9_ForLVcaculation/passLV/RNA${i}.passLV > ./Step10_OligoCount/RNA${i}.tmp

sort ./Step10_OligoCount/DNA${i}.tmp | uniq -c | awk '{print $2 " " $1}' > ./Step10_OligoCount/DNA${i}.Barcode.count
sort ./Step10_OligoCount/RNA${i}.tmp | uniq -c | awk '{print $2 " " $1}' > ./Step10_OligoCount/RNA${i}.Barcode.count

join -1 1 -2 1 -t ' ' ./Step10_OligoCount/DNA${i}.Barcode.count /data1/gaoqian/MPRA/GQ_MPRA/Step11_Barcode_Specific/SpecificBarcode.sort_Oligo.strict.tsv > ./Step10_OligoCount/DNA${i}_Oligo.Barcodecount
join -1 1 -2 1 -t ' ' ./Step10_OligoCount/RNA${i}.Barcode.count /data1/gaoqian/MPRA/GQ_MPRA/Step11_Barcode_Specific/SpecificBarcode.sort_Oligo.strict.tsv > ./Step10_OligoCount/RNA${i}_Oligo.Barcodecount

awk '{print $3, $2}' ./Step10_OligoCount/DNA${i}_Oligo.Barcodecount | sort -T ./ > ./Step10_OligoCount/DNA${i}.sort.Barcodecount
awk '{sums[$1]+=$2} END {for (group in sums) {print group, sums[group]}}' ./Step10_OligoCount/DNA${i}.sort.Barcodecount > ./Step10_OligoCount/DNA${i}.sum
awk '{print $3, $2}' ./Step10_OligoCount/RNA${i}_Oligo.Barcodecount | sort -T ./ > ./Step10_OligoCount/RNA${i}.sort.Barcodecount
awk '{sums[$1]+=$2} END {for (group in sums) {print group, sums[group]}}' ./Step10_OligoCount/RNA${i}.sort.Barcodecount > ./Step10_OligoCount/RNA${i}.sum
done


# specific barcode for each oligo
cat ./Step10_OligoCount/*_Oligo.Barcodecount > ./Step10_OligoCount/MPRApoolincell.tmp
awk '{print $3, $1}' ./Step10_OligoCount/MPRApoolincell.tmp | sort | uniq > ./Step10_OligoCount/MPRApoolincell.sort.uniq.tmp
awk '{print $1}'  ./Step10_OligoCount/MPRApoolincell.sort.uniq.tmp | uniq -c | awk '{print $2 " " $1}' > ./Step10_OligoCount/Specific_Barcode_Counts_for_OligoID.MPRApoolincell
rm ./Step10_OligoCount/MPRApoolincell.sort.uniq.tmp

for i in {1..16}
do
awk '{print $3}' ./Step10_OligoCount/DNA${i}_Oligo.Barcodecount | sort | uniq -c | awk '{print $2 " " $1}' > ./Step10_OligoCount/DNA${i}.SpeciBarcodeCounts.Oligo.MPRApoolincell
awk '{print $3}' ./Step10_OligoCount/RNA${i}_Oligo.Barcodecount | sort | uniq -c | awk '{print $2 " " $1}' > ./Step10_OligoCount/RNA${i}.SpeciBarcodeCounts.Oligo.MPRApoolincell
done

for i in {1..16}
do
cat ./Step10_OligoCount/DNA${i}_Oligo.Barcodecount ./Step10_OligoCount/RNA${i}_Oligo.Barcodecount > ./Step10_OligoCount/DNARNA${i}_MPRApoolincell.tmp
awk '{print $3, $1}' ./Step10_OligoCount/DNARNA${i}_MPRApoolincell.tmp | sort -T ./ | uniq > ./Step10_OligoCount/DNARNA${i}_MPRApoolincell.sort.uniq.tmp
awk '{print $1}' ./Step10_OligoCount/DNARNA${i}_MPRApoolincell.sort.uniq.tmp | uniq -c | awk '{print $2 " " $1}' > ./Step10_OligoCount/DNARNA${i}_Specific_Barcode_Counts_for_OligoID.MPRApoolincell
rm ./Step10_OligoCount/DNARNA${i}_MPRApoolincell.sort.uniq.tmp
done
###
rm ./Step10_OligoCount/*.tmp

#--------------------------------------------------Step11_Diffexpression
mkdir Step11_Diffexpression
mv ./Step10_OligoCount/*.sum ./Step11_Diffexpression
cd ./Step11_Diffexpression

R
library(data.table)
library(readxl)
sample1 = data.frame(fread('../Step10_OligoCount/DNARNA1_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample2 = data.frame(fread('../Step10_OligoCount/DNARNA2_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample3 = data.frame(fread('../Step10_OligoCount/DNARNA3_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample4 = data.frame(fread('../Step10_OligoCount/DNARNA4_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample5 = data.frame(fread('../Step10_OligoCount/DNARNA5_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample6 = data.frame(fread('../Step10_OligoCount/DNARNA6_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample7 = data.frame(fread('../Step10_OligoCount/DNARNA7_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample8 = data.frame(fread('../Step10_OligoCount/DNARNA8_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))

sample9 = data.frame(fread('../Step10_OligoCount/DNARNA9_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample10 = data.frame(fread('../Step10_OligoCount/DNARNA10_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample11 = data.frame(fread('../Step10_OligoCount/DNARNA11_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample12 = data.frame(fread('../Step10_OligoCount/DNARNA12_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample13 = data.frame(fread('../Step10_OligoCount/DNARNA13_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample14 = data.frame(fread('../Step10_OligoCount/DNARNA14_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample15 = data.frame(fread('../Step10_OligoCount/DNARNA15_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))
sample16 = data.frame(fread('../Step10_OligoCount/DNARNA16_Specific_Barcode_Counts_for_OligoID.MPRApoolincell',header=F))

ID = read_xlsx('./oligo.xlsx')
ID = as.data.frame(ID[,c('seq_ID','group')])
setnames(ID, 'seq_ID', 'V1')

mat = merge(ID,sample1, 'V1', all.x=T)
mat = merge(mat,sample2,'V1',all.x=T)
mat = merge(mat,sample3,'V1',all.x=T)
mat = merge(mat,sample4,'V1',all.x=T)
mat = merge(mat,sample5,'V1',all.x=T)
mat = merge(mat,sample6,'V1',all.x=T)
mat = merge(mat,sample7,'V1',all.x=T)
mat = merge(mat,sample8,'V1',all.x=T)

colnames(mat) = c('oligoid','sample1','sample2','sample3','sample4','sample5','sample6','sample7','sample8')
mat$min = apply(mat[, 3:ncol(mat)], 1, min)

mat_m20 = mat[which(mat$min >=20),]
mat_l20 = mat[which(mat$min <20),]

library(tidyr)
mat_m20$variant = mat_m20$oligoid 
mat_m20 = separate(mat_m20,variant,into=c('variant','allele'),sep='_')
mat_l20$variant = mat_l20$oligoid 
mat_l20 = separate(mat_l20,variant,into=c('variant','allele'),sep='_')

mat_m20 = mat_m20[!mat_m20$variant %in% mat_l20$variant,]
fwrite(mat_m20[,c('V1','group')], './variant_passSpecificBarcodeCountinCell.oligoid')

q()
n

