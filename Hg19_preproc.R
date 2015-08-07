setwd('/Volumes/raid0/_SeqPlots_paper/SubPlotsHuman/raw')

k36 <- dir(pattern = '36me3')
k36r1  <- import.bw(BigWigFileList(k36)[[1]], as = 'RleList')
k36r2  <- import.bw(BigWigFileList(k36)[[2]], as = 'RleList')
k36r3  <- import.bw(BigWigFileList(k36)[[3]], as = 'RleList')


k36a <- k36r1 + k36r2 + k36r3
k36a <- round(k36a, 2)
export.bw(k36a, 'H3k36me3_comb.bw')

# gr <- import.bw('H3K36me3.bw')
gr <- import.bw('H3k36me3_comb.bw')
mi <- weighted.mean(gr$score, width(gr)); s <- sqrt(sum(width(gr) * (gr$score - mi)^2) / sum(as.numeric(width(gr))) ); gr$score <- ( gr$score - mi ) / s;
gr$score <- round(gr$score, 2)


export.bw(gr, 'H3K36me3_zsc.bw')


# mi <- weighted.mean(mean(k36a), elementLengths(k36a)) 
# s <- weighted.mean(sd(k36a), elementLengths(k36a)) 
# k36zsc <- (k36a - mi) / s;
# export.bw(k36zsc, 'H3k36me3_zsc.bw')


# export.bw(gr, 'H3K36me3.bw')

k4 <- dir(pattern = '4me3')
k4r1  <- import.bw(BigWigFileList(k4)[[1]], as='RleList')
k4r2  <- import.bw(BigWigFileList(k4)[[2]], as='RleList')

k4a <- k4r1 + k4r2
k4a <- round(k4a, 2)


export.bw(k4a, 'H3K4me3_comb.bw')

 mi <- weighted.mean(mean(k4a), elementLengths(k4a)) 
s <- weighted.mean(sd(k4a), elementLengths(k4a)) 
k36zsc <- (k4a - mi) / s;
export.bw(k4a, 'k4a_zsc.bw')

gr <- import.bw('H3K4me3_comb.bw')
mi <- weighted.mean(gr$score, width(gr)); s <- sqrt(sum(width(gr) * (gr$score - mi)^2) / sum(as.numeric(width(gr))) ); gr$score <- ( gr$score - mi ) / s;
gr$score <- round(gr$score, 2)


export.bw(gr, 'H3K36me3_zsc.bw')





gr <- import.bw("/Volumes/raid0/_SeqPlots_paper/SubPlotsHuman/tracks/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bw")
mi <- weighted.mean(gr$score, width(gr)); s <- sqrt(sum(width(gr) * (gr$score - mi)^2) / sum(as.numeric(width(gr))) ); gr$score <- ( gr$score - mi ) / s;
gr$score <- round(gr$score, 2)

export.bw(gr, 'MNase_zsc.bw')


###### Comand line ######

bigWigMerge GSM*3k36* merge.bedGraph 
bedGraphToBigWig merge.bedGraph hg19.chrom.sizes merge.bw
mi=$(bigWigInfo merge.bw | grep mean | sed -n 's/mean: //pg')
sd=$(bigWigInfo merge.bw | grep std | sed -n 's/std: //pg')
cat merge.bedGraph | awk -v mi=$mi -v sd=$sd '{$4=($4-mi)/sd;print}' > zsc.bedGraph    
bedGraphToBigWig zsc.bedGraph hg19.chrom.sizes H3k36me3_bigWigMerge_zsc.bw
rm merge* zsc.bedGraph

bigWigMerge GSM*3k4* merge.bedGraph 
bedGraphToBigWig merge.bedGraph hg19.chrom.sizes merge.bw
mi=$(bigWigInfo merge.bw | grep mean | sed -n 's/mean: //pg')
sd=$(bigWigInfo merge.bw | grep std | sed -n 's/std: //pg')
cat merge.bedGraph | awk -v mi=$mi -v sd=$sd '{$4=($4-mi)/sd;print}' > zsc.bedGraph    
bedGraphToBigWig zsc.bedGraph hg19.chrom.sizes H3k4me3_bigWigMerge_zsc.bw
rm merge* zsc.bedGraph


mi=$(bigWigInfo GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bw | grep mean | sed -n 's/mean: //pg')
sd=$(bigWigInfo GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bw | grep std | sed -n 's/std: //pg')
bigWigToBedGraph GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bw stdout | awk -v mi=$mi -v sd=$sd '{$4=($4-mi)/sd;print}' > zsc.bedGraph    
bedGraphToBigWig zsc.bedGraph hg19.chrom.sizes MNase_bigWigMerge_zsc.bw
rm merge* zsc.bedGraph

f=GSM822270_hg19_wgEncodeOpenChromChipGm12878Pol2Sig.bigWig
mi=$(bigWigInfo $f | grep mean | sed -n 's/mean: //pg')
sd=$(bigWigInfo $f | grep std | sed -n 's/std: //pg')
bigWigToBedGraph $f stdout | awk -v mi=$mi -v sd=$sd '{$4=($4-mi)/sd;print}' > zsc.bedGraph    
bedGraphToBigWig zsc.bedGraph hg19.chrom.sizes Pol2_zsc.bw

bigWigMerge GSM803485_hg19_wgEncodeHaibTfbsGm12878Pol24h8Pcr1xRawRep2.bigWig \
GSM803485_hg19_wgEncodeHaibTfbsGm12878Pol24h8Pcr1xRawRep1.bigWig \
GSM803355_hg19_wgEncodeHaibTfbsGm12878Pol2Pcr2xRawRep2.bigWig \
GSM803355_hg19_wgEncodeHaibTfbsGm12878Pol2Pcr2xRawRep1.bigWig merge.bedGraph 
bedGraphToBigWig merge.bedGraph hg19.chrom.sizes merge.bw
mi=$(bigWigInfo merge.bw | grep mean | sed -n 's/mean: //pg')
sd=$(bigWigInfo merge.bw | grep std | sed -n 's/std: //pg')
cat merge.bedGraph | awk -v mi=$mi -v sd=$sd '{$4=($4-mi)/sd;print}' > zsc.bedGraph    
bedGraphToBigWig zsc.bedGraph hg19.chrom.sizes Pol2_bigWigMerge_zsc.bw
rm merge* zsc.bedGraph



f=GSM733677_hg19_wgEncodeBroadHistoneGm12878H3k9acStdSig.bigWig
mi=$(bigWigInfo $f | grep mean | sed -n 's/mean: //pg')
sd=$(bigWigInfo $f | grep std | sed -n 's/std: //pg')
bigWigToBedGraph $f stdout | awk -v mi=$mi -v sd=$sd '{$4=($4-mi)/sd;print}' > zsc.bedGraph    
bedGraphToBigWig zsc.bedGraph hg19.chrom.sizes H3K9ac_zsc.bw
rm merge* zsc.bedGraph

f=GSM733767_hg19_wgEncodeBroadHistoneGm12878H2azStdSig.bigWig
mi=$(bigWigInfo $f | grep mean | sed -n 's/mean: //pg')
sd=$(bigWigInfo $f | grep std | sed -n 's/std: //pg')
bigWigToBedGraph $f stdout | awk -v mi=$mi -v sd=$sd '{$4=($4-mi)/sd;print}' > zsc.bedGraph    
bedGraphToBigWig zsc.bedGraph hg19.chrom.sizes H2A.Z_zsc.bw
rm merge* zsc.bedGraph

f=GSM733771_hg19_wgEncodeBroadHistoneGm12878H3k27acStdSig.bigWig
mi=$(bigWigInfo $f | grep mean | sed -n 's/mean: //pg')
sd=$(bigWigInfo $f | grep std | sed -n 's/std: //pg')
bigWigToBedGraph $f stdout | awk -v mi=$mi -v sd=$sd '{$4=($4-mi)/sd;print}' > zsc.bedGraph    
bedGraphToBigWig zsc.bedGraph hg19.chrom.sizes H3K27ac_zsc.bw
rm merge* zsc.bedGraph
 


bigWigMerge -threshold=-1e10  *.bw stdout 

| awk '{$4=$4/2;print}' | wigToBigWig -clip stdin <(bigWigInfo -chroms $(ls *bw | grep -v average | head -1) | grep -P "\t" | awk '{print $1"\t"$3}') stdout > average.bw

