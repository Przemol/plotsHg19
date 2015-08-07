require(seqplots)
require(rtracklayer)
require(magrittr)

setwd('/Volumes/raid0/_SeqPlots_paper/SubPlotsHuman/')
source('fn_plotTopBar.R')
plotTopBar <- plotTopBar
devtools::load_all("/Users/przemol/code/seqplotsR")


tracks <- dir('tracks', pattern = 'bw', full.names = TRUE)
features <- dir('features', pattern = 'bed', full.names = TRUE)

##### Top figures #############
m1 <- getPlotSetArray(
    tracks = tracks,
    features = grep('band1', features, value = TRUE),
    refgenome = 'hg19',
    bin = 10L, xmin = 1500L, xmax = 3000L
)

pdf('a2.pdf', width = 100.0, height = 100.0, onefile = FALSE, paper = 'a4r')
plotAverage(
    m1[,c(3,1,2)], labels = c('Pol2', 'H3K36me3', 'MNaseZ'), 
    legend_pos = 'topleft', legend_ext = TRUE, legend_ext_pos = 'topright',
    cex.legend = 18, cex.main = 36, cex.lab = 28, cex.axis = 22,
    main = "Top expression quintile", xlab = "TSS", ylab = "Locally Z-scored ChIP signal",
    plotScale = "", #'zscore', 'log2'
)
dev.off()


pdf('b2.pdf', width = 100.0, height = 100.0, onefile = FALSE, paper = 'a4r')
plotAverage(
    m[c(1,3,5),2], labels = c("TSS_top_20%", "TSS_second_20%", "TSS_middle_20%", "TSS_fourth_20%", "TSS_bottom_20%")[c(1,3,5)], 
    legend_pos = 'topleft', legend_ext = TRUE, legend_ext_pos = 'topright',
    cex.legend = 17, cex.main = 36, cex.lab = 28, cex.axis = 22,
    main = "H3K4me3 ChIP-seq signal", xlab = "TSS", ylab = "Z-scored ChIP signal",
    plotScale = ""
)
dev.off()

##### Bottom figures C and d  #############
ms <- MotifSetup()
ms$addBigWig(tracks[2])$addBigWig(tracks[3])$addBigWig(tracks[1])
ms$addMotif('CG', revcomp = FALSE, name = "CpG")

M <- getPlotSetArray(
    tracks = ms,
    features = features[1],
    refgenome = 'hg19',
    bin = 10L, xmin = 1000L, xmax = 1500L
)
plotTopBar(M[,c(3,1,2,4)], 'c.pdf')


pdf('d_bw.pdf', width = 100.0, height = 100.0, onefile = FALSE, paper = 'a4r')
cls <- plotHeatmap(
    M[,c(1,3,2,4)], raster = TRUE, include = c(T,F,F,F), sortrows = 'decreasing', 
    clusters = 2,  labels = c('H3K4me3', 'H3K36me3', 'MNase160', 'CpG'), 
    cex.legend = 18, cex.main = 36, cex.lab = 28, cex.axis = 22, 
    colvec = list(NA, NA, NA, c('white', 'darkgreen', 'darkgreen')),
    #o_min = c(-1, -2, -1, 8),
    #o_max = c(3, 8, 4, 16)
)
dev.off()

##### Bottom figures D and F  #############
export.bed(import.bed(features[1])[cls$ClusterID==1], 'flt.bed')
flt <- 'flt.bed'
MF <- getPlotSetArray(
    tracks = ms,
    features = flt,
    refgenome = 'hg19',
    bin = 10L, xmin = 5000L, xmax = 5000L
)
plotTopBar(MF[,c(3,1,2,4)], 'e.pdf')

pdf(
    'f_bw.pdf', width = 100.0, height = 100.0, onefile = FALSE, paper = 'a4r'
)
plotHeatmap(
    MF, raster = TRUE, include = c(F,T,T,F), sortrows = 'decreasing',
    clstmethod = 'ssom', ssomt1 = 2, ssomt2 = 2,
    clusters = 3,  labels = c('H3K4me3', 'H3K36me3', 'MNase', 'CpG'),
    cex.legend = 18, cex.main = 36, cex.lab = 28, cex.axis = 22,
    #clspace = c('darkblue', 'white', 'darkred'),
    colvec = list(NA, NA, NA, c('white', 'darkgreen', 'darkgreen')),
    #clspace = c('white', 'grey', 'black'),
    #o_min = c(-1, -2, -1, 8),
    #o_max = c(3, 8, 4, 16)
)
dev.off()




