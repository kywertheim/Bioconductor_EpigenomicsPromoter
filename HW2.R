#Load the libraries needed for the following tasks.
library(BSgenome)
library(AnnotationHub)
library(GenomicRanges)
library(rtracklayer)

#Question 1.
#What is the GC content of "chr22" in the "hg19" build of the human genome?

#Check the available genomes in BSgenome and load the hg19 genome.
available.genomes()
library(BSgenome.Hsapiens.UCSC.hg19)

#Check the sequences in the dataset, especially those on chromosome 22.
seqnames(BSgenome.Hsapiens.UCSC.hg19)
BSgenome.Hsapiens.UCSC.hg19$chr22

#Find out the fraction of bases formed by each base type on chromosome 22.
Fraction_A <- letterFrequency(BSgenome.Hsapiens.UCSC.hg19$chr22,"A",as.prob=TRUE)
Fraction_C <- letterFrequency(BSgenome.Hsapiens.UCSC.hg19$chr22,"C",as.prob=TRUE)
Fraction_T <- letterFrequency(BSgenome.Hsapiens.UCSC.hg19$chr22,"T",as.prob=TRUE)
Fraction_G <- letterFrequency(BSgenome.Hsapiens.UCSC.hg19$chr22,"G",as.prob=TRUE)
Fraction_N <- letterFrequency(BSgenome.Hsapiens.UCSC.hg19$chr22,"N",as.prob=TRUE)

#Count the number of bases of each base type on chromosome 22.
Count_A <- letterFrequency(BSgenome.Hsapiens.UCSC.hg19$chr22,"A",as.prob=FALSE)
Count_C <- letterFrequency(BSgenome.Hsapiens.UCSC.hg19$chr22,"C",as.prob=FALSE)
Count_T <- letterFrequency(BSgenome.Hsapiens.UCSC.hg19$chr22,"T",as.prob=FALSE)
Count_G <- letterFrequency(BSgenome.Hsapiens.UCSC.hg19$chr22,"G",as.prob=FALSE)
Count_N <- letterFrequency(BSgenome.Hsapiens.UCSC.hg19$chr22,"N",as.prob=FALSE)

#Calculate the CG content of chromosome 22.
CGcontent <- (Count_C + Count_G)/(Count_A + Count_C + Count_T + Count_G)

#Question 2.
#In the previous assessment we studied H3K27me3 "narrowPeak" regions from the H1 cell line (recall that the Roadmap ID for this cell line is "E003"). We want to examine whether the GC content of the regions influence the signal; in other words whether the reported results appear biased by GC content.
#What is mean GC content of H3K27me3 "narrowPeak" regions from Epigenomics Roadmap from the H1 stem cell line on chr 22?

#Create an annotation hub.
ah <- AnnotationHub()

#Search for the Epigenomics Roadmap data on H3K27me3 histone modification for the H1 cell line.
ah_h1H3K27me3 <- query(ah, c('H3K27me3', 'H1 cell', 'roadmap'))

#Download the narrow-peak data.
ah_h1H3K27me3_narrow <- ah_h1H3K27me3[['AH29892']]

#Check that there are many chromosomes in this dataset.
seqlevels(ah_h1H3K27me3_narrow)

#Split the dataset into groups defined by chromosome.
ah_h1H3K27me3_narrow_split <- split(ah_h1H3K27me3_narrow, seqnames(ah_h1H3K27me3_narrow))

#Apply the filter to subset the regions lying on chromosome 22.
h1H3K27me3_chr22 <- unlist(ah_h1H3K27me3_narrow_split['chr22'])

#Obtain the sequences corresponding to these regions.
h1H3K27me3_chr22_seq <- Views(BSgenome.Hsapiens.UCSC.hg19, h1H3K27me3_chr22)

#Calculate the CG content of each region and average the results.
CGcontent_h1H3K27me3_chr22 <- letterFrequency(h1H3K27me3_chr22_seq, 'CG', as.prob=TRUE)
mean(CGcontent_h1H3K27me3_chr22)

#Question 3.
#What is the correlation between GC content and "signalValue" of these regions (on chr22)?

#Calculate the correlation.
cor(h1H3K27me3_chr22$signalValue, CGcontent_h1H3K27me3_chr22)

#Question 4.
#The "narrowPeak" regions are presumably reflective of a ChIP signal in these regions. To confirm this, we want to obtain the "fc.signal" data from AnnotationHub package on the same cell line and histone modification. This data represents a vector of fold-change enrichment of ChIP signal over input.
#What is the correlation between the "signalValue" of the "narrowPeak" regions and the average "fc.signal" across the same regions?

#Search for the fc.signal data and download the data in the relevant regions.
ah_h1H3K27me3_fc <- query(ah, c('H3K27me3', 'H1 cell', 'fc.signal'))
ah_h1H3K27me3_fc_regions <- ah_h1H3K27me3_fc[['AH32033']]

#Find out the length of chromosome 22 and create a GRange to represent it.
length(BSgenome.Hsapiens.UCSC.hg19$chr22)
chr22 <- GRanges(seqnames = "chr22", ranges = IRanges(start =1, end = 51304566))

#Subset the data, keeping those on chromosome 22 only.
ah_h1H3K27me3_fc_dummy <- import(ah_h1H3K27me3_fc_regions, which = chr22, as = "Rle")
ah_h1H3K27me3_fc_chr22 <- ah_h1H3K27me3_fc_dummy$chr22

#View the subsetted data by applying the ranges representing the chromosome.
fc.signal_chr22 <- Views(ah_h1H3K27me3_fc_chr22, start=start(h1H3K27me3_chr22), end=end(h1H3K27me3_chr22))

#Calculate the correlation.
cor(mean(fc.signal_chr22), h1H3K27me3_chr22$signalValue)

#Question 5.
#How many bases on chr22 have an fc.signal greater than or equal to 1?
sum(ah_h1H3K27me3_fc_chr22 >= 1)

#Question 6.
#The H1 stem cell line is an embryonic stem cell line, a so-called pluripotent cell. Many epigenetic marks change upon differentiation. We will examine this. We choose the cell type with Roadmap ID "E055" which is foreskin fibroblast primary cells.
#We will use the "fc.signal" for this cell type for the H3K27me3 mark, on chr22. We now have a signal track for E003 and a signal track for E055. We want to identify regions of the genome which gain H3K27me3 upon differentiation. These are regions which have a higher signal in E055 than in E003. To do this properly, we would need to standardize (normalize) the signal across the two samples; we will ignore this for now.
#Identify the regions of the genome where the signal in E003 is 0.5 or lower and the signal in E055 is 2 or higher.

#Search for the fc.signal data for this cell line (E055) and download the data in the relevant regions.
ah_E055H3K27me3_fc <- query(ah, c('H3K27me3', 'E055', 'fc.signal'))
ah_E055H3K27me3_fc_regions <- ah_E055H3K27me3_fc[['AH32470']]

#Subset the data, keeping those on chromosome 22 only.
ah_E055H3K27me3_fc_dummy <- import(ah_E055H3K27me3_fc_regions, which = chr22, as = "Rle")
ah_E055H3K27me3_fc_chr22 <- ah_E055H3K27me3_fc_dummy$chr22

#Identify the regions in E003 where the signal is at most 0.5 and the regions in E055 where the signal is at least 2.
#Present the results as IRanges rather than Views.
WeakFc_E003 <- as(slice(ah_h1H3K27me3_fc_chr22, upper = 0.5), 'IRanges')
StrongFc_E055 <- as(slice(ah_E055H3K27me3_fc_chr22, lower = 2), 'IRanges')

#Find the regions where the signal is higher in E055 than in E003 and then count the number of bases.
sum(width(intersect(WeakFc_E003, StrongFc_E055)))

#Question 7.
#CpG Islands are dense clusters of CpGs. The classic definition of a CpG Island compares the observed to the expected frequencies of CpG dinucleotides as well as the GC content.
#Specifically, the observed CpG frequency is just the number of "CG" dinucleotides in a region. The expected CpG frequency is defined as the frequency of C multiplied by the frequency of G divided by the length of the region.
#What is the average observed-to-expected ratio of CpG dinucleotides for CpG Islands on chromosome 22?

#Search for the CpG islands in the hg19 genome and download the data.
ah_hg19_CpG <- query(ah, c("hg19", "CpG"))
ah_hg19_CpG_regions <- ah_hg19_CpG[['AH5086']]

#Subset the data, keeping those on chromosome 22 only and download the sequences in these regions.
ah_hg19chr22_CpG_regions <- subset(ah_hg19_CpG_regions, seqnames == 'chr22')
ah_hg19chr22_CpG_seq <- Views(BSgenome.Hsapiens.UCSC.hg19, ah_hg19chr22_CpG_regions)

#Calculate the average observed-to-expected ratio of CpG dinucleotides in the sequences.
CG_observed <- dinucleotideFrequency(ah_hg19chr22_CpG_seq)[,7]
CG_expected <- letterFrequency(ah_hg19chr22_CpG_seq, 'C')*letterFrequency(ah_hg19chr22_CpG_seq, 'G')/width(ah_hg19chr22_CpG_seq)
mean(CG_observed/CG_expected)

#Question 8.
#A TATA box is a DNA element of the form "TATAAA". Around 25% of genes should have a TATA box in their promoter. We will examine this statement.
#How many TATA boxes are there on chr 22 of build hg19 of the human genome?

#Search for the pattern in the forward strand.
TATAbox_forward <- countPattern("TATAAA", BSgenome.Hsapiens.UCSC.hg19$chr22)

#Search for the pattern in the reverse strand.
TATAbox_reverse <- countPattern("TATAAA", reverseComplement(BSgenome.Hsapiens.UCSC.hg19$chr22))

#Add the two counts.
TATAbox_forward + TATAbox_reverse

#Question 9.
#How many promoters of transcripts on chromosome 22 containing a coding sequence, contains a TATA box on the same strand as the transcript?
#I failed to obtain the answer expected.

#Load the genomic features of hg19.
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Locate the coding sequences and keep those on chromosome 22 only.
codingreg <- cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = 'tx')
codingreg_chr22 <- keepSeqlevels(codingreg, 'chr22', pruning.mode = 'coarse')

#Find the transcripts on chromosome 22.
transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
transcripts_chr22 <- subsetByOverlaps(transcripts, chr22, ignore.strand = TRUE)

#Find the promoters of these transcripts 
promoters_chr22 <- promoters(transcripts_chr22, upstream=900, downstream = 100)

#Keep the promoters containing coding sequences only.
promoters_chr22_seq <- findOverlaps(promoters_chr22, codingreg_chr22)

#Count the number of promoters containing at least one TATA box.
#Make sure the promoters are not counted more than once as one promoter can contain multiple coding regions.
count <- 0

for (i in unique(queryHits(promoters_chr22_seq)))
{
dummy <- Views(BSgenome.Hsapiens.UCSC.hg19, promoters_chr22[i])
if (vcountPattern("TATAAA", DNAStringSet(dummy))>0)
  {
  count <- count + 1
  }
}

#Question 10.
#It is possible for two promoters from different transcripts to overlap, in which case the regulatory features inside the overlap might affect both transcripts. This happens frequently in bacteria.
#How many bases on chr22 are part of more than one promoter of a coding sequence?

#Find out the lengths of the transcripts associated with hg19.
transcript_lengths <- transcriptLengths(TxDb.Hsapiens.UCSC.hg19.knownGene,  with.cds_len = TRUE)

#Only retain the transcripts with coding sequences.
transcript_lengths <- transcript_lengths[transcript_lengths$cds_len > 0, ]

#Identify the promoters of the subsetted transcripts on chromosome 22.
promoters_chr22_seq2 <- promoters_chr22[mcols(promoters_chr22)$tx_id %in% transcript_lengths$tx_id]

#Find out the coverage of each base and identify those with more than one read, i.e. such a base is a part of more than one promoter of a coding sequence.
#Then, sum over each chromosome.
sum(coverage(promoters_chr22_seq2)>1)