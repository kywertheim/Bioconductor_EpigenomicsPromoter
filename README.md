# Bioconductor_EpigenomicsPromoter
Context: I wrote the R script named 'HW2.R' for the second assignment of the Coursera course entitled 'Bioconductor for Genomic Data Science'.

About:
1. The R script calculates the GC content values of different parts of chromosome 22, including the whole chromosome, the regions with H3K27me3 histone modifications, and its CpG islands.
2. It quantifies how the GC content, signal value, and CHIP signal are correlated in the regions with H3K27me3 histone modifications.
3. It identifies the genomic regions that gain H3K27me3 upon differentiation by spotting heightened CHIP signals in the differentiated dataset.
4. It counts the TATA boxes on chromosome 22 (both strands).
5. It counts the promoters that contain TATA boxes and transcribe the transcripts with coding sequences on chromosome 22.
6. It counts the bases forming more than one promoter of a transcript with coding sequences on chromosome 22.

Software:
1. R version 4.0.4 (2021-02-15).
2. Bioconductor version 3.12.
3. AnnotationHub.
4. GenomicRanges.
5. rtracklayer.
6. BSgenome.
