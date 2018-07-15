QualityTrim <- function() {

}

myFilterAndTrim <- function(fl, destination=sprintf("%s_subset", fl)) +{
  ## open input stream
  fl <- "/home/kuan-hao/Documents/bioinformatics/howard/gene_data/raw_fastq.gz/ERR188044_1.fastq.gz"
  rfq <- readFastq(fl)
  quality(rfq)
  alphabet(quality(rfq))
  rfq1 = trimEnds(rfq, "1")
  quality(rfq1)
  repeat {
    ## input chunk
    fq <- yield(stream)
    if (length(fq) == 0)
      break
    ## trim and filter, e.g., reads cannot contain 'N'...
    fq <- fq[nFilter()(fq)]  # see ?srFilter for pre-defined filters
    ## trim as soon as 2 of 5 nucleotides has quality encoding less
    ## than "4" (phred score 20)
    fq <- trimTailw(fq, 2, "4", 2)
    ## drop reads that are less than 36nt
    fq <- fq[width(fq) >= 36]
    ## append to destination
    writeFastq(fq, destination, "a")
  }
}
