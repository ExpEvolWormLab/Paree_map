require(data.table)
library("R.utils")
library('parallel')

infastq = as.vector(list.files(path="/mnt/data1/paree/Data/2bRadseq/2bRadseq2/output_stacks", full.names = T, pattern = "*.fq"))
library(readxl)
#cross =read_excel("Documents/Data/2bRadSeq/barcode/2bRad_samples.xlsx")
#x = which(cross$Sample %in% c('C', 'D', 'G', 'H'))
#infastq=infastq[x]
outfastq =unlist(lapply(infastq, function(i) strsplit(i, ".fq")))
outfastq =unlist(lapply(outfastq, function(i) strsplit(i, "output_stacks/")))
outfastq = paste(outfastq[seq(1,length(outfastq),2)], outfastq[seq(2,length(outfastq),2)], sep = "trimmed" )
outfastq = paste(outfastq, rep("trimmed.fq", length(outfastq)), sep = "_" )

np = 20

# read =
#1:@@711_1_11101_21321_1137/1
#2: TNCAAAGAGTTTGCAGGTTCCTCGTGATTGTGCATANGATCGGAAGAGCN ...
#3: +
#4: AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA...

#output of btrim: read if perfect BcgI match, else number of mismatch

btrim <- function(read) {
  # pass a single fastq read block
  read = as.character(unlist(read))
  # get rid of 3' adaptor
  read[2] = substr(read[2], 1, 36)
  read[4] = substr(read[4], 1, 36)
  sp = strsplit(read[2], '')[[1]]
  # weak primer nucleotides at 36bp fragment +2 & -2
  # site is NW.10N.CGA.6N.TGC.10N.WN | NW.10N.GCA.6N.TCG.10N.WN
  wl = sp[2] %in% c('A', 'T')
  wr = sp[35] %in% c('A', 'T')
  if(wl & wr){
    checkf = c(sp[13:15]==c('C', 'G', 'A'), sp[22:24]==c('T', 'G', 'C'))
    checkr = c(sp[13:15]==c('G', 'C', 'A'), sp[22:24]==c('T', 'C', 'G'))
    if(all(checkf) | all(checkr)){
      # return(T)
      return(read)
    } else {
      # Hamming distance from recognition sites
      return(6-max(c(sum(checkf), sum(checkr))))
      # cat(c(as.numeric(wl), as.numeric(wr), paste(sp[13:15], collapse=''), paste(sp[22:24], collapse='')), '\n')
    }
  }
  return(F)
}


for(i in 1:length(infastq)) {
df = fread(infastq[i], header=F)
nr = nrow(df)/4
rfilter = mclapply(seq(1, nrow(df), 4), mc.cores = np, function(i) btrim(df[i:(i+3),]))
rm(df)
# failed first due to one or both weak sites missing in adapters
weak_missing = sum(sapply(rfilter, function(x) length(x)==1 & x[1]==F))
# weak sites present, but imperfect recognition sites
site_missing = length(unlist(sapply(rfilter, function(x) if(length(x)==1 & x[1]>0) x)))
one_site_missing = sum((unlist(sapply(rfilter, function(x) if(length(x)==1 & x[1]>0) x)))==1)
# pass
goodr = sapply(rfilter, function(x) length(x)==4)
goodr = rfilter[goodr]
rm(rfilter)
goodo = do.call(c, goodr)
writeLines(goodo, outfastq[i])

sink("stat_trim2.txt", append = T)
cat(sprintf("%s\n%s reads\n%.2f pass\n%s failed initial weak site check\n%s pass but have imperfect RE sites (%s within Hamming 1)\n",
            basename(infastq[i]), nr, length(goodr)/nr, weak_missing, site_missing, one_site_missing))
cat("\n")
sink()

rm(goodr)
rm(goodo)
gc()

}
s

