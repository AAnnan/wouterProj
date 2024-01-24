# Suppress library loading messages
suppressPackageStartupMessages({
  library("GenomicRanges")
  library("data.table")
  library("parallel")
})
options(scipen = 999) #makes sure numbers are not in sci notation
#numCores <- detectCores() #use all cores

### Arguments
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check the arguments provided
if (length(args) != 5) {
  stop("Usage: Rscript IOM.R <directory> <peakWidth> <inBed> <outBed> <numCores>")
}

# Set
#working directory
wd <- args[1]
setwd(wd)
#padding around summit
peakWidth <- as.numeric(args[2])
#BED summit file from MACS
inBed <- args[3]
#output file names
outBed <- args[4]
#number of cores
numCores <- args[5]

###Import MACS summit bed file
summits = fread(file=inBed, header=F)
colnames(summits) = c('chr','start','end','name','score')
gr_summits = GRanges(summits) + peakWidth/2

###Parallel iterative overlap peak merging
ClusteredOverlap = function(gr_s){
  #Get list of overlapping region
  gr_summits_revmap = mcols(reduce(gr_s,with.revmap=T))$revmap
  
  #From a set of GRanges, return the one with highest score
  keep_max_score <- function(idx_rev,gr) {
    if (length(idx_rev)==1) {return(gr[idx_rev])} 
    else {return(gr[idx_rev][which.max(gr[idx_rev]$score)])}
  }
  
  #Loop overlapping GRanges
  l_gr_summits = lapply(gr_summits_revmap,keep_max_score,gr=gr_s)
  
  #Concatenate
  red_gr_summits = do.call("c", l_gr_summits)
  rm(l_gr_summits)
  
  return(red_gr_summits)
}
IterativeOverlap = function(gr_s){
  if (length(gr_s)==1) {return(gr_s)}
  
  #ini
  n=1
  gr_res <- GRanges(seqnames = character(),
                    ranges = IRanges(start = integer(),
                                     end = integer()),
                    strand = character(),
                    name = character(),
                    score = numeric())
  
  #Sort by decreasing score
  gr_s = gr_s[order(-mcols(gr_s)$score)]
  
  #Go through peak list
  while (length(gr_s)) {
    
    #First scoring peak is always kept for the final list
    # and removed from the working list
    peak = gr_s[1]
    gr_res = do.call("c", list(gr_res, peak))
    gr_s = gr_s[-1]
    
    #get overlappign peaks
    overlap = subjectHits(findOverlaps(peak,gr_s,ignore.strand=T))
    
    #If no overlap continue
    if (!length(overlap)) {
      n = n+1
      next}
    
    #If overlap remove it 
    gr_s = gr_s[-overlap]
    n = n+1
  }
  return(gr_res)
}
ClItOverlapParallel = function(gr_s) {
  # Get list of overlapping region
  gr_summits_revmap = mcols(reduce(gr_s, with.revmap = TRUE))$revmap
  
  # From a set of overlapping GRanges, return the result of iterative overlap
  keep_max_score <- function(idx_rev, gr) {
    return(IterativeOverlap(gr[idx_rev]))
  }
  
  # Parallelize the loop using mclapply
  l_gr_summits = mclapply(gr_summits_revmap, keep_max_score, gr = gr_s, mc.cores = numCores)
  
  # Concatenate
  red_gr_summits = do.call("c", l_gr_summits)
  rm(l_gr_summits)
  
  return(red_gr_summits)
}
red_gr_summits = ClItOverlapParallel(gr_summits)

### Get Blacklist region
#Nature Paper (Amamiya, Kundaje, Boyle et al 2019)
#https://www.nature.com/articles/s41598-019-45839-z
#Download from GitHub (link from the paper)
blBedURL <- "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz"
blBed <- basename(blBedURL)
#Download and extract blacklist file
download.file(blBedURL, blBed, mode = "wb", quiet = TRUE)
system(sprintf("gzip -d %s", blBed),ignore.stdout=T)
blBed <- substr(blBed, 1, nchar(blBed) - 3)
#Import BL file
bl = fread(file=blBed, header=F)
colnames(bl) = c('chr','start','end','reason')
#Remove BL regions from Final Peak BED
gr_bl = GRanges(bl)
red_gr_summits_filt <- subsetByOverlaps(red_gr_summits, gr_bl, invert = TRUE)

#Write Final Peak BED
fwrite(as.data.frame(red_gr_summits_filt)[c(1,2,3,6,7)], file = outBed,
       sep = "\t", quote = F, row.names = F, col.names = F)

#Build Peak SAF
red_summits = as.data.frame(fread(file=outBed, header=F))[c(4,1,2,3,5)]
colnames(red_summits) = c('GeneID','Chr','Start','End','Score')
red_summits$Strand = rep('.',nrow(red_summits))
red_summits$End = red_summits$End - 1
red_summits = red_summits[c(1,2,3,4,6,5)]
#Write Peak SAF
fwrite(red_summits, file = paste0(substr(outBed,1,nchar(outBed)-3),'saf'),
       sep = "\t", quote = F, row.names = F, col.names = T)

#remove the blacklist file
file.remove(blBed)
#Done
cat("IOM done:", inBed, "\n")
