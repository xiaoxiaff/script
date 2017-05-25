#!/usr/bin/env Rscript
library(polyester)
library(Biostrings)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<5) {
  stop("Usage: Rscript --vanilla [script location] [number_of_transcripts] [readlen] [error_rate] [coverage] [output_dir]", call.=FALSE)
}

#args = c("10","100","0.001","20","/Users/liyuanqi/Google_Drive/UCLA_MSCS/Quarter3/CS229S/Project")

################# input parameters ############################

number_of_samples = 10

# max = 918
number_of_transcripts = as.numeric(args[1])

# default = 100
readlen = as.numeric(args[2])

# between 0~1
error_rate = as.numeric(args[3])

coverage = as.numeric(args[4])

# '/Users/liyuanqi/Google_Drive/UCLA_MSCS/Quarter3/CS229S/Project'
outdir = args[5]

print(number_of_samples)
print(number_of_transcripts)
print(readlen)
print(error_rate)
print(outdir)

###############################################################

# FASTA annotation
fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet(fasta_file)

# subset the FASTA file to first [number_of_transcripts] transcripts
small_fasta = fasta[1:number_of_transcripts]
small_fasta_file_location = paste(outdir, "/chr22_small.fa", sep="")
writeXStringSet(small_fasta, small_fasta_file_location, width=100000)

names = names(small_fasta)

lapply(names, write, paste(outdir,"/transcript_names.txt",sep=""), append=TRUE)

# coverage ----> reads per transcript = transcriptlength/readlength * coverage
# here all transcripts will have ~equal FPKM
readspertx = round(coverage * width(small_fasta) / readlen)
lapply(readspertx, write, paste(outdir,"/num_of_reads.txt",sep=""), append=TRUE)

fold_changes = matrix(c(rep(1,number_of_transcripts),rep(1,number_of_transcripts)), nrow=number_of_transcripts)

# simulation call:
simulate_experiment(small_fasta_file_location,
                    reads_per_transcript=readspertx,
                    readlen=readlen,
                    num_reps=c(number_of_samples/2,number_of_samples/2),
                    fold_changes=fold_changes,
                    outdir=paste(outdir,"/simulated_reads",sep=""),
                    paired=TRUE,
                    seed=7) 
