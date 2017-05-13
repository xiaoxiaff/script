library(polyester)
library(Biostrings)

outdir = '/Users/liyuanqi/Google\ Drive/UCLA_MSCS/Quarter3/CS229S/Project/simulated_reads'

# FASTA annotation
fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet(fasta_file)


# subset the FASTA file to first 20 transcripts
small_fasta = fasta[1:20]
writeXStringSet(small_fasta, 'chr22_small.fa')

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
# readspertx = round(20 * width(small_fasta) / 100)
readspertx = round(20 * width(small_fasta) / 100)

fold_changes = matrix(c(rep(1,20),rep(1,20)), nrow=20)

# simulation call:
simulate_experiment('chr22_small.fa', reads_per_transcript=readspertx, 
                    num_reps=c(10,10), fold_changes=fold_changes, outdir=outdir, paired=TRUE) 
