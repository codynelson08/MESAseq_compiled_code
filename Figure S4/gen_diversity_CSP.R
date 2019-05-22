#calculate Pi, PiS, PiN for csp data

load('diversity_piNS_3.RData')

##  load data. 
haps = read.csv('csp_haps_seqs_demog_merged.csv')
haps = haps [,2:ncol(haps)]

##  Get the sequences in a DNAStringSet.
seqs = names(haps)
seqs = seqs[nchar(seqs) > 100]
##  Convert to another format. 
library('Biostrings')
library(ape)
seqs = DNAStringSet(seqs)
seqs


## Create haplotypes ids. ->  as column names +  ids in the dna string. 
haps_ids = paste('H', 1:length(seqs), sep='')

names(seqs) = haps_ids #  Change names stringset. 
##  Change names of the columns. 
names(haps)[1:length(seqs)] = haps_ids

##  convert the sequences. 
seqs2 = reverseComplement(seqs)
seqs2 = subseq(seqs2, 3, 287)
translate(seqs2)

##  Calculate genetic diversity.
results = read.csv("csp_haps_demog_merged_nocontrols.csv")
results = results [,2:124]
hap_freq = as.matrix(results[,haps_ids])
results$piS = results$piN = results$pi = NA

for(i in 1:nrow(results)){
  ##  Per row. 
  tmp = hap_freq[i,]
  tmp
  results[i, ] %>% class
  ix = which( !tmp ==0 )
  freq = as.vector(tmp[ix]) * 10000  #  Hap frequencies. 
  seq = seqs2[ix] # Sequences
  #pi = nt_diversity_hap(seq, freq)
  pi = calculate_piNS(seq, freq)
  results$pi[i] = pi[1]
  results$piS[i] = pi[2]
  results$piN[i] = pi[3]

}

write.csv(results, 'csp_haps_demog_merged_nocontrols_pi.csv')






