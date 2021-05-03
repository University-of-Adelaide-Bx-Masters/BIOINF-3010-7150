#!/bin/bash

# Get stats with all the variants ----
echo "allele_frq_cutoff"
echo "mapping_identity_rate"

echo "0" # no cutoff in the first  time.
vg map -x z.xg -g z.gcsa -G z.sim --compare -j | jq .identity | sed s/null/0/ | awk '{i+=$1; n+=1} END {print i/n}'

# Looping through frequency cuttoffs ----
for frq in 0.01 0.1 0.5 # DEFINE FREQUENCY CUTOFFS
do 

  echo $frq

  ## SUBSAMPLE VCF FILE
  #bcftools filter -i "AF > ${frq}" z.vcf.gz > z.min_af${frq}_filtered.vcf

  ## BUILD & INDEX THE GRAPH
  #vg construct -r z.fa -v z.min_af${frq}_filtered.vcf -m 32 > z.min_af${frq}_filtered.vg 2>/dev/null
  #vg index -x z.min_af${frq}_filtered.xg z.min_af${frq}_filtered.vg
  #vg index -g z.min_af${frq}_filtered.gcsa -k 16 z.min_af${frq}_filtered.vg

  ## MAP AND GET STAT
  vg map -x z.min_af${frq}_filtered.xg -g z.min_af${frq}_filtered.gcsa -G z.sim --compare -j | jq .identity | sed s/null/0/ | awk '{i+=$1; n+=1} END {print i/n}'
  
  
done
