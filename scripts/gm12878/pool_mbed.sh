#!/bin/bash
for mod in cpg gpc; do
  echo $mod
  out=pooled/mbed/190603_GM12878_nanoNOMe.pooled.$mod.meth.bed.gz
  find . -name "GM12878*$mod*meth.bed.gz" -exec gunzip -c {} \; |\
    sort -k1,1 -k2,2n -T ./ |\
    bgzip > $out
  tabix -p bed $out

done
