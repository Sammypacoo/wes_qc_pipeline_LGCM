zcat results/cobertura.thresholds.bed.gz | \
awk 'BEGIN {OFS="\t"; print "#chrom","start","end","region","10X","30X","%10X","%30X"} 
     {
       len = $3 - $2
       pct10 = (len > 0) ? ($5 / len) * 100 : 0
       pct30 = (len > 0) ? ($6 / len) * 100 : 0
       print $1, $2, $3, $4, $5, $6, pct10, pct30
     }' > results/cobertura.thresholds_percentual.tsv