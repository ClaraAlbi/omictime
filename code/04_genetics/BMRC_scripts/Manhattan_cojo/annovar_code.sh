# perl annovar/table_annovar.pl ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/cojo_CA.avinput annovar/humandb/ \
#     -buildver hg19 \
#     -out myanno \
#     -remove \
#     -protocol refGene \
#     -operation g \
#     -nastring . \
#     -polish \
#     -arg '-neargene 50000'

# awk -F'\t' 'BEGIN{OFS="\t"}
# {
#   chrom=$3;
#   if(substr(chrom,1,3)!="chr") chrom="chr"chrom;
#   txStart=$5; txEnd=$6; gene=$13;
#   # ensure start is non-negative
#   if(txStart<0) txStart=0;
#   print chrom, txStart, txEnd, gene
# }' annovar/humandb/hg19_refGene.txt > genes.bed

# awk 'BEGIN{OFS="\t"}{
#   chrom=$1; pos=$2;
#   if(substr(chrom,1,3)!="chr") chrom="chr"chrom;
#   start=pos-1; end=pos;
#   print chrom, start, end, $0
# }' ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/cojo_CA.avinput > cojo_CA.bed

bedtools window -a cojo_CA.bed -b code/genes.bed -w 100000 > cojo_CA_100kb.tsv