wget https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/bigZips/genes/galGal6.ensGene.gtf.gz
gunzip galGal6.ensGene.gtf.gz

python integrated_quantification_sam_only.py \
    -g galGal6.ensGene.gtf \
    -t galGal6_TE.gtf \
    -l SRR13267657.sam \
    -o output
