idat2gtc:

/opt/dragena-linux-x64-DAv1.0.0/dragena/dragena genotype call --bpm-manifest manifest-cluster/GSACyto_20044998_A1.bpm --cluster-file manifest-cluster/2003.egt --idat-folder idats/ --output-folder output_idat2gtc/

gtc2vcf:

/opt/dragena-linux-x64-DAv1.0.0/dragena/dragena genotype gtc-to-vcf --bpm-manifest gtcler/GSA-24v1-0_C1.bpm --csv-manifest gtcler/GSA-24v1-0_C1.csv --genome-fasta-file /staging/references/hg19/hg19.fa --gtc-folder output_idat2gtc/ --output-folder /output/vcf

vcf2bedgraph:
bcftool query ........


vcf2penncnv:

printf "Name\tChr\tPosition\t$sample.GType\t$sample.Log R Ratio\t$sample.B Allele Freq"> $sample.txt; bcftools query -s $sample -f '%ID\t%CHROM\t%POS[\t%GT\t%LRR\t%BAF]\n' merge.vcf.gz | sed 's/1\/1/BB/g;s/0\/0/AA/g;s/0\/1/AB/g;s/.\/./NC/g' >> $sample.txt

penncnv detect cnv:

./penncnv_latest.sif ./detect_cnv.pl -test -hmm lib/hhall.hmm -pfb test/out.pfb test/208056550034_R11C02.txt -log logfile -out test/208056550034_R11C02.out.cnv

detect X

./penncnv_latest.sif ./detect_cnv.pl -test -hmm lib/hhall.hmm -pfb test/out.pfb test/208056550004_R04C02.txt --confidence --chrx --sexfile ./sexfile -log logfile -out test/208056550004_R04C02.out.chrx.cnv

detect Y:

./penncnv_latest.sif ./detect_cnv.pl -test -hmm lib/hhall.hmm -pfb test/out.pfb test/208056550004_R04C02.txt --confidence --chry --sexfile ./sexfile -log logfile -out test/208056550004_R04C02.out.chry.cnv

penncnv clean cnv:

./penncnv_latest.sif ./clean_cnv.pl combineseg --signalfile test/out.pfb > test/test42/out.cleanedcnv


outpenncnv2bed:

pyhton3 penncnv2bed.py out.rawcnv out.rawcnv.bed


classification:
./classifycnv_1.0.sif python3 /ClassifyCNV/ClassifyCNV.py --infile out.rawcnv.bed --GenomeBuild hg19 --precise
