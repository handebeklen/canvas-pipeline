#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bpm = 'manifest-cluster/GSACyto_20044998_A1.bpm'
params.csv = 'gtcler/GSA-24v1-0_C1.csv'
params.egt = 'manifest-cluster/2003.egt'
params.genome_fasta = '/staging/references/hg19/hg19.fa'
params.sexfile = './sexfile'
params.penncnv_pfb = 'test/out.pfb'
params.output_dir = 'outputs'
params.samplesheet = 'samplesheets/chip001.tsv'
params.idat_folder = 'idats/chip001'

workflow {
    idat_folder = Channel.fromPath(params.idat_folder)
    
    idat2gtc(idat_folder)
    gtc2vcf(idat2gtc.out)
    vcf2bedgraph(gtc2vcf.out)
    vcf2penncnv(vcf2bedgraph.out)
    penncnv_detect_all(vcf2penncnv.out)
    penncnv_clean_cnv(penncnv_detect_all.out)
    outpenncnv2bed(penncnv_clean_cnv.out)
    classification(outpenncnv2bed.out)
}

process idat2gtc {
    publishDir "${params.output_dir}/${params.idat_folder.name}/gtcs"

    input:
    path idat_folder

    output:
    path gtcs

    script:
    """
    /opt/dragena-linux-x64-DAv1.0.0/dragena/dragena genotype call --bpm-manifest ${params.bpm} --cluster-file ${params.egt} --idat-folder ${idat_folder} --output-folder ${gtcs}
    """
}

process gtc2vcf {
    publishDir "${params.output_dir}/${params.idat_folder.name}/vcfs"

    input:
    path gtcs

    output:
    path vcfs

    script:
    """
    /opt/dragena-linux-x64-DAv1.0.0/dragena/dragena genotype gtc-to-vcf --bpm-manifest ${params.bpm} --csv-manifest ${params.csv} --genome-fasta-file ${params.genome_fasta} --gtc-folder ${gtcs} --output-folder ${vcfs}
    """
}

process vcf2bedgraph {
    publishDir "${params.output_dir}/${params.idat_folder.name}/bedgraphs"

    input:
    path vcfs

    output:
    path bedgraphs

    script:
    """
    bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%REF\t%ALT\n' ${vcfs} > ${bedgraphs}
    """
}

process vcf2penncnv {
    publishDir "${params.output_dir}/${params.idat_folder.name}/beds"

    input:
    path bedgraphs

    output:
    path cnv_input

    script:
    """
    printf "Name\tChr\tPosition\tGType\tLog R Ratio\tB Allele Freq" > ${cnv_input}
    bcftools query -f '%ID\t%CHROM\t%POS[\t%GT\t%LRR\t%BAF]\n' ${bedgraphs} | sed 's/1\\/1/BB/g;s/0\\/0/AA/g;s/0\\/1/AB/g;s/.\\/./NC/g' >> ${cnv_input}
    """
}

process penncnv_detect_all {
    publishDir "${params.output_dir}/${params.idat_folder.name}/qcs"

    input:
    path cnv_input

    output:
    path combined_cnv

    script:
    """
    ./penncnv_latest.sif ./detect_cnv.pl -test -hmm lib/hhall.hmm -pfb ${params.penncnv_pfb} ${cnv_input} -log ${combined_cnv}.log -out ${combined_cnv}.out.cnv
    ./penncnv_latest.sif ./detect_cnv.pl -test -hmm lib/hhall.hmm -pfb ${params.penncnv_pfb} ${cnv_input} --confidence --chrx --sexfile ${params.sexfile} -log ${combined_cnv}.chrx.log -out ${combined_cnv}.chrx.out.cnv
    ./penncnv_latest.sif ./detect_cnv.pl -test -hmm lib/hhall.hmm -pfb ${params.penncnv_pfb} ${cnv_input} --confidence --chry --sexfile ${params.sexfile} -log ${combined_cnv}.chry.log -out ${combined_cnv}.chry.out.cnv

    cat ${combined_cnv}.out.cnv ${combined_cnv}.chrx.out.cnv ${combined_cnv}.chry.out.cnv > ${combined_cnv}.combined.out.cnv
    """
}

process penncnv_clean_cnv {
    publishDir "${params.output_dir}/${params.idat_folder.name}/beds"

    input:
    path combined_cnv

    output:
    path clean_cnv

    script:
    """
    ./penncnv_latest.sif ./clean_cnv.pl combineseg --signalfile ${params.penncnv_pfb} > ${clean_cnv}
    """
}

process outpenncnv2bed {
    publishDir "${params.output_dir}/${params.idat_folder.name}/beds"

    input:
    path clean_cnv

    output:
    path raw_cnv_bed

    script:
    """
    python3 penncnv2bed.py ${clean_cnv} ${raw_cnv_bed}
    """
}

process classification {
    publishDir "${params.output_dir}/${params.idat_folder.name}/pdfs"

    input:
    path raw_cnv_bed

    output:
    path pdfs

    script:
    """
    ./classifycnv_1.0.sif python3 /ClassifyCNV/ClassifyCNV.py --infile ${raw_cnv_bed} --GenomeBuild hg19 --precise --outfile ${pdfs}
    """
}
