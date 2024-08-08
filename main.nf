#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bpm = 'manifest-cluster/GSACyto_20044998_A1.bpm'
params.csv = 'gtcler/GSA-24v1-0_C1.csv'
params.egt = 'manifest-cluster/2003.egt'
params.fasta = '/staging/references/hg19/hg19.fa'
params.pfb = 'test/out.pfb'
params.output_dir = 'outputs'
params.samplesheet = 'samplesheets/chip001.tsv'
params.idat_folder = 'idats/chip001'

chip_id = file(params.idat_folder).baseName
output_dir = "${params.output_dir}/${chip_id}"

workflow {
    idat_folder = Channel.fromPath(params.idat_folder)
    bpm = Channel.fromPath(params.bpm)
    csv = Channel.fromPath(params.csv)
    egt = Channel.fromPath(params.egt)
    pfb = Channel.fromPath(params.pfb)
    fasta = Channel.fromPath(params.fasta)

    fasta_dir = file(params.fasta).getParent()
    print(output_dir)

    idat2gtc(idat_folder, bpm, egt)
    vcfs = gtc2vcf(idat2gtc.out, bpm, csv, fasta_dir, fasta)



    vcfs .flatten() .map {
        vcf -> [vcf.getSimpleName(), vcf]
    }
    .set {named_vcfs}

    vcf2bedgraph(named_vcfs)
    vcf2penncnv(named_vcfs)
    penncnv_detect(vcf2penncnv.out, pfb)
    penncnv_clean_cnv(penncnv_detect.out, pfb)
    outpenncnv2bed(penncnv_clean_cnv.out)
    classification(outpenncnv2bed.out)

}

process idat2gtc {
    publishDir "${output_dir}/gtcs"

    input:
    path(idat_folder)
    path(bpm)
    path(egt)

    output:
    path "gtcs/"

    script:
    """
    mkdir gtcs/
    /opt/dragena-linux-x64-DAv1.0.0/dragena/dragena \
    genotype call \
    --bpm-manifest ${bpm} \
    --cluster-file ${egt} \
    --idat-folder ${idat_folder} \
    --output-folder gtcs/
    """
}

process gtc2vcf {
    publishDir "${output_dir}/vcfs"

    input:
    path gtcs
    path bpm
    path csv
    path fasta_dir
    path fasta

    output:
    path("*.vcf.gz")

    script:
    """
    /opt/dragena-linux-x64-DAv1.0.0/dragena/dragena genotype gtc-to-vcf \
        --gtc-folder ${gtcs} \
        --bpm-manifest ${bpm} \
        --csv-manifest ${csv} \
        --genome-fasta-file "${fasta_dir}/${fasta}" \
        --output-folder ./
    """
}

process vcf2bedgraph {
    tag "$sampleId"
    publishDir "${output_dir}/bedgraphs"

    input:
    tuple val(sampleId), path(vcf)

    output:
    path "*.bedgraph.gz"

    script:
    """
    bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%REF\t%ALT\n' ${vcf} | gzip -c > "${sampleId}.bedgraph.gz"
    """
}

process vcf2penncnv {
    tag "$sampleId"
    publishDir "${output_dir}/penncnv_inputs/"

    input:
    tuple val(sampleId), path(vcf)

    output:
    tuple val(sampleId), path("*.txt")

    script:
    """
    printf "Name\tChr\tPosition\tGType\tLog R Ratio\tB Allele Freq\n" > "${sampleId}.txt"
    bcftools query -f '%ID\t%CHROM\t%POS[\t%GT\t%LRR\t%BAF]\n' ${vcf} | sed 's/1\\/1/BB/g;s/0\\/0/AA/g;s/0\\/1/AB/g;s/.\\/./NC/g' >> "${sampleId}.txt"
    """
}

process penncnv_detect {
    container "/mnt/dragen/pipelines/canvas/penncnv_latest.sif"
    publishDir "${output_dir}/cnvs"

    input:
    tuple val(sampleId), path(txt)
    path pfb

    output:
    tuple val(sampleId), path("*.txt")

    script:
    """
    /home/user/PennCNV/detect_cnv.pl -test -hmm /home/user/PennCNV/lib/hhall.hmm -pfb ./${pfb} ${txt} --confidence         -log ${sampleId}.log      -out ${sampleId}.autosomal.out.cnv
    /home/user/PennCNV/detect_cnv.pl -test -hmm /home/user/PennCNV/lib/hhall.hmm -pfb ./${pfb} ${txt} --confidence --chrx  -log ${sampleId}.chrx.log -out ${sampleId}.chrx.out.cnv
    /home/user/PennCNV/detect_cnv.pl -test -hmm /home/user/PennCNV/lib/hhall.hmm -pfb ./${pfb} ${txt} --confidence --chry  -log ${sampleId}.chry.log -out ${sampleId}.chry.out.cnv

    cat ${sampleId}.autosomal.out.cnv ${sampleId}.chrx.out.cnv ${sampleId}.chry.out.cnv > ${sampleId}.cnv.txt
    """
}

process penncnv_clean_cnv {
    container "/mnt/dragen/pipelines/canvas/penncnv_latest.sif"
    publishDir "${output_dir}/cnvs"

    input:
    tuple val(sampleId), path(cnv)
    path pfb

    output:
    path "*.cleaned.txt"

    script:
    """
    /home/user/PennCNV/clean_cnv.pl combineseg --signalfile ${pfb} ${cnv} > "${sampleId}.cleaned.txt"
    """
}

process outpenncnv2bed {
    publishDir "${output_dir}/beds"

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
    publishDir "${output_dir}/classifyCNV"

    input:
    path raw_cnv_bed

    output:
    path pdfs

    script:
    """
    ./classifycnv_1.0.sif python3 /ClassifyCNV/ClassifyCNV.py --infile ${raw_cnv_bed} --GenomeBuild hg19 --precise --outfile ${pdfs}
    """
}
