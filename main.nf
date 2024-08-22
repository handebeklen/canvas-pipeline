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
params.band = "canvas-pipeline/"
params.tex_template = "asdadfa"

chip_id = file(params.idat_folder).baseName
output_dir = "${params.output_dir}/${chip_id}"

workflow {
    idat_folder = Channel.fromPath(params.idat_folder)
    bpm = Channel.fromPath(params.bpm)
    csv = Channel.fromPath(params.csv)
    egt = Channel.fromPath(params.egt)
    pfb = file(params.pfb)
    fasta = Channel.fromPath(params.fasta)
    band = Channel.fromPath(params.band)
    tex_template = Channel.fromPath(params.tex_template)

    fasta_dir = file(params.fasta).getParent()
    template_dir = file(params.tex_template).getParent()

    idat2gtc(idat_folder, bpm, egt)
    makesexfile(idat2gtc.out.sample_summary)
    vcfs = gtc2vcf(idat2gtc.out.gtc_dir, bpm, csv, fasta_dir, fasta)

    vcfs .flatten() .map {
        vcf -> [vcf.getSimpleName(), vcf]
    }
    .set {named_vcfs}

    vcf2penncnv(named_vcfs)

    makebedgraphs(vcf2penncnv.out)

    penncnv_detect(vcf2penncnv.out, pfb)
    penncnv_clean_cnv(penncnv_detect.out, pfb)

    beds = penncnv2bed(penncnv_clean_cnv.out, makesexfile.out.first())
    classification(penncnv2bed.out)

    addiscn(penncnv_clean_cnv.out, band.first())

    bedgraphswithbeds = makebedgraphs.out.combine(beds, by:0)
    makeplots(bedgraphswithbeds)

    cnvswithscores = classification.out.scoresheet
        .combine(addiscn.out, by:0)
        .combine(makeplots.out, by:0)

    maketemplate(cnvswithscores, tex_template)
    makereport(maketemplate.out, template_dir)
}


process makesexfile {
    memory "1 GB"
    cpus 1
    input:
    path sample_summary

    output:
    path "sex_file"

    shell:
    '''
    awk -F"," '{printf "%s\\t%s\\n", $1, $7}' !{sample_summary} > sex_file
    '''
}

process idat2gtc {
    publishDir "${output_dir}/", mode: "copy"

    input:
    path(idat_folder)
    path(bpm)
    path(egt)

    output:
    path "gtcs/", emit: gtc_dir
    path "gtcs/gt_sample_summary.csv", emit: sample_summary

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
    publishDir "${output_dir}/vcfs", mode: "copy"

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

process vcf2penncnv {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/penncnv_inputs/", mode: "copy"

    input:
    tuple val(sampleId), path(vcf)

    output:
    tuple val(sampleId), path("${sampleId}.txt")

    script:
    """
    printf "Name\tChr\tPosition\tGType\tLog R Ratio\tB Allele Freq\n" > "${sampleId}.txt"
    bcftools query -f '%ID\t%CHROM\t%POS[\t%GT\t%LRR\t%BAF]\n' ${vcf} | sed 's/1\\/1/BB/g;s/0\\/0/AA/g;s/0\\/1/AB/g;s/.\\/./NC/g' >> "${sampleId}.txt"
    """
}

process makebedgraphs {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/bedgraphs", mode: "copy"

    input:
    tuple val(sampleId), path(txt)

    output:
    tuple val(sampleId), path("*BAF.bedgraph.gz"), path("*LRR.bedgraph.gz")


    shell:
    '''
    awk 'BEGIN {print "track type=bedGraph name="!{sampleId} LRR" maxHeightPixels=200:200:200 graphType=points viewLimits=-2:2 windowingFunction=none color=0,0,0 altColor=0,0,0"}{printf "%s\\t%s\\t%s\\t%s\\n", $2, $3, $3, $5}' !{txt} | sed 1d | sed '/\\.$/d' | gzip -c > !{sampleId}.LRR.bedgraph.gz
    awk 'BEGIN {print "track type=bedGraph name="!{sampleId} BAF" maxHeightPixels=200:200:200 graphType=points viewLimits=-2:2 windowingFunction=none color=0,0,0 altColor=0,0,0"}{printf "%s\\t%s\\t%s\\t%s\\n", $2, $3, $3, $6}' !{txt} | sed 1d | sed '/\\.$/d' | gzip -c > !{sampleId}.BAF.bedgraph.gz
    '''
}

process penncnv_detect {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    container "/mnt/dragen/pipelines/canvas/penncnv_latest.sif"
    publishDir "${output_dir}/cnvs", mode: "copy"

    input:
    tuple val(sampleId), path(txt)
    path pfb

    output:
    tuple val(sampleId), path("${sampleId}.cnv.txt")

    script:
    """
    /home/user/PennCNV/detect_cnv.pl -test -hmm /home/user/PennCNV/lib/hhall.hmm -pfb ./${pfb} ${txt} --confidence         -log ${sampleId}.log      -out ${sampleId}.autosomal.out.cnv
    /home/user/PennCNV/detect_cnv.pl -test -hmm /home/user/PennCNV/lib/hhall.hmm -pfb ./${pfb} ${txt} --confidence --chrx  -log ${sampleId}.chrx.log -out ${sampleId}.chrx.out.cnv
    /home/user/PennCNV/detect_cnv.pl -test -hmm /home/user/PennCNV/lib/hhall.hmm -pfb ./${pfb} ${txt} --confidence --chry  -log ${sampleId}.chry.log -out ${sampleId}.chry.out.cnv

    cat ${sampleId}.autosomal.out.cnv ${sampleId}.chrx.out.cnv ${sampleId}.chry.out.cnv > ${sampleId}.cnv.txt
    """
}

process penncnv_clean_cnv {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    container "/mnt/dragen/pipelines/canvas/penncnv_latest.sif"
    publishDir "${output_dir}/cnvs", mode: "copy"

    input:
    tuple val(sampleId), path(txt)
    path pfb

    output:
    tuple val(sampleId), path("${sampleId}.cleaned.txt")


    script:
    """
    /home/user/PennCNV/clean_cnv.pl combineseg --signalfile ${pfb} ${txt} > ${sampleId}.cleaned.txt
    """
}

process penncnv2bed {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/beds", mode: "copy"

    input:
    tuple val(sampleId), path(txt)
    path sex_file

    output:
    tuple val(sampleId), path ("${sampleId}.bed")

    script:
    """
    python3 ${projectDir}/penncnv2bed.py --input_file ${txt} --output_file ${sampleId}.bed --sex_file ${sex_file}
    """
}

process classification {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    container "/mnt/dragen/pipelines/canvas/classifycnv_1.0.sif"
    publishDir "${output_dir}/ClassifyCNV/", mode: "copy"

    input:
    tuple val(sampleId), path(bed)

    output:
    tuple val(sampleId), path("${sampleId}/Scoresheet.txt"), emit: scoresheet
    path "${sampleId}/Intermediate_files/*.bed"

    script:
    """
    python3 /ClassifyCNV/ClassifyCNV.py --infile ${bed} --GenomeBuild hg19 --precise --outdir "${sampleId}"
    """
}

process addiscn {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/cnvs/", mode: "copy"

    input:
    tuple val(sampleId), path(txt)
    path(band)

    output:
    tuple val(sampleId), path("${sampleId}_iscn.txt")

    script:
    """
    python3 ${projectDir}/getISCN.py ${txt} ${band} "${sampleId}_iscn.txt"

    """
}

process makeplots {
    tag "$sampleId"
    container "/mnt/dragen/pipelines/canvas/bedgraph-visualizer_latest.sif"
    publishDir "${output_dir}/plots", mode: "copy"

    input:
    tuple val(sampleId), path(BAF), path(LRR), path(bed)

    output:
    tuple val(sampleId), path("${sampleId}/")

    script:
    """
    sed 's/^chr//' ${bed} > regions
    Rscript /usr/src/app/bedgraph-visualizer.R region_plot ${BAF} ${LRR} regions ${sampleId}

    Rscript /usr/src/app/bedgraph-visualizer.R genome_plot ${BAF} ${LRR} ${sampleId}
    """
}

process maketemplate {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/texs", mode: "copy"

    input:
    tuple val(sampleId), path(scoresheet), path(cnv), path(plot_dir)
    each path(tex_template)

    output:
    tuple val(sampleId), path("${sampleId}.tex"), path(plot_dir)

    script:
    """
    python3 ${projectDir}/tex_template_compile.py --scoresheet_file ${scoresheet} --cnv_file ${cnv} --tex_template ${tex_template} --output_file "${sampleId}.tex" --plot_dir ${plot_dir}
    """
}


process makereport {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    container '/mnt/dragen/pipelines/canvas/texlive_latest.sif'
    publishDir "${output_dir}/pdfs", mode: "copy"

    input:
    tuple val(sampleId), path(tex), path(plot_dir)
    path(template_dir)

    output:
    tuple val(sampleId), path ("${sampleId}.pdf")

    script:
    """
    cp "${template_dir}"/* ./
    pdflatex ${tex}
    pdflatex ${tex}
    """
