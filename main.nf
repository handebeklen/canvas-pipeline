#!/usr/bin/env nextflow

nextflow.enable.dsl2

params.bpm = 'manifest-cluster/GSACyto_20044998_A1.bpm'
params.csv = 'gtcler/GSA-24v1-0_C1.csv'
params.egt = 'manifest-cluster/2003.egt'
params.genome_fasta = '/staging/references/hg19/hg19.fa'
params.sexfile = './sexfile'
params.penncnv_pfb = 'test/out.pfb'
params.output_dir = 'outputs'
params.samplesheet = 'samplesheets/chip001.tsv'  // Specify the default or input samplesheet
params.idat_folder = 'idats/chip001' // Specify the IDAT folder

workflow {
    idat_folder = Channel.fromPath(params.idat_folder)
    
    idat_folder
        .map { dir ->
            def chip = dir.name.split('/')[1]
            def gtcs = "${params.output_dir}/chip${chip}/gtcs"
            idat2gtc(dir, gtcs)
        }
        .set { gtcs }

    gtcs
        .map { gtc_folder ->
            def vcfs = "${params.output_dir}/chip001/vcfs"
            gtc2vcf(gtc_folder, vcfs)
        }
        .set { vcfs }

    vcfs
        .map { vcf_folder ->
            def bedgraphs = "${params.output_dir}/chip001/bedgraphs"
            vcf2bedgraph(vcf_folder, bedgraphs)
        }
        .set { bedgraphs }

    bedgraphs
        .map { bedgraph_folder ->
            def cnv_input = "${params.output_dir}/chip001/beds"
            vcf2penncnv(bedgraph_folder, cnv_input)
        }
        .set { cnv_input }

    cnv_input
        .map { cnv_file ->
            def combined_cnv = "${params.output_dir}/chip001/qcs"
            penncnv_detect_all(cnv_file, combined_cnv)
        }
        .set { combined_cnv_outputs }

    combined_cnv_outputs
        .map { combined_cnv_file ->
            def clean_cnv = "${params.output_dir}/chip001/beds"
            penncnv_clean_cnv(combined_cnv_file, clean_cnv)
        }
        .set { clean_cnv_outputs }

    clean_cnv_outputs
        .map { clean_cnv_file ->
            def raw_cnv_bed = "${params.output_dir}/chip001/beds"
            outpenncnv2bed(clean_cnv_file, raw_cnv_bed)
        }
        .set { raw_cnv_bed_outputs }

    raw_cnv_bed_outputs
        .map { raw_cnv_bed_file ->
            def pdfs = "${params.output_dir}/chip001/pdfs"
            classification(raw_cnv_bed_file, pdfs)
        }
}

process idat2gtc {
    input:
    path idat_folder
    path gtcs

    output:
    path gtcs

    script:
    """
    /opt/dragena-linux-x64-DAv1.0.0/dragena/dragena genotype call --bpm-manifest ${params.bpm} --cluster-file ${params.egt} --idat-folder ${idat_folder} --output-folder ${gtcs}
    """
}

process gtc2vcf {
    input:
    path gtc_folder
    path vcfs

    output:
    path vcfs

    script:
    """
    /opt/dragena-linux-x64-DAv1.0.0/dragena/dragena genotype gtc-to-vcf --bpm-manifest ${params.bpm} --csv-manifest ${params.csv} --genome-fasta-file ${params.genome_fasta} --gtc-folder ${gtc_folder} --output-folder ${vcfs}
    """
}

process vcf2bedgraph {
    input:
    path vcf_folder
    path bedgraphs

    output:
    path bedgraphs

    script:
    """
    bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%REF\t%ALT\n' ${vcf_folder} > ${bedgraphs}
    """
}

process vcf2penncnv {
    input:
    path vcf_folder
    path cnv_input

    output:
    path cnv_input

    script:
    """
    printf "Name\tChr\tPosition\tGType\tLog R Ratio\tB Allele Freq" > ${cnv_input}
    bcftools query -f '%ID\t%CHROM\t%POS[\t%GT\t%LRR\t%BAF]\n' ${vcf_folder} | sed 's/1\\/1/BB/g;s/0\\/0/AA/g;s/0\\/1/AB/g;s/.\\/./NC/g' >> ${cnv_input}
    """
}

process penncnv_detect_all {
    input:
    path cnv_input
    path combined_cnv

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
    input:
    path combined_cnv_file
    path clean_cnv

    output:
    path clean_cnv

    script:
    """
    ./penncnv_latest.sif ./clean_cnv.pl combineseg --signalfile ${params.penncnv_pfb} > ${clean_cnv}
    """
}

process outpenncnv2bed {
    input:
    path clean_cnv_file
    path raw_cnv_bed

    output:
    path raw_cnv_bed

    script:
    """
    python3 penncnv2bed.py ${clean_cnv_file} ${raw_cnv_bed}
    """
}

process classification {
    input:
    path raw_cnv_bed
    path pdfs

    output:
    path pdfs

    script:
    """
    ./classifycnv_1.0.sif python3 /ClassifyCNV/ClassifyCNV.py --infile ${raw_cnv_bed} --GenomeBuild hg19 --precise --outfile ${pdfs}
    """
}

