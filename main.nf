#!/usr/bin/env nextflow


/*
========================================================================================
                         nf-gigaassay
========================================================================================
 pipeline to process & identify differentially abundant transcripts in SARS CoV-2 IsoSeq3 data
 #### Homepage / Documentation
TODO
----------------------------------------------------------------------------------------
*/

//TODO package each of the scripts into a script in bind and execute with find xargs 
// look at top answer below for guidance; find each file and execute seperately
//https://stackoverflow.com/questions/40700230/find-xargs-execute-chain-of-commands-for-each-file

params.samplesheet = "/wynton/group/krogan/mgordon/projects/100423_VRamani_isoseq3/isoseq3.SARSCoV2/docs/new.samplesheet.csv" 
params.subseq = 'CTTTCGATCTCTTGTAGATCTGTTCTC' // SARS CoV2 leader sequence; using this to extract reads from fastq files
params.scripts = "$projectDir/bin" // scripts for pipeline here
params.host_reference = '/wynton/group/krogan/apelin/genomes/Homo_sapiens/ensembl_release-101/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz' // Host genome
 
def helpMessage() {

    log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run nf-gigaassay  -profile singularity --reads 'path/to/*.fastq.gz' --reference 'path/to/refgenome.fa' --subseq 'CTTTCGATCTCTTGTAGATCTGTTCTC'
    Mandatory arguments:
      --reads [str]                 Full path to directory containing the input reads for analysis
      -profile [str]                Configuration profile to use. Currently supports conda, docker & singularity 

    References:
      --reference [file]            Fasta reference genome file used to perform read alignment       
      --subseq [file]               Substring used for read extraction from file
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

log.info """\
    ISOSEQ3 VARIANT PROCESSING P I P E L I N E
    ===================================
    leaderSequence : ${params.subseq}
    """
    .stripIndent(true)

/*
 * Indexing host genome
 */

process PBMM2_INDEX_HOST {
    tag "Indexing host reference"
    label 'process_high'
    publishDir "${params.outdir}/alignment/pbmm2/index", mode: 'copy'

    conda "bioconda::pbmm2==1.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmm2:1.9.0--h9ee0642_0':
        'biocontainers/pbmm2:1.9.0--h9ee0642_0' }"

    input:
    path fasta

    output:
    path("*.mmi"), emit: host_indx

    script:
    """
    hname=\$(echo $fasta | sed 's/.fa.gz/.mmi/g')
    fname=\$(echo $fasta | sed 's/.fa.gz/.fasta/g')

    gunzip -c  $fasta > \$fname

    pbmm2 index \\
        \$fname  \\
        \$hname \\
        --preset "ISOSEQ"
    """
}

/*
 * Convert unaligned bam to fastq
 */

process BAMTOOLS_CONVERT {
    tag "Converting $sample_id unaligned bam to fastq"
    label 'process_medium'
    publishDir "${params.outdir}/preproc/fastq", mode: 'copy'

    conda "bioconda::bamtools=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamtools:2.5.2--hdcf5f25_2' :
        'biocontainers/bamtools:2.5.2--hdcf5f25_2' }"

    input:
    tuple val(sample_id), path(bam), path(pbi), path(genome), path(fasta) // need to specify all?

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: fastq

    script:
    """
    bamtools convert \\
        -format fastq \\
        -in $bam \\
        -out ${sample_id}.fastq

    gzip ${sample_id}.fastq
    """
}

/*
* QC Assessment of the LR fastq
*/

process NANOPLOT {
    tag "QC of $sample_id fastq"
    label 'process_medium'
    publishDir "${params.outdir}/QC/nanoplot", mode: 'copy'

    conda "bioconda::nanoplot=1.41.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.41.6--pyhdfd78af_0' :
        'biocontainers/nanoplot:1.41.6--pyhdfd78af_0' }"

    input:
    tuple val(sample_id), path(bam), path(pbi), path(genome), path(fasta) // need to specify all?

    output:
    tuple val(sample_id), path("*.html")                , emit: html
    tuple val(sample_id), path("*.png") , optional: true, emit: png
    tuple val(sample_id), path("*.txt")                 , emit: txt
    tuple val(sample_id), path("*.log")                 , emit: log

    script:
    """
    NanoPlot \\
        -t $task.cpus \\
        --prefix $sample_id \\
        --fasta $fasta
    """
}

/*
 * Identify & extract reads with the SARS-CoV2 leader sequence
 */

process BBTOOLS_BBDUK_FILTER {
    tag "Extracting leader sequences from fasta.."
    label 'process_medium'
    publishDir "${params.outdir}/preproc/extract", mode: 'copy'

    conda "bioconda::bbmap=39.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    val subseq
    tuple val(sample_id), path(bam), path(pbi), path(genome), path(fasta)

    output:
    tuple val(sample_id), path("*.match.fasta"), path(genome), emit: match_fa
    tuple val(sample_id), path("*.nomatch.fasta"), emit: nomatch_fa
    path("*.txt")

    script:
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g') 

    bbduk.sh \\
        -Xmx\$maxmem \\
        in=$fasta \\
        literal=$subseq \\
        K=27 \\
        mm=f \\
        out=${sample_id}.nomatch.fasta \\
        outm=${sample_id}.match.fasta \\
        lhist=${sample_id}.lhist.txt \\
        stats=${sample_id}.filtering.txt
    """
}

/*
 * trim the leader sequence from the reads to keep just the transcripts
 */

// process BBTOOLS_BBDUK_TRIM {
//     tag "Trimming SARS-CoV2 leader sequences"
//     label 'process_medium'
//     publishDir "${params.outdir}/preproc/trimmed", mode: 'copy'

//     conda "bioconda::bbmap=39.01"
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
//         'biocontainers/bbmap:39.01--h5c4e2a8_0' }"

//     input:
//     tuple val(sample_id), path(fasta), path(genome) // need to specify all?


//     output:
//     tuple val(sample_id), path("*.trim.match.fasta.gz"), path(genome) emit: trimmed_fa
//     path("*.txt")

//     script:
//     """
//     maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g') 

//     bbduk.sh \\
//         -Xmx\$maxmem \\
//         in=$fasta \\
//         literal=$subseq \\
//         ktrim=l \\ 
//         out=${sample_id}.trim.match.fasta.gz \\
//         stats=${sample_id}.trimming.txt
//     """

// }

/*
 * Align the filtered reads using minimap2 
 * use ISOSEQ3 preset adjust for intron length size (just make genome size to capture everything)
 * dropping cost for non-canonical splicing 
 */

// process PBMM2_ALIGN_VIRUS {
//     tag "Aligning $sample_id sequences to viral reference"
//     label 'process_high'
//     publishDir "${params.outdir}/alignment/pbmm2/variants", mode: 'copy'

//     conda "bioconda::pbmm2==1.9.0"
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/pbmm2:1.9.0--h9ee0642_0':
//         'biocontainers/pbmm2:1.9.0--h9ee0642_0' }"

//     input:
//     tuple val(sample_id), path(fasta), path(genome)

//     output:
//     tuple val(sample_id), path("*.align.bam"),  emit: variant_bam
//     path("*.alignment.log")

//     script:
//     """
//     variant=\$(echo $genome | sed 's/UK5-Calu3-//g;s/.fa//g') 

//     echo \$variant

//     pbmm2 align \\
//         $genome \\
//         $fasta \\
//         ${sample_id}.align.bam \\
//         --sort \\
//         --preset "ISOSEQ" \\
//         -g 30000 -G 30000 \\ 
//         --rg '@RG\tID:\$variant\tSM:${sample_id}' \\ 
//         --log-level INFO > ${sample_id}.variant.alignment.log 
//     """
// }

/*
 * Align the filtered reads using minimap2 to host 
 * use ISOSEQ3 preset
 */

// process PBMM2_ALIGN_HOST {
//     tag "Aligning $sample_id sequences to viral reference"
//     label 'process_high'
//     publishDir "${params.outdir}/alignment/pbmm2/host", mode: 'copy'

//     conda "bioconda::pbmm2==1.9.0"
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/pbmm2:1.9.0--h9ee0642_0':
//         'biocontainers/pbmm2:1.9.0--h9ee0642_0' }"

//     input:
//     path index
//     tuple val(sample_id), path(fasta), path(genome)

//     output:
//     tuple val(sample_id), path("*.align.bam"),  emit: host_bam
//     path("*.alignment.log")

//     script:
//     """
//     hname=\$(basename $genome| sed 's/.mmi//g')

//     pbmm2 align \\
//         $genome \\
//         $fasta \\
//         ${sample_id}.host.align.bam \\
//         --sort \\
//         --preset "ISOSEQ" \\
//         --rg '@RG\tID:\$hname\tSM:${sample_id}' \\ 
//         --log-level INFO > ${sample_id}.host.alignment.log 
//     """
// }


workflow {
    // parse csv for sample info 
    // must match names and order in samplesheet

    input_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header:true)
        .map { row-> tuple(row.Sample_Name, file(row.Sample_BAM), file(row.Sample_PBI), file(row.Genome), file(row.Sample_FASTA)) } 

    // leader sequence channel
    leader_ch = Channel
        .value(params.subseq)

    // host reference
    href_ch = Channel
        .fromPath(params.host_reference, checkIfExists: true)

    // index the host genome  for faster mapping
    hindex_ch    = PBMM2_INDEX_HOST(href_ch)

    // convert to fastq file
    fq_ch         = BAMTOOLS_CONVERT(input_ch)

    // qc plots of the samples
    qc_plots_ch   = NANOPLOT(input_ch)

    // extracting reads with leader sequences from fastq
    filtered_fq_ch   = BBTOOLS_BBDUK_FILTER(leader_ch, input_ch)

    // trim the 5' end to remove leader sequence prior to alignment
    //trimmed_fq_ch   = BBTOOLS_BBDUK_TRIM(leader_ch, filtered_fq_ch.match_fa)

    // map viral reads to the variant genome
    //viral_bam_ch = PBMM2_ALIGN_VIRUS(trimmed_fq_ch.trimmed_fa)
   // map host reads to the host genome
    //host_bam_ch  = PBMM2_ALIGN_HOST(hindex_ch.host_indx, filtered_fq_ch.nomatch_fa)


}

