#!/usr/bin/env nextflow 
nextflow.enable.dsl=2

//INPUTs
params.infastq = "/scicomp/home-pure/tsz0/Projects/Homopolymer_detection_work/mitsko_sterne_day0.fastq"
//params.infiles = params.indir + "*.fastq"
//params.outdir = "$HOME/Projects/Homopolymer_detection_work/homopolymer_analysis"
params.outfile = "$HOME/Projects/Homopolymer_detection_work/sterne_day0_homopolymer_counts.txt"
params.reference = "$HOME/Projects/Homopolymer_detection_work/Sterne_pilon.fasta"
params.homopolymers = "$HOME/BDRD_group_drive/Passage_Strain_Analysis/Variant_Analysis/Sterne/sterne_5n.txt"

process HOPE {
    cpus 1
	tag{"hope ${fastq}"}
	//publishDir( "${params.outdir}", mode: 'copy')
	
	input:
	val fastq

	output:
	//tuple path("${fastq.baseName}.bam"), path("${fastq.baseName}_homopolymers_out.txt"), emit: hope_out
    path("${fastq.baseName}_homopolymers_out.txt"), emit: hope_out
	script:
	"""
	minimap2 -t ${task.cpus} -A 2 -B 10 -a ${params.reference} ${fastq} | samtools sort -o ${fastq.baseName}.bam
	samtools index ${fastq.baseName}.bam
    ~/programs/hope/hope.py -t ${task.cpus} -a ${params.reference} -b ${fastq.baseName}.bam -f ${params.homopolymers} -o ${fastq.baseName}_homopolymers_
	"""	
}


workflow {
    //reads_files = Channel.fromPath( "/scicomp/home-pure/tsz0/Projects/Homopolymer_detection_work/test_stern0/*.fastq" )
	//read in files & chunk into 10,000 sequences per entry
    read_files = Channel
                        .fromPath(params.infastq)
                        .splitFastq( by: 10000, file: true)
                        .view()
    //run hope.py, concatenate the results
	HOPE(read_files).collectFile(
        name: params.outfile,
        keepHeader: true,
    )
}

workflow.onComplete {
   println ( workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}
