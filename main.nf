#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/toolprofiler
========================================================================================
 nf-core/toolprofiler Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/toolprofiler
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
### -> for the helper / checks / header code:
### -> this is all thanks to @drpatelh and the work done
### -> for the RNAseq pipeline  
*/


////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run lescailab/toolprofiler --input samplesheet.csv -profile test,docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)


////////////////////////////////////////////////////
/*          INCLUDE MODULES FOR TESTING           */
////////////////////////////////////////////////////

def bwa_index_options = [:]
def bwa_mem_options = [:]
def markduplicates_options = [:]
def salmon_index_options = [:]
def salmon_quant_options = [:]

include { BWA_INDEX } from './modules/nf-core-mod/software/bwa/index' addParams(bwa_index_options)
include { BWA_MEM } from './modules/nf-core-mod/software/bwa/mem' addParams(bwa_mem_options)
include { PICARD_MARKDUPLICATES } from './modules/nf-core-mod/software/picard/markduplicates' addParams(markduplicates_options)
include { SALMON_INDEX } from './modules/nf-core-mod/software/salmon/index' addParams(salmon_index_options)
include { SALMON_QUANT } from './modules/nf-core-mod/software/salmon/quant' addParams(salmon_quant_options)


workflow DNA {
    take:
        inputSample

    if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Genome fasta file not specified!' }
    ch_bwa_index = Channel.empty()
    ch_bwa_index = BWA_INDEX(ch_fasta).index
    BWA_MEM(inputSample, ch_bwa_index, ch_fasta)
    PICARD_MARKDUPLICATES(BWA_MEM.out.bam)

}


workflow RNA {
    take:
        inputSample
    // Check mandatory parameters
    if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Genome fasta file not specified!' }
    if (params.transcript_fasta) { ch_transcript_fasta = file(params.transcript_fasta) } else { exit 1, 'Transcript fasta file not specified!' }
    if (!params.gtf) { exit 1, "No GTF annotation specified!" }
    
    ch_salmon_index   = Channel.empty()
    ch_salmon_index   = SALMON_INDEX( ch_fasta, ch_transcript_fasta ).index
    ch_dummy = Channel.empty()
    SALMON_QUANT( reads, ch_salmon_index, gtf, ch_dummy, false)
    
}


workflow {

    /*
    ################ INPUT FILE PARSING ####################
    */
    input = file(params.input)
    inputSample = Channel.empty()
    inputSample = readInputFile(input, params.single_end)
    reference = file(params.reference)

    if(params.mode == 'dna'){
        DNA(inputSample)
    }
    else {
        RNA(inputSample)
    }



}





// ############## WARNING !!! #######################################
// the part below is going to be transferred to a module / class soon
// ############## UTILITIES AND SAMPLE LOADING ######################

// ### preliminary check functions

def checkExtension(file, extension) {
    file.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def checkFile(filePath, extension) {
  // first let's check if has the correct extension
  if (!checkExtension(filePath, extension)) exit 1, "File: ${filePath} has the wrong extension. See --help for more information"
  // then we check if the file exists
  if (!file(filePath).exists()) exit 1, "Missing file in TSV file: ${filePath}, see --help for more information"
  // if none of the above has thrown an error, return the file
  return(file(filePath))
}

// the function expects a tab delimited sample sheet, with a header in the first line
// the header will name the variables and therefore there are a few mandatory names
// sampleID to indicate the sample name or unique identifier
// read1 to indicate read_1.fastq.gz, i.e. R1 read or forward read
// read2 to indicate read_2.fastq.gz, i.e. R2 read or reverse read
// any other column should fulfill the requirements of modules imported in main
// the function also expects a boolean for single or paired end reads from params

def readInputFile(tsvFile, single_end) {
    Channel.from(tsvFile)
        .splitCsv(header:true, sep: '\t')
        .map { row ->
            def meta = [:]
            def reads = []
            def sampleinfo = []
            meta.sampleID = row.sampleID
            if (single_end) {
              reads = checkFile(row.read1, "fastq.gz")
            } else {
              reads = [ checkFile(row.read1, "fastq.gz"), checkFile(row.read2, "fastq.gz") ]
            }
            sampleinfo = [ meta, reads ]
            return sampleinfo
        }
}