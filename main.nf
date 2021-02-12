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

include { BWA_INDEX } from './modules/nf-core-mod/software/bwa/index' params(params)
include { BWA_MEM } from './modules/nf-core-mod/software/bwa/mem' params(params)
include { PICARD_MARKDUPLICATES } from './modules/nf-core-mod/software/picard/markduplicates' params(params)
include { SALMON_INDEX } from './modules/nf-core-mod/software/salmon/index' params(params)
include { SALMON_QUANT } from './modules/nf-core-mod/software/salmon/quant' params(params)


workflow DNA {

}


workflow RNA {

    
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
        DNA(inputSample, reference)
    }
    else {
        RNA(inputSample, reference)
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