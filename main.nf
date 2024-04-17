#!/usr/bin/env nextflow
/*
=============================================================================
|                        dsl2-template-nf                                   |
=============================================================================
|    #### Homepage / Documentation                                          |
|    https://github.com/lifebit-ai/template4users/blob/main/docs/README.md |
-----------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/* --------------------
| Summary              |
--------------------- */

def summary = [:]

if (workflow.revision) summary['Pipeline Release'] = workflow.revision

summary['Launch dir']                                  = workflow.launchDir
summary['Working dir']                                 = workflow.workDir
summary['Script dir']                                  = workflow.projectDir
summary['User']                                        = workflow.userName
summary['Output dir']                                  = params.outdir

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

/* --------------------
| Help Message         |
--------------------- */

def helpMessage() {
    if ( workflow.userName != "ec2-user" ) {
        log.info comptext()
    }
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:

    Mandatory:
    --input                Input file

    Resource Options:
    --max_cpus            Maximum number of CPUs (int)
                        (default: $params.max_cpus)  
    --max_memory          Maximum memory (memory unit)
                        (default: $params.max_memory)
    --max_time            Maximum time (time unit)
                        (default: $params.max_time)

    """.stripIndent()
}

// Print completion text to stdout
if ( workflow.userName != "ec2-user" ) {
    print comptext()
}

/* --------------------
| Process definitions  |
--------------------- */

process OBTAIN_PIPELINE_METADATA {
    label "process_micro"
    publishDir "${params.tracedir}", mode: "copy"

    input:
    val repository
    val commit
    val revision
    val script_name
    val script_path
    val project_dir
    val launch_dir
    val work_dir
    val user_name
    val command_line
    val config_paths
    val propath
    val container
    val container_engine
    val raci_owner
    val domain_keywords

    output:
    path("pipeline_metadata_report.tsv"), emit: ch_pipeline_metadata_report

    shell:
    '''
    echo "Repository\t!{repository}"                  > temp_report.tsv
    echo "Commit\t!{commit}"                         >> temp_report.tsv
    echo "Revision\t!{revision}"                     >> temp_report.tsv
    echo "Script name\t!{script_name}"               >> temp_report.tsv
    echo "Script path\t!{script_path}"               >> temp_report.tsv
    echo "Project directory\t!{project_dir}"         >> temp_report.tsv
    echo "Launch directory\t!{launch_dir}"           >> temp_report.tsv
    echo "Work directory\t!{work_dir}"               >> temp_report.tsv
    echo "User name\t!{user_name}"                   >> temp_report.tsv
    echo "Command line\t!{command_line}"             >> temp_report.tsv
    echo "Configuration path(s)\t!{config_paths}"    >> temp_report.tsv
    echo "Propath\t!{propath}"                       >> temp_report.tsv
    echo "Container\t!{container}"                   >> temp_report.tsv
    echo "Container engine\t!{container_engine}"     >> temp_report.tsv
    echo "RACI owner\t!{raci_owner}"                 >> temp_report.tsv
    echo "Domain keywords\t!{domain_keywords}"       >> temp_report.tsv

    awk 'BEGIN{print "Metadata_variable\tValue"}{print}' OFS="\t" temp_report.tsv > pipeline_metadata_report.tsv
    '''
}

// EDIT HERE
process process1 {
    label "process_large"
    label "process1"
    echo true
    publishDir "${params.outdir}", mode: "copy"

    input:
    // Define your input files e.g. path(input_paths_file)


    output:
    // Define your Output files and variables e.g. tuple val(params.add1), path("*.txt")
   

    script:
    """
    # Write your script of commands here
    echo "[INFO] Running the process!"
    // ------------------------------------------
    // e.g. bcftools view \
    //    -Oz --output ${output.vcf.gz} \
    //    --regions chr1 \
    //    --threads ${params.large_cpus} \
    //    $vcf_paths_file
    echo "[INFO] Completed process1"
    // ------------------------------------------
    """
}
// TOO HERE

// EDIT workflow name
workflow process1 {

    main:

  // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }

    /*---------------------------------------------------
    | Setting up introspection variables and channels    |
    ----------------------------------------------------*/
    ch_repository = Channel.of(workflow.manifest.homePage)
    ch_commit_id = Channel.of(workflow.commitId ?: "Not available is this execution mode. Please run 'nextflow run ${workflow.manifest.homePage} [...]' instead of 'nextflow run main.nf [...]'")
    ch_revision = Channel.of(workflow.manifest.version)
    ch_script_name = Channel.of(workflow.scriptName)
    ch_script_file = Channel.of(workflow.scriptFile)
    ch_project_dir = Channel.of(workflow.projectDir)
    ch_launch_dir = Channel.of(workflow.launchDir)
    ch_work_dir = Channel.of(workflow.workDir)
    ch_user_name = Channel.of(workflow.userName)
    ch_command_line = Channel.of(workflow.commandLine)
    ch_config_files = Channel.of(workflow.configFiles)
    ch_profile = Channel.of(workflow.profile)
    ch_container = Channel.of(workflow.container)
    ch_container_engine = Channel.of(workflow.containerEngine)

    /*------------------------------------------------------------------
    | Setting up additional variables used for documentation purposes   |
    -------------------------------------------------------------------*/
    ch_raci_owner = Channel.of(params.raci_owner) 
    ch_domain_keywords = Channel.of(params.domain_keywords)

    /*------------------------
    | Setting up input data   |
    -------------------------*/
    // Define channels from repository files
    project_dir = workflow.projectDir

    // Define Channels from input
    // Edit channels, for multiple conditional channels use if and else statements
    // EDIT HERE
    ch_input_file = Channel.fromPath("${params.input}")
    // TOO HERE
    

    /*--------------
    | Processes    |
    --------------*/
    // Do not delete this process
    // Create introspection report
    if (params.obtain_pipeline_metadata) {
        OBTAIN_PIPELINE_METADATA(
            ch_repository,
            ch_commit_id,
            ch_revision,
            ch_script_name,
            ch_script_file,
            ch_project_dir,
            ch_launch_dir,
            ch_work_dir,
            ch_user_name,
            ch_command_line,
            ch_config_files,
            ch_profile,
            ch_container,
            ch_container_engine,
            ch_raci_owner,
            ch_domain_keywords
        )
    }
    
// Run process1 with our input channel parameter
 // EDIT HERE
        process1(ch_input_file)

}

workflow {

   process1()

}
 // TOO HERE

// Trace report
user_name = workflow.userName

if (user_name == "ubuntu" || user_name == "ec2-user") {
    workflow.onComplete {
        def trace_timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
        trace_report = file("/home/${user_name}/nf-out/trace.txt")
        trace_report.copyTo("results/pipeline_info/execution_trace_${trace_timestamp}.txt")
    }
}

// ANSII string of completion text
def comptext() {
    text  = """
    .................................................................................................
    .................................................................................................
    .......CONGRATULATIONS FOR WRITING & RUNNING YOUR FIRST PIPELINE ON THE LIFEBIT PLATFORM!........
    .................................................................................................
    .................................................................................................
    """.stripIndent()
    return text
}
