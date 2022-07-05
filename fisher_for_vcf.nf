/*
#########################################################################
##                                                                     ##
##     fisher_for_vcf.nf                                               ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Computational Biology Department                                ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################
*/


//////// Arguments of nextflow run

params.modules = ""

//////// end Arguments of nextflow run


//////// Variables

// from the nextflow.config file
config_file = file("${projectDir}/fisher_for_vcf.config")
log_file = file("${launchDir}/.nextflow.log")
cute_file=file(cute_path) // converted to file directly to use it as a constant
// end from the nextflow.config file

// from parameters
modules = params.modules // remove the dot -> can be used in bash scripts
// end from parameters


//////// end Variables


//////// Variables from config.file that need to be modified

sample_path_test = file("${sample_path}") // to test if exist below
tbi_path_test = file("${sample_path}.tbi") // to test if exist below
ped_path_test = file("${ped_path}") // to test if exist below
chr_path_test = file("${chr_path}") // to test if exist below

//////// end Variables from config.file that need to be modified


//////// Channels

//// used once

vcf_ch = Channel.fromPath("${sample_path}", checkIfExists: false) // I could use true, but I prefer to perform the check below, in order to have a more explicit error message
tbi_ch = Channel.fromPath("${sample_path}.tbi", checkIfExists: false) // I could use true, but I prefer to perform the check below, in order to have a more explicit error message
ped_ch = Channel.fromPath("${ped_path}", checkIfExists: false) // I could use true, but I prefer to perform the check below, in order to have a more explicit error message
chr_ch = Channel.fromPath("${chr_path}", checkIfExists: false) // I could use true, but I prefer to perform the check below, in order to have a more explicit error message
if(region == 'None'){
    region_ch = Channel.from("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chrY", "chrX", "chrM") // .split(",") split according to comma and create a tuple
}else{
    //String[] tempo
    tempo = region.split(",") // .split(",") split according to comma and create an array https://www.tutorialspoint.com/groovy/groovy_split.htm
    region_ch = Channel.from(tempo) 
}

//// end used once

//////// end Channels



//////// Checks

if(system_exec == 'local' || system_exec == 'slurm'){
    def file_exists1 = sample_path_test.exists()
    if( ! file_exists1){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE: ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    def file_exists2 = tbi_path_test.exists()
    if( ! file_exists2){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID .tbi FILE ASSOCIATED TO sample_path PARAMETER IN nextflow.config FILE: ${sample_path}.tbi\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    def file_exists3 = ped_path_test.exists()
    if( ! file_exists3){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ped_path PARAMETER IN nextflow.config FILE: ${ped_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    def file_exists4 = ped_path_test.exists()
    if( ! file_exists4){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID chr_path PARAMETER IN nextflow.config FILE: ${chr_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
}else{
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID system_exec PARAMETER IN nextflow.config FILE: ${system_exec}\n\n========\n\n"
}

//////// end Checks



//////// Processes


process WorkflowVersion { // create a file with the workflow version in out_path
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy'
    cache 'false'

    output:
    file "Run_info.txt"

    script:
    """
    echo "Project (empty means no .git folder where the main.nf file is present): " \$(git -C ${projectDir} remote -v | head -n 1) > Run_info.txt # works only if the main script run is located in a directory that has a .git folder, i.e., that is connected to a distant repo
    echo "Git info (empty means no .git folder where the main.nf file is present): " \$(git -C ${projectDir} describe --abbrev=10 --dirty --always --tags) >> Run_info.txt # idem. Provide the small commit number of the script and nextflow.config used in the execution
    echo "Cmd line: ${workflow.commandLine}" >> Run_info.txt
    echo "execution mode": ${system_exec} >> Run_info.txt
    modules=$modules # this is just to deal with variable interpretation during the creation of the .command.sh file by nextflow. See also \$modules below
    if [[ ! -z \$modules ]] ; then
        echo "loaded modules (according to specification by the user thanks to the --modules argument of main.nf)": ${modules} >> Run_info.txt
    fi
    echo "Manifest's pipeline version: ${workflow.manifest.version}" >> Run_info.txt
    echo "result path: ${out_path}" >> Run_info.txt
    echo "nextflow version: ${nextflow.version}" >> Run_info.txt
    echo -e "\\n\\nIMPLICIT VARIABLES:\\n\\nlaunchDir (directory where the workflow is run): ${launchDir}\\nprojectDir (directory where the main.nf script is located): ${projectDir}\\nworkDir (directory where tasks temporary files are created): ${workDir}" >> Run_info.txt
    echo -e "\\n\\nUSER VARIABLES:\\n\\nout_path: ${out_path}\\nsample_path: ${sample_path}" >> Run_info.txt
    """
}
//${projectDir} nextflow variable
//${workflow.commandLine} nextflow variable
//${workflow.manifest.version} nextflow variable
//Note that variables like ${out_path} are interpreted in the script block



process fisher {
    label 'python' // see the withLabel: bash in the nextflow config file 
    //publishDir path: "${out_path}", mode: 'copy', overwrite: false
    cache 'true'

    //no channel input here for the vcf, because I do not transform it
    input:
    tuple val(region2), file(vcf) from region_ch.combine(vcf_ch) // parallelization expected for each value of region_ch
    //file vcf from vcf_ch
    file ped from ped_ch.first()
    file tbi from tbi_ch.first()
    //val region2 from region_ch

    output:
    file "*.tsv" into fisher_ch1 // multi channel

    script:
    """
    #!/bin/bash -ue
    fisher_lod.py ${vcf} ${ped} "${region2}"
    """
}

fisher_ch1.collectFile(name: "fisher.tsv", skip:1, keepHeader:true).into{fisher_ch2 ; fisher_ch3}
fisher_ch2.subscribe{it -> it.copyTo("${out_path}")}


process miamiplot {
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}", mode: 'copy', pattern: "{*.png}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{miami_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    //no channel input here for the vcf, because I do not transform it
    input:
    file fisher from fisher_ch3
    file chr from chr_ch
    val x_lim
    val y_lim1
    val y_lim2
    file cute_file

    output:
    file "*.png"
    file "miami_report.txt"

    script:
    """
    #!/bin/bash -ue
    miami.R ${fisher} ${chr} "${x_lim}" "${y_lim1}" "${y_lim2}" "${cute_file}" "miami_report.txt"
    """

}



process Backup {
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', overwrite: false // since I am in mode copy, all the output files will be copied into the publishDir. See \\wsl$\Ubuntu-20.04\home\gael\work\aa\a0e9a739acae026fb205bc3fc21f9b
    cache 'false'

    input:
    file config_file
    file log_file

    output:
    file "${config_file}" // warning message if we use file config_file
    file "${log_file}" // warning message if we use file log_file
    file "Log_info.txt"

    script:
    """
    echo -e "full .nextflow.log is in: ${launchDir}\nThe one in the result folder is not complete (miss the end)" > Log_info.txt
    """
}


//////// end Processes
