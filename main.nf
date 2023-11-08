nextflow.enable.dsl=1
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


//////// Options of nextflow run

params.modules = ""

//////// end Options of nextflow run


//////// Variables

// from the nextflow.config file
config_file = file("${projectDir}/fisher_for_vcf.config") // file() create a path object necessary o then create the file
log_file = file("${launchDir}/.nextflow.log")

// files objects created in order to use .exists() to test the path
chr = file(chr_path)
ped = file(ped_path)
cute = file(cute_path) // converted to file path directly to use it as a constant
out = file(out_path)
// end from the nextflow.config file

// from parameters (options of the nexflow command line)
modules = params.modules // remove the dot -> can be used in bash scripts
// end from parameters (options of the nexflow command line)


//////// end Variables


//////// Variables from config.file that need to be modified


if(x_lim == 'whole' || (x_lim == 'region' && region == 'none')){ // for the miami plot
    x_lim_val = "chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chrY, chrX, chrM" // I have added both notations "chr23", "chr24", "chr25" or "chrY", "chrX", "chrM" because either one or the other can be used in a VCF file 
}else if(x_lim == 'region'){
    x_lim_val = region
}else{
    x_lim_val = x_lim // value for the miami plot
}

//////// end Variables from config.file that need to be modified


//////// Channels


Channel.fromPath("${sample_path}", checkIfExists: false).into{vcf_ch1 ; vcf_ch2 ; vcf_ch3} // I could use true, but I prefer to perform the check below, in order to have a more explicit error message
if(region == 'none'){  // for combine below for parallelization of the fisher process
    region_ch = Channel.from("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chrY", "chrX", "chrM") // I have added both notations "chr23", "chr24", "chr25" or "chrY", "chrX", "chrM" because either one or the other can be used in a VCF file 
}else{
    if(region =~ /,/){
        tempo = region.replaceAll(':.+,', ',')
    }else{
        tempo = region
    }
    tempo2 = tempo.replaceAll(':.+$', '')
    tempo3 = tempo2.replaceAll(' ', '')
    tempo4 = tempo3.split(",") // .split(",") split according to comma and create an array https://www.tutorialspoint.com/groovy/groovy_split.htm
    region_ch = Channel.from(tempo4) 
}


//////// end Channels



//////// Checks

sample_path_test = file("${sample_path}") // because is a channel
tbi_test = file("${sample_path}.tbi")

def file_exists1 = sample_path_test.exists()
if( ! file_exists1){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE: ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
}else if(sample_path =~ /.*\.gz$/){
    def file_exists2 = tbi_test.exists()
    if( ! file_exists2){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID .tbi FILE ASSOCIATED TO sample_path PARAMETER IN nextflow.config FILE: ${sample_path}.tbi\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\nOTHERWISE, USE tabix -p vcf <NAME>.vcf TO INDEX THE .gz FILE\n\n========\n\n"
    }else{
        tbi = file("${sample_path}.tbi")
    }
}
def file_exists3 = ped.exists()
if( ! file_exists3){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ped_path PARAMETER IN nextflow.config FILE: ${ped_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
}
def file_exists4 = chr.exists()
if( ! file_exists4){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID chr_path PARAMETER IN nextflow.config FILE: ${chr_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
}

if( ! region in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID region PARAMETER IN nextflow.config FILE:\n${region}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! tsv_extra_fields in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tsv_extra_fields PARAMETER IN nextflow.config FILE:\n${tsv_extra_fields}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! x_lim in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID x_lim PARAMETER IN nextflow.config FILE:\n${x_lim}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! vgrid in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID vgrid PARAMETER IN nextflow.config FILE:\n${vgrid}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! top_y_column in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID top_y_column PARAMETER IN nextflow.config FILE:\n${top_y_column}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! bottom_y_column in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID bottom_y_column PARAMETER IN nextflow.config FILE:\n${bottom_y_column}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! color_column in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID color_column PARAMETER IN nextflow.config FILE:\n${color_column}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! y_lim1 in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_lim1 PARAMETER IN nextflow.config FILE:\n${y_lim1}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! y_lim2 in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_lim2 PARAMETER IN nextflow.config FILE:\n${y_lim2}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! y_reverse1 in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_reverse1 PARAMETER IN nextflow.config FILE:\n${y_reverse1}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! y_reverse2 in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_reverse2 PARAMETER IN nextflow.config FILE:\n${y_reverse2}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! y_threshold1 in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_threshold1 PARAMETER IN nextflow.config FILE:\n${y_threshold1}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! y_threshold2 in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_threshold2 PARAMETER IN nextflow.config FILE:\n${y_threshold2}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! y_log1 in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_log1 PARAMETER IN nextflow.config FILE:\n${y_log1}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! y_log2 in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_log2 PARAMETER IN nextflow.config FILE:\n${y_log2}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
def file_exists5 = cute.exists()
if( ! file_exists5){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN nextflow.config FILE:\n${cute_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
}

// below : those variable are already used in the config file. Thus, to late to check them. And not possible to check inside the config file
// system_exec
// out_ini
print("\n\nRESULT DIRECTORY: ${out_path}")
print("\n\nWARNING: PARAMETERS ALREADY INTERPRETED IN THE .config FILE:")
print("    system_exec: ${system_exec}")
print("    out_path: ${out_path_ini}")
print("    queue: ${queue}")
print("    qos: ${qos}")
print("    add_options: ${add_options}")
print("\n\n")


//////// end Checks



//////// Processes


process WorkflowVersion { // create a file with the workflow version in out_path
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir "${out}/reports", mode: 'copy'
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
        echo "loaded modules (according to specification by the user thanks to the --modules argument of main.nf): ${modules}" >> Run_info.txt
    fi
    echo "Manifest's pipeline version: ${workflow.manifest.version}" >> Run_info.txt
    echo "result path: ${out}" >> Run_info.txt
    echo "nextflow version: ${nextflow.version}" >> Run_info.txt
    echo -e "\\n\\nIMPLICIT VARIABLES:\\n\\nlaunchDir (directory where the workflow is run): ${launchDir}\\nprojectDir (directory where the main.nf script is located): ${projectDir}\\nworkDir (directory where tasks temporary files are created): ${workDir}" >> Run_info.txt
    echo -e "\\n\\nUSER VARIABLES:\\n\\nout_path: ${out}\\nsample_path: ${sample_path}" >> Run_info.txt
    """
}
//${projectDir} nextflow variable
//${workflow.commandLine} nextflow variable
//${workflow.manifest.version} nextflow variable
//Note that variables like ${out} are interpreted in the script block


process vcf_subfield_title {
    label 'r_ext' // see the withLabel: bash in the nextflow config file
    publishDir "${out}/reports", mode: 'copy', pattern: "{*.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    file vcf from vcf_ch1
    file cute

    output:
    file "vcf_info_field_titles.txt" into vcf_info_field_titles_ch
    file "vcf_csq_subfield_titles.txt" into vcf_csq_subfield_titles_ch
    file "vcf_subfield_title_report.txt"

    script:
    """
    #!/bin/bash -ue
    vcf_subfield_title.R ${vcf} "${cute}" "vcf_subfield_title_report.txt"
    """
}


process fisher {
    label 'python' // see the withLabel: bash in the nextflow config file 
    publishDir "${out}/reports", mode: 'copy', pattern: "{fisher_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    tuple val(region2), file(vcf) from region_ch.combine(vcf_ch2) // parallelization expected for each value of region_ch
    file ped
    if(sample_path =~ /.*\.gz$/){file tbi}
    file vcf_info_field_titles from vcf_info_field_titles_ch.first()
    file vcf_csq_subfield_titles from vcf_csq_subfield_titles_ch.first()
    val tsv_extra_fields

    output:
    file "*.tsv" into fisher_ch1 // multi channel
    file "*.txt"

    script:
    """
    #!/bin/bash -ue
    if grep -Eq "FISHER" ${vcf_info_field_titles} ; then # does not work without [[]]
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nTHE VCF FILE HAS ALREADY FISHER COMPUTATION PERFORMED, AS FIELDS ARE:\\n\$(cat ${vcf_info_field_titles})\\n\\n========\\n\\n"
        exit 1
    else
        fisher_lod.py ${vcf} ${ped} "${region2}" ${vcf_info_field_titles} "${tsv_extra_fields}" ${vcf_csq_subfield_titles} "fisher_report.txt"
    fi
    """
}

fisher_ch1.collectFile(name: "fisher.tsv", skip:1, keepHeader:true).into{fisher_ch2 ; fisher_ch3 ; fisher_ch4 ; fisher_ch5}
//fisher_ch2.subscribe{it -> it.copyTo("${out}")} // will be published below, after zipping


process miami_plot {
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir "${out}", mode: 'copy', pattern: "{*.png}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out}/reports", mode: 'copy', pattern: "{miami_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    file fisher from fisher_ch3
    file chr
    val x_lim_val
    val vgrid
    val top_y_column
    val bottom_y_column
    val color_column
    val dot_border_color
    val y_lim1
    val y_lim2
    val y_reverse1
    val y_reverse2
    val y_threshold1
    val y_threshold2
    val y_log1
    val y_log2
    file cute

    output:
    file "*.png"
    file "miami_report.txt"

    script:
    """
    #!/bin/bash -ue
    miami.R ${fisher} ${chr} "${x_lim_val}" "${vgrid}" "${top_y_column}" "${bottom_y_column}" "${color_column}" "${dot_border_color}" "${y_lim1}" "${y_lim2}" "${y_reverse1}" "${y_reverse2}" "${y_threshold1}" "${y_threshold2}" "${y_log1}" "${y_log2}" "${cute}" "miami_report.txt"
    """
}


process tsv2vcf {
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir "${out}", mode: 'copy', overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    file vcf from vcf_ch3
    file fisher from fisher_ch4

    output:
    file "res_fisher.*"

    script:
    """
    #!/bin/bash -ue
    PREHEADER='##fileformat=VCFv4.2;build by fisher_for_vcf.nf\\n##WARNING: This file is not a true VCF since FORMAT AND sample (indiv) columns are not present'
    HEADER='#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO'
    echo -e \$PREHEADER > res_fisher.vcf
    FILENAME=\$(basename -- "${vcf}") # recover a file name without path
    FILE_EXTENSION="\${FILENAME##*.}" #  ## means "delete the longest regex starting at the beginning of the tested string". If nothing, delete nothing. Thus ##*. means delete the longest string finishing by a dot. Use # instead of ## for "delete the shortest regex starting at the beginning of the tested string"
    if [[ "\${FILE_EXTENSION}" =~ gz ]] ; then
        zcat ${vcf} | awk '{
            if(\$0 ~ "^##.*"){
                print \$0
            }else{
                exit 0
            }
        }' >> res_fisher.vcf
    else
        awk '{
            if(\$0 ~ "^##.*"){
                print \$0
            }else{
                exit 0
            }
        }' ${vcf} >> res_fisher.vcf
    fi
    awk -v var1=\$HEADER 'BEGIN{FS="\\t" ; OFS="" ; ORS=""}
        NR==1{
            print "##WARNING: 5 first names of the header of the initial file: "\$1" "\$2" "\$3" "\$4" "\$5"\\n" ;
            print "##WARNING: if the 5 first columns of the .tsv file are not CHROM POS REF ALT INFO, then the .vcf file produced by this process is not good\\n" ;
            print "##INFO=<ID=FISHER,Number=.,Type=String,Description=\""Fisher exact tests based on the presence/absence of the variant in the affected/unaffected indiv. Format: " ;
            for(i=6;i<=NF;i++){print \$i ; if(i < NF){print "|"}} ;
            print "\"">\\n" ;
            print var1"\\n"
        }
        NR > 1{
            gsub("[\\\\[\\\\]\\'"'"']", "", \$4)
            print \$1"\\t"\$2"\\t.\\t"\$3"\\t"\$4"\\t.\\t.\\t"\$5";FISHER=" ;
            for(i=6;i<=NF;i++){print \$i ; if(i < NF){print "|"}} ;
            print "\\n"
        }
    ' ${fisher} >> res_fisher.vcf
    bgzip -f -l 9 res_fisher.vcf > res_fisher.vcf.gz # htslib command, -l 9 best compression, -c to standard output, -f to force without asking
    tabix -p vcf res_fisher.vcf.gz # htslib command
    """
}


process tsv_compress {
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}", mode: 'copy', overwrite: false
    cache 'true'

    //no channel input here for the vcf, because I do not transform it
    input:
    file tsv from fisher_ch5
    // see the scope for the use of affected_patients which is already a variable from .config file

    output:
    file "res_fisher.*"


    script:
    """
    #!/bin/bash -ue
    gzip -cf9 ${tsv} > res_fisher.tsv.gz # htslib command, -l 9 best compression, -c to standard output, -f to force without asking
    """
    // write ${} between "" to make a single argument when the variable is made of several values separated by a space. Otherwise, several arguments will be considered
}


process Backup {
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir "${out}/reports", mode: 'copy', overwrite: false // since I am in mode copy, all the output files will be copied into the publishDir. See \\wsl$\Ubuntu-20.04\home\gael\work\aa\a0e9a739acae026fb205bc3fc21f9b
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
    echo -e "full .nextflow.log is in: ${launchDir}\\nThe one in the result folder is not complete (miss the end)" > Log_info.txt
    """
}


//////// end Processes
