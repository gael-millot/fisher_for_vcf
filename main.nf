nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     fisher_for_vcf.nf                                               ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################
*/




//////// Processes


process workflowParam { // create a file with the workflow parameters in out_path
    label 'bash'
    publishDir "${out_path}/reports", mode: 'copy', overwrite: false
    cache 'false'

    input:
    val modules

    output:
    path "Run_info.txt"

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


process header {
    label 'bcftools' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{header_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    path vcf // no parralelization
    path tbi

    output:
    path "vcf_header.txt", emit: vcf_header_ch

    script:
    """
    bcftools head ${vcf} > vcf_header.txt | tee extract_report.txt
    """
}


process vcf_info {
    label 'r_ext' // see the withLabel: bash in the nextflow config file
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{vcf_info_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out_path}", mode: 'copy', pattern: "{vcf_csq_subfield_titles.txt,vcf_info_field_titles.txt}", overwrite: false
    cache 'true'

    input:
    path vcf_header // no parall
    path cute

    output:
    path "vcf_info_field_titles.txt", emit: vcf_info_field_titles_ch
    path "vcf_csq_subfield_titles.txt", emit: vcf_csq_subfield_titles_ch
    path "vcf_info_report.txt"

    script:
    """
    #!/bin/bash -ue
    vcf_header.R ${vcf_header} "${cute}" "vcf_info_report.txt"
    """
}


process extract {
    label 'bcftools' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{extract_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    path vcf // no parralelization
    path tbi
    val region

    output:
    path "extracted.vcf", emit: extracted_vcf_ch
    path "extracted.vcf.tbi", emit: extracted_tbi_ch

    script:
    """
    #!/bin/bash -ue
    AVAIL_THREADS=\$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print \$4}')
    if [[ "${region}" == "none" ]] ; then
        ln -s ${vcf} "extracted.vcf" | tee extract_report.txt # warning can still be zipped
        ln -s ${tbi} "extracted.vcf.tbi" # tbi always exists (as NULL if no compression)
    else
        bcftools filter --regions ${region} --write-index ${vcf} --output extracted.vcf | tee extract_report.txt
    fi
    """
    // bcftools view --no-header --threads \${AVAIL_THREADS} --write-index --output-type z --output extracted.vcf.zip ${vcf} | tee extract_report.txt # --write-index cannot be used if \${vcf} is output of slivar
    // | bcftools view --no-header --threads \${AVAIL_THREADS} --write-index --output-type z --output extracted.vcf.zip - | tee extract_report.txt # --write-index cannot be used if \${vcf} is output of slivar
}


process vcf_split {
    label 'bcftools' // see the withLabel: bash in the nextflow config file
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{vcf_split_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    path vcf // no parall
    val thread_nb

    output:
    path "split_vcf*", emit: split_vcf_ch // parallel
    path "vcf_split_report.txt"

    script:
    """
    #!/bin/bash -ue
    if [[ "${thread_nb}" == "NULL" ]] ; then
        ln -s ${vcf} "split_vcf" | tee vcf_split_report.txt
        ln -s ${tbi} "split_vcf.tbi" # tbi always exists (as NULL if no compression)
    else
        TOTAL_LINE_NB=\$(bcftools query -l ${vcf} | wc -l)
        echo -e "TOTAL LINE NUMBER IN FILE: \${TOTAL_LINE_NB}\\n\\n" > vcf_split_report.txt
        LINE_NB=\$(( (\$TOTAL_LINE_NB+(${thread_nb}-1))/${thread_nb} ))
        echo -e "NB OF FILES REQUIRED BY THE USER (thread_nb PARAMETER OF THE nextflow.config FILE) : ${thread_nb}\\n\\n" >> vcf_split_report.txt
        echo -e "LINE NUMBER PER FILE (AROUND): \${LINE_NB}\\n\\n" >> vcf_split_report.txt
        zcat ${vcf} | split --lines=\$LINE_NB --additional-suffix ".vcf" - "split_vcf" | tee -a vcf_split_report.txt 
        # --number=l/\${thread_nb} splits the files in \${thread_nb} files for parral without cutting inside lines (it is L, not one). The hyphen is for the piped file. --number=l/${thread_nb} does not work when using pipe because split needs to know the line number before splitting, which is not possible with pipe. See https://www.gnu.org/software/coreutils/manual/html_node/split-invocation.html#split-invocation
        # --filter='pigz' --additional-suffix ".vcf.gz" to rezip the file, --filter='pigz' is faster than --filter='gzip'
    fi
    """
}


process fisher {
    label 'python' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{fisher_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    path vcf // parallelization expected
    path header
    path ped
    path vcf_info_field_titles
    path vcf_csq_subfield_titles
    val tsv_extra_fields
    val filter_indiv_DP
    val filter_indiv_GQ
    val thread_nb

    output:
    path "*.tsv", emit: fisher_ch1 // multi channel
    path "*.txt"

    script:
    """
    #!/bin/bash -ue
    if grep -Eq "FISHER" ${vcf_info_field_titles} ; then # does not work without [[]]
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nTHE VCF FILE HAS ALREADY FISHER COMPUTATION PERFORMED, AS FIELDS ARE:\\n\$(cat ${vcf_info_field_titles})\\n\\n========\\n\\n"
        exit 1
    fi
    if [[ "${thread_nb}" == "NULL" ]] ; then
        echo -e "BEFORE: \$(zcat ${vcf} | wc -l)\\n"
        add_fisher.py ${vcf} ${ped} ${vcf_info_field_titles} "${tsv_extra_fields}" ${vcf_csq_subfield_titles} ${filter_indiv_DP} ${filter_indiv_GQ} | tee fisher_report.txt
        echo -e "AFTER: \$(cat fisher.tsv | wc -l)\\n"
    else
        echo -e "BEFORE: \$(cat ${vcf} | wc -l)\\n"
        cat ${header} ${vcf} > ./TEMPO.txt # assembling header and each chunk of vcf, as required by VCF tool
        add_fisher.py ./TEMPO.txt ${ped} ${vcf_info_field_titles} "${tsv_extra_fields}" ${vcf_csq_subfield_titles} ${filter_indiv_DP} ${filter_indiv_GQ} | tee fisher_report.txt
        echo -e "AFTER: \$(cat fisher.tsv | wc -l)\\n"
    fi
    """
}



process tsv_compress {
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}", mode: 'copy', overwrite: false
    cache 'true'

    //no channel input here for the vcf, because I do not transform it
    input:
    path tsv
    // see the scope for the use of affected_patients which is already a variable from .config file

    output:
    path "res_fisher.*"

    script:
    """
    #!/bin/bash -ue
    gzip -cf9 ${tsv} > res_fisher.tsv.gz # htslib command, -l 9 best compression, -c to standard output, -f to force without asking
    """
    // write ${} between "" to make a single argument when the variable is made of several values separated by a space. Otherwise, several arguments will be considered
}

// bcftools annotate if we want the fisher inside a vcf
// https://www.biostars.org/p/122690/



process miami_plot {
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}", mode: 'copy', pattern: "{*.png}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{miami_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    path fisher
    path chr
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
    path cute

    output:
    path "*.png"
    path "miami_report.txt"

    script:
    """
    #!/bin/bash -ue
    miami.R ${fisher} ${chr} "${x_lim_val}" "${vgrid}" "${top_y_column}" "${bottom_y_column}" "${color_column}" "${dot_border_color}" "${y_lim1}" "${y_lim2}" "${y_reverse1}" "${y_reverse2}" "${y_threshold1}" "${y_threshold2}" "${y_log1}" "${y_log2}" "${cute}" "miami_report.txt"
    """
}



process backup {
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', overwrite: false // since I am in mode copy, all the output files will be copied into the publishDir. See \\wsl$\Ubuntu-20.04\home\gael\work\aa\a0e9a739acae026fb205bc3fc21f9b
    cache 'false'

    input:
    path config_file
    path log_file

    output:
    path "${config_file}" // warning message if we use file config_file
    path "${log_file}" // warning message if we use file log_file
    path "Log_info.txt"

    script:
    """
    echo -e "full .nextflow.log is in: ${launchDir}\\nThe one in the result folder is not complete (miss the end)" > Log_info.txt
    """
}


//////// end Processes


//////// Workflow


workflow {

    //////// Options of nextflow run

    print("\n\nINITIATION TIME: ${workflow.start}")

    //////// end Options of nextflow run


    //////// Options of nextflow run

    // --modules (it is just for the process workflowParam)
    params.modules = "" // if --module is used, this default value will be overridden
    // end --modules (it is just for the process workflowParam)

    //////// end Options of nextflow run


    //////// Variables

    modules = params.modules // remove the dot -> can be used in bash scripts
    config_file = file("${projectDir}/nextflow.config") // file() create a path object necessary o then create the file
    log_file = file("${launchDir}/.nextflow.log")

    // from parameters (options of the nexflow command line)
    modules = params.modules // remove the dot -> can be used in bash scripts
    // end from parameters (options of the nexflow command line)


    //////// end Variables


    //////// Checks

    if( ! (sample_path.class == String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN repertoire_profiler.config FILE:\n${sample_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(sample_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }else if(sample_path =~ /.*\.gz$/){
        if( ! (file("${sample_path}.tbi").exists()) ){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID .tbi FILE ASSOCIATED TO sample_path PARAMETER IN nextflow.config FILE: ${sample_path}.tbi\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\nOTHERWISE, USE tabix -p vcf <NAME>.vcf TO INDEX THE .gz FILE\n\n========\n\n"
        }else{
            tbi_file = file("${sample_path}.tbi")
        }
    }else{
        tbi_file = file("NULL")
    }
    if( ! (ped_path.class == String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ped_path PARAMETER IN repertoire_profiler.config FILE:\n${ped_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(ped_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ped_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${ped_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (chr_path.class == String ) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID chr_path PARAMETER IN repertoire_profiler.config FILE:\n${chr_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(chr_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID chr_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${chr_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! region.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID region PARAMETER IN nextflow.config FILE:\n${region}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! tsv_extra_fields.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tsv_extra_fields PARAMETER IN nextflow.config FILE:\n${tsv_extra_fields}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! filter_indiv_DP.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID filter_indiv_DP PARAMETER IN nextflow.config FILE:\n${filter_indiv_DP}\nMUST BE A SINGLE CHARACTER STRING OF A FLOAT\n\n========\n\n"
    }else if( ! (filter_indiv_DP =~ /^[0123456789]+\.*[0123456789]*$/)){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID filter_indiv_DP PARAMETER IN nextflow.config FILE:\n${filter_indiv_DP}\nMUST BE A SINGLE CHARACTER STRING OF A FLOAT\n\n========\n\n"
    }
    if( ! filter_indiv_GQ.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID filter_indiv_GQ PARAMETER IN nextflow.config FILE:\n${filter_indiv_GQ}\nMUST BE A SINGLE CHARACTER STRING OF A FLOAT\n\n========\n\n"
    }else if( ! (filter_indiv_GQ =~ /^[0123456789]+\.*[0123456789]*$/)){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID filter_indiv_GQ PARAMETER IN nextflow.config FILE:\n${filter_indiv_GQ}\nMUST BE A SINGLE CHARACTER STRING OF A FLOAT\n\n========\n\n"
    }
    if( ! thread_nb.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID thread_nb PARAMETER IN nextflow.config FILE:\n${thread_nb}\nMUST BE A SINGLE CHARACTER STRING OF AN INTEGER VALUE\n\n========\n\n"
    }else if( ! (thread_nb =~ /^[0123456789]+$/ || thread_nb == "NULL")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID thread_nb PARAMETER IN nextflow.config FILE:\n${thread_nb}\nMUST BE A SINGLE CHARACTER STRING OF AN INTEGER VALUE\n\n========\n\n"
    }
    if( ! miami_plot =~ /^(TRUE)|(FALSE)$/ ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID miami_plot PARAMETER IN nextflow.config FILE:\n${miami_plot}\nMUST BE \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! x_lim.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID x_lim PARAMETER IN nextflow.config FILE:\n${x_lim}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! vgrid =~ /^(TRUE)|(FALSE)$/ ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID vgrid PARAMETER IN nextflow.config FILE:\n${vgrid}\nMUST BE \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! top_y_column.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID top_y_column PARAMETER IN nextflow.config FILE:\n${top_y_column}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! bottom_y_column.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID bottom_y_column PARAMETER IN nextflow.config FILE:\n${bottom_y_column}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! color_column.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID color_column PARAMETER IN nextflow.config FILE:\n${color_column}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! y_lim1.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_lim1 PARAMETER IN nextflow.config FILE:\n${y_lim1}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (y_lim1 =~ /^\-{0,1}[0123456789]+\.*[0123456789]* \-{0,1}[0123456789]+\.*[0123456789]*$/ || y_lim1 == "NULL")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_lim1 PARAMETER IN nextflow.config FILE:\n${y_lim1}\nMUST BE A SINGLE CHARACTER STRING OF TWO FLOATS SEPARATED BY A SINGLE SPACE\n\n========\n\n"
    }
    if( ! y_lim2.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_lim2 PARAMETER IN nextflow.config FILE:\n${y_lim2}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (y_lim2 =~ /^\-{0,1}[0123456789]+\.*[0123456789]* \-{0,1}[0123456789]+\.*[0123456789]*$/ || y_lim2 == "NULL")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_lim2 PARAMETER IN nextflow.config FILE:\n${y_lim2}\nMUST BE A SINGLE CHARACTER STRING OF TWO FLOATS SEPARATED BY A SINGLE SPACE\n\n========\n\n"
    }
    if( ! y_reverse1 =~ /^(TRUE)|(FALSE)$/ ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_reverse1 PARAMETER IN nextflow.config FILE:\n${y_reverse1}\nMUST BE \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! y_reverse2 =~ /^(TRUE)|(FALSE)$/ ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_reverse2 PARAMETER IN nextflow.config FILE:\n${y_reverse2}\nMUST BE \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! y_threshold1.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_threshold1 PARAMETER IN nextflow.config FILE:\n${y_threshold1}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (y_threshold1 =~ /^\-{0,1}[0123456789]+\.*[0123456789]*$/ || y_threshold1 == "NULL")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_threshold1 PARAMETER IN nextflow.config FILE:\n${y_threshold1}\nMUST BE A SINGLE CHARACTER STRING OF A FLOAT\n\n========\n\n"
    }
    if( ! y_threshold2.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_threshold2 PARAMETER IN nextflow.config FILE:\n${y_threshold2}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (y_threshold2 =~ /^\-{0,1}[0123456789]+\.*[0123456789]*$/ || y_threshold2 == "NULL")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_threshold2 PARAMETER IN nextflow.config FILE:\n${y_threshold2}\nMUST BE A SINGLE CHARACTER STRING OF A FLOAT\n\n========\n\n"
    }
    if( ! y_log1 =~ /^(TRUE)|(FALSE)$/ ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_log1 PARAMETER IN nextflow.config FILE:\n${y_log1}\nMUST BE \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! y_log2 =~ /^(TRUE)|(FALSE)$/ ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID y_log2 PARAMETER IN nextflow.config FILE:\n${y_log2}\nMUST BE \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if(miami_plot == "TRUE" || region != "none"){
        if( ! (cute_path.class == String) ){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN repertoire_profiler.config FILE:\n${cute_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
        }else if( ! (file(cute_path).exists()) ){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${cute_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
        }
    }

    // below : those variable are already used in the config file. Thus, to late to check them. And not possible to check inside the config file
    // system_exec
    // out_ini
    print("\n\nRESULT DIRECTORY: ${out_path}")
    if("${system_exec}" != "local"){
        print("    queue: ${queue}")
        print("    qos: ${qos}")
        print("    add_options: ${add_options}")
    }
    print("\n\n")


    //////// end Checks


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

    vcf_ch = Channel.fromPath(sample_path) 

    if(region != 'none'){
        if(region =~ / /){ 
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID region PARAMETER IN nextflow.config FILE:\n${region}\nMUST NOT CONTAIN SPACE. EXAMPLE OF WRITING: \"chr7:0-10000,chr10\"\n\n========\n\n"
        }
        if(region =~ /,/){ // split if comma present
            region = region.split(",")
            tempo = region.replaceAll(':.+,', ',') // for checking
        }else{
            tempo = region
        }
        tempo2 = tempo.replaceAll(':.+$', '')
        if( ! tempo2 =~ /(,)|(chr(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|Y|X|M))/){ // I have added both notations "chr23", "chr24", "chr25" or "chrY", "chrX", "chrM" because either one or the other can be used in a VCF file 
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID region PARAMETER IN nextflow.config FILE:\n${region}\nMUST BE A SINGLE CHARACTER STRING CONTAINING ONLY THESE CHROMO WRITING + OPTIONAL POSITIONS\nchr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chrY, chrX, chrM\n\n========\n\n"
        }
    }



    //////// end Channels


    //////// files import

    // in variable because a single file. If "NULL", will create a empty file, present in work folders, but that cannot be correctly linked. Thus, if the file has to be redirected into a channel inside a process, it will not work. Thus, in the first process using meta_file, I hard copy the NULL file if required (see below)
    ped_file = file(ped_path) // in variable because a single file
    chr_file = file(chr_path) // in variable because a single file
    cute_file = file(cute_path) // in variable because a single file

    //////// end files import



    //////// Main

    workflowParam(
        modules
    )

    header(
        vcf_ch, 
        tbi_file
    )

    vcf_info(
        header.out.vcf_header_ch,
        cute_file
    )

    if(region != "none"){
        extract(
            vcf_ch, 
            tbi_file, 
            region
        )
        extracted_vcf_ch2 = extract.out.extracted_vcf_ch
        extracted_tbi_ch2 = extract.out.extracted_tbi_ch
    }else{
        extracted_vcf_ch2 = vcf_ch
    }

extracted_vcf_ch2.view()

    if(thread_nb != "NULL"){
        vcf_split(
            extracted_vcf_ch2,
            thread_nb
        )
        split_ch = vcf_split.out.split_vcf_ch.flatten() // to convert the list into multiple channels
    }else{
        split_ch = extracted_vcf_ch2
    }

    fisher(
        split_ch, 
        header.out.vcf_header_ch.first(), 
        ped_file,
        vcf_info.out.vcf_info_field_titles_ch.first(),
        vcf_info.out.vcf_csq_subfield_titles_ch.first(),
        tsv_extra_fields, 
        filter_indiv_DP, 
        filter_indiv_GQ, 
        thread_nb
    )

    fisher_ch2 = fisher.out.fisher_ch1.collectFile(name: "fisher.tsv", keepHeader: true, skip: 1)

    tsv_compress(
        fisher_ch2
    )

    if(miami_plot == "TRUE"){
        miami_plot(
            fisher_ch2,
            chr_file,
            x_lim_val,
            vgrid,
            top_y_column,
            bottom_y_column,
            color_column,
            dot_border_color,
            y_lim1,
            y_lim2,
            y_reverse1,
            y_reverse2,
            y_threshold1,
            y_threshold2,
            y_log1,
            y_log2,
            cute_file
        )
    }

    backup(
        config_file, 
        log_file
    )


}
