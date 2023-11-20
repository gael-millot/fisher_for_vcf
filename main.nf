nextflow.enable.dsl=2
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


process vcf_subfield_title {
    label 'r_ext' // see the withLabel: bash in the nextflow config file
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{*.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    path vcf // no parall
    path cute

    output:
    path "vcf_info_field_titles.txt", emit: vcf_info_field_titles_ch
    path "vcf_csq_subfield_titles.txt", emit: vcf_csq_subfield_titles_ch
    path "vcf_subfield_title_report.txt"

    script:
    """
    #!/bin/bash -ue
    vcf_subfield_title.R ${vcf} "${cute}" "vcf_subfield_title_report.txt"
    """
}


process fisher {
    label 'python' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{fisher_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    tuple val(region2), path(vcf) // parallelization expected for each value of region_ch
    path ped
    path tbi
    path vcf_info_field_titles
    path vcf_csq_subfield_titles
    val tsv_extra_fields

    output:
    path "*.tsv", emit: fisher_ch1 // multi channel
    path "*.txt"

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


process tsv2vcf {
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}", mode: 'copy', overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'true'

    input:
    path vcf
    path fisher

    output:
    path "res_fisher.*"

    script:
    """
    #!/bin/bash -ue
    PREHEADER='##fileformat=VCFv4.2;build by main.nf\\n##WARNING: This file is not a true VCF since FORMAT AND sample (indiv) columns are not present'
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


    //////// Variables from config.file that need to be modified


    if(x_lim == 'whole' || (x_lim == 'region' && region == 'none')){ // for the miami plot
        x_lim_val = "chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chrY, chrX, chrM" // I have added both notations "chr23", "chr24", "chr25" or "chrY", "chrX", "chrM" because either one or the other can be used in a VCF file 
    }else if(x_lim == 'region'){
        x_lim_val = region
    }else{
        x_lim_val = x_lim // value for the miami plot
    }

    //////// end Variables from config.file that need to be modified


    //////// Checks

    if( ! (sample_path in String) ){
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
    if( ! (ped_path in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ped_path PARAMETER IN repertoire_profiler.config FILE:\n${ped_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(ped_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID ped_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${ped_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (chr_path in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID chr_path PARAMETER IN repertoire_profiler.config FILE:\n${chr_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(chr_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID chr_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${chr_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
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
    if( ! (cute_path in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN repertoire_profiler.config FILE:\n${cute_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(cute_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${cute_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
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

    //////// Channels

    vcf_ch = Channel.fromPath(sample_path) 

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

    vcf_subfield_title(
        vcf_ch,
        cute_file
    )

    fisher(
        region_ch.combine(vcf_ch),
        ped_file,
        tbi_file,
        vcf_subfield_title.out.vcf_info_field_titles_ch.first(),
        vcf_subfield_title.out.vcf_csq_subfield_titles_ch.first(),
        tsv_extra_fields
    )

    fisher_ch2 = fisher.out.fisher_ch1.collectFile(name: "fisher.tsv", skip: 1, keepHeader: true)

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

    tsv2vcf(
        vcf_ch,
        fisher_ch2
    )

    tsv_compress(
        fisher_ch2
    )

    backup(
        config_file, 
        log_file
    )


}
