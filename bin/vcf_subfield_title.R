#!/usr/bin/env Rscript

#########################################################################
##                                                                     ##
##     vcf_subfield_title.R                                            ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Computational Biology Department                                ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################




################################ Aim


################################ End Aim


################################ Introduction

# https://stackoverflow.com/questions/59668347/rmarkdown-turn-off-title

################################ End Introduction


################################ Acknowlegments


################################ End Acknowlegments


################################ Initialization


# R version checking
if(version$version.string != "R version 4.1.2 (2021-11-01)"){
    stop(paste0("\n\n================\n\nERROR IN vcf_subfield_title.R\n", version$version.string, " IS NOT THE 4.1.2 RECOMMANDED\n\n================\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "vcf_subfield_title"


################################ End Initialization


################################ Parameters that need to be set by the user


################################ End Parameters that need to be set by the user


################################ Config import


tempo.cat <- "KIND OF RUN (SCRIPT, COPY-PASTE OR SOURCE): "
if(interactive() == FALSE){ # if(grepl(x = commandArgs(trailingOnly = FALSE), pattern = "R\\.exe$|\\/R$|Rcmd\\.exe$|Rcmd$|Rgui\\.exe$|Rgui$|Rscript\\.exe$|Rscript$|Rterm\\.exe$|Rterm$")){ # detection of script usage
    run.way <- "SCRIPT"
    cat(paste0("\n\n", tempo.cat, run.way, "\n"))
    command <- paste0(commandArgs(trailingOnly = FALSE), collapse = ",") # recover the full command
    args <- commandArgs(trailingOnly = TRUE) # recover arguments written after the call of the R script
    if(any(is.na(args))){
        stop(paste0("\n\n================\n\nERROR IN vcf_subfield_title.R\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "vcf",  
        "cute", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the vcf_subfield_title.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN vcf_subfield_title.R\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
    }
    for(i1 in 1:length(tempo.arg.names)){
        assign(tempo.arg.names[i1], args[i1])
    }
    rm(tempo.arg.names, args, i1)
}else if(sys.nframe() == 0L){ # detection of copy-paste/direct execution (for debugging). With script it is also 0, with source, it is 4
    run.way <- "COPY-PASTE"
    cat(paste0("\n\n", tempo.cat, run.way, "\n"))
}else{
    run.way <- "SOURCE" # using source(), sys.nframe() is 4
    cat(paste0("\n\n", tempo.cat, run.way, "\n"))
}
rm(tempo.cat)


################################ End Config import

################################ Test

# fisher <- "C:/Users/gael/Documents/Git_projects/fisher_for_vcf/dataset/fisher.tsv"
# chr.path <- "C:/Users/gael/Documents/Git_projects/fisher_for_vcf/dataset/hg19_grch37p5_chr_size_cumul.txt"
# region <- "chr1:0-50000, chr3:0-150000"
# x.lim <- "region"
# bottom.y.column <- 
# color.column <- 
# y.lim1 <- 5
# y.lim2 <- 3
# cute <- "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.4.0/cute_little_R_functions.R" 
# log <- "report.txt"



################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    if(run.way == "SCRIPT"){"command"}, 
    "vcf", 
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN vcf_subfield_title.R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN vcf_subfield_title.R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (vcf_subfield_title.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (vcf_subfield_title.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
    }
}
char.length <- nchar(param.list)
space.add <- max(char.length) - char.length + 5
param.ini.settings <- character(length = length(param.list))
for(i in 1:length(param.list)){
    param.ini.settings[i] <- paste0("\n", param.list[i], paste0(rep(" ", space.add[i]), collapse = ""), paste0(get(param.list[i]), collapse = ",")) # no env = sys.nframe(), inherit = FALSE in get() because look for function in the classical scope
}


################################ End Recording of the initial parameters


################################ Functions


# Functions are built such that they should have no direct use of Global objects (going through the R scope), and only use function arguments
# 1) Cute little function is sourced for the moment into the .GlobalEnv environment, but may be interesting to put it into a new environement just above .GlobalEnv environment. See https://stackoverflow.com/questions/9002544/how-to-add-functions-in-an-existing-environment
# 2) Argument names of each function must not be a name of Global objects (error message otherwise)
# 3) Argument name of each function ends with "_fun" in the first function, "_2fun" in the second, etc. This prevent conflicts with the argument partial names when using these functions, notably when they are imbricated


################ import functions from cute little functions toolbox

if(length(cute) != 1){
    stop(paste0("\n\n============\n\nERROR IN vcf_subfield_title.R\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN vcf_subfield_title.R\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN vcf_subfield_title.R\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else{
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN vcf_subfield_title.R: CODE HAS TO BE MODIFIED\n\n============\n\n")
    stop(tempo.cat, call. = FALSE)
}


# required cute function checking
req.function <- c(
    "fun_check", 
    "fun_pack", 
    "fun_report"
)
tempo <- NULL
for(i1 in req.function){
    if(length(find(i1, mode = "function")) == 0L){
        tempo <- c(tempo, i1)
    }
}
if( ! is.null(tempo)){
    tempo.cat <- paste0("ERROR IN vcf_subfield_title.R\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking

################ end import functions from cute little functions toolbox

################ local function: package import

# R Packages required
req.package.list <- c(
    "lubridate"
)
for(i in 1:length(req.package.list)){suppressMessages(library(req.package.list[i], character.only = TRUE))}
# fun_pack(req.package = req.package.list, load = TRUE, lib.path = NULL) # packages are imported even if inside functions are written as package.name::function() in the present code

################ end local function: package import

################ other functions

vcf.info.fields <- function(path){
    # AIM
    # scan a vcf file header and return the subfields of the INFO field (i.e., everything that starts by "##INFO=<ID=" in the vcf header)
    # ARGUMENTS
    # vcf: path of a non zipped vcf file
    # RETURN
    # A character vector of the INFO subfields
    # REQUIRED PACKAGES
    # None
    # REQUIRED FUNCTIONS FROM CUTE_LITTLE_R_FUNCTION
    # none
    # EXAMPLES
    # vcf.info.fields("C:\\Users\\gael\\Documents\\Git_projects\\fisher_for_vcf\\dataset\\Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10000.vcf")
    # DEBUGGING
    # path = "C:\\Users\\gael\\Documents\\Git_projects\\fisher_for_vcf\\dataset\\Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10000.vcf" # for function debugging
    goon <- TRUE
    count <- 0
    info.subfield.name <- NULL
    while(goon == TRUE){
        count <- count + 1
        tempo <- scan(path, what = "list", sep = '\n', skip = count, nlines = 1, quiet = TRUE)
        if(grepl(tempo, pattern = "^##INFO=<ID=.*")){
            tempo.right <- sub(x = tempo, pattern = "^##INFO=<ID=", replacement = "")
            tempo.left <- sub(x = tempo.right, pattern = ",.*$", replacement = "")
            info.subfield.name <- c(info.subfield.name, tempo.left)
        }
        if( ! grepl(tempo, pattern = "^#.*")){ # because a vcf always starts by #
            goon <- FALSE
            if(is.null(info.subfield.name)){
                stop(paste0("\n\n============\n\nERROR IN the vcf.csq.subfields function of vcf_subfield_title.R\nNO \"##INFO=<ID=\" DESCRIPTION FOUND IN THE HEADER OF:\n", path, "\n\n============\n\n"), call. = FALSE)
            }
        }
        # add special split for CSQ
    }
    # add check with bottom.y.column
    return(info.subfield.name)
}

vcf.csq.subfields <- function(path){
    # AIM
    # scan a vcf file header and return the subfields of the INFO CSQ subfield (i.e., subfields described in "##INFO=<ID=CSQ" in the vcf header)
    # ARGUMENTS
    # vcf: path of a non zipped vcf file
    # RETURN
    # A character vector of the INFO subfields
    # REQUIRED PACKAGES
    # None
    # REQUIRED FUNCTIONS FROM CUTE_LITTLE_R_FUNCTION
    # none
    # EXAMPLES
    # vcf.info.fields("C:\\Users\\gael\\Documents\\Git_projects\\fisher_for_vcf\\dataset\\Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10000.vcf")
    # DEBUGGING
    # path = "C:\\Users\\gael\\Documents\\Git_projects\\fisher_for_vcf\\dataset\\Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10000.vcf" # for function debugging
    goon <- TRUE
    count <- 0
    info.subfield.name <- NULL
    while(goon == TRUE){
        count <- count + 1
        tempo <- scan(path, what = "list", sep = '\n', skip = count, nlines = 1, quiet = TRUE)
        if(grepl(tempo, pattern = "^##INFO=<ID=CSQ.*")){
            tempo.right <- sub(x = tempo, pattern = "^##INFO=<ID=.*Format: ", replacement = "")
            tempo.left <- sub(x = tempo.right, pattern ='">$', replacement = "")
            info.subfield.name <- unlist(strsplit(tempo.left, split = "\\|"))
        }
        if( ! grepl(tempo, pattern = "^#.*")){ # because a vcf always starts by #
            goon <- FALSE
            if(is.null(info.subfield.name)){
                stop(paste0("\n\n============\n\nERROR IN the vcf.csq.subfields function of vcf_subfield_title.R\nNO \"##INFO=<ID=CSQ\" DESCRIPTION FOUND IN THE HEADER OF:\n", path, "\n\n============\n\n"), call. = FALSE)
            }
        }
        # add special split for CSQ
    }
    # add check with bottom.y.column
    return(info.subfield.name)
}


################ end other functions

################################ End Functions


################################ Pre-ignition checking


# reserved words
# end reserved words
# argument primary checking
arg.check <- NULL #
text.check <- NULL #
checked.arg.names <- NULL # for function debbuging: used by r_debugging_tools
ee <- expression(arg.check <- c(arg.check, tempo$problem) , text.check <- c(text.check, tempo$text) , checked.arg.names <- c(checked.arg.names, tempo$object.name))
tempo <- fun_check(data = vcf, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = log, class = "vector", typeof = "character", length = 1) ; eval(ee)
if(any(arg.check) == TRUE){ # normally no NA
    stop(paste0("\n\n================\n\n", paste(text.check[arg.check], collapse = "\n"), "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between == #
}
# end argument primary checking
# second round of checking and data preparation
# management of NA arguments
# end management of NA arguments
# management of NULL arguments
tempo.arg <-c(
    "vcf",
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN vcf_subfield_title.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of NULL arguments
# code that protects set.seed() in the global environment
# end code that protects set.seed() in the global environment
# warning initiation
ini.warning.length <- options()$warning.length
options(warning.length = 8170)
warn <- NULL
# warn.count <- 0 # not required
# end warning initiation
# other checkings
# end other checkings
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking

if( ! file.exists(vcf)){
    stop(paste0("\n\n============\n\nERROR IN vcf_subfield_title.R\nFILE INDICATED IN THE vcf PARAMETER DOES NOT EXISTS: ", vcf, "\n\n============\n\n"), call. = FALSE)
}


################################ End pre-ignition checking


################################ Main code

res1 <- vcf.info.fields(vcf)
res2 <- vcf.csq.subfields(vcf)
cat(res1, file= "./vcf_info_field_titles.txt", append = FALSE, sep =" ") # add a sep
cat(res2, file= "./vcf_csq_subfield_titles.txt", append = FALSE, sep =" ") # add a sep


################################ end Main code


################ Ignition


fun_report(data = paste0("\n\n################################################################ vcf_subfield_title PROCESS\n\n"), output = log, path = "./", overwrite = FALSE)
ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING\n\n"), output = log, path = "./", overwrite = FALSE)


################ End ignition


############ main




############ end main


################ Seeding inactivation


set.seed(NULL)


################ end Seeding inactivation


################ Environment saving


fun_report(data = paste0("\n\n################################ RUNNING END"), output = log, path = "./", overwrite = FALSE)
end.date <- Sys.time()
end.time <- as.numeric(end.date)
total.lapse <- round(lubridate::seconds_to_period(end.time - ini.time))
fun_report(data = paste0("\n\nEND TIME: ", end.date), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\nTOTAL TIME LAPSE: ", total.lapse), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\nALL DATA SAVED IN all_objects.RData"), output = log, path = "./", overwrite = FALSE)


################ end Environment saving


################ Warning messages


fun_report(data = paste0("\n\n################################ RECAPITULATION OF WARNING MESSAGES"), output = log, path = "./", overwrite = FALSE)
if( ! is.null(warn)){
    fun_report(data = paste0("\n\n", warn), output = log, path = "./", overwrite = FALSE)
}else{
    fun_report(data = paste0("\n\nNO WARNING MESSAGE TO REPORT"), output = log, path = "./", overwrite = FALSE)
}


################ end Warning messages


################ Parameter printing


fun_report(data = paste0("\n\n################################ INITIAL SETTINGS OF PARAMETERS"), output = log, path = "./", overwrite = FALSE)
fun_report(data = param.ini.settings, output = log, path = "./", overwrite = FALSE, , vector.cat = TRUE)
fun_report(data = paste0("\n\n################################ R SYSTEM AND PACKAGES"), output = log, path = "./", overwrite = FALSE)
tempo <- sessionInfo()
tempo$otherPkgs <- tempo$otherPkgs[order(names(tempo$otherPkgs))] # sort the packages
tempo$loadedOnly <- tempo$loadedOnly[order(names(tempo$loadedOnly))] # sort the packages
fun_report(data = tempo, output = log, path = "./", overwrite = FALSE, , vector.cat = TRUE)
fun_report(data = paste0("\n\n################################ JOB END\n\nTIME: ", end.date, "\n\nTOTAL TIME LAPSE: ", total.lapse, "\n"), output = log, path = "./", overwrite = FALSE)


################ end Parameter printing


################################ End Main code







