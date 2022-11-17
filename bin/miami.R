#!/usr/bin/env Rscript

#########################################################################
##                                                                     ##
##     miami.R                                                         ##
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
    stop(paste0("\n\n================\n\nERROR IN miami.R\n", version$version.string, " IS NOT THE 4.1.2 RECOMMANDED\n\n================\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "miami"


################################ End Initialization


################################ Parameters that need to be set by the user


whole <- "chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chrY, chrX, chrM" # true values of x_lim if x_lim = "whole", as imported through nextflow


################################ End Parameters that need to be set by the user


################################ Config import


tempo.cat <- "KIND OF RUN (SCRIPT, COPY-PASTE OR SOURCE): "
if(interactive() == FALSE){ # if(grepl(x = commandArgs(trailingOnly = FALSE), pattern = "R\\.exe$|\\/R$|Rcmd\\.exe$|Rcmd$|Rgui\\.exe$|Rgui$|Rscript\\.exe$|Rscript$|Rterm\\.exe$|Rterm$")){ # detection of script usage
    run.way <- "SCRIPT"
    cat(paste0("\n\n", tempo.cat, run.way, "\n"))
    command <- paste0(commandArgs(trailingOnly = FALSE), collapse = ",") # recover the full command
    args <- commandArgs(trailingOnly = TRUE) # recover arguments written after the call of the R script
    if(any(is.na(args))){
        stop(paste0("\n\n================\n\nERROR IN miami.R\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "fisher", 
        "chr.path", 
        "x.lim", 
        "vgrid", 
        "top.y.column",
        "bottom.y.column",
        "color.column",
        "dot.border.color", 
        "y.lim1", 
        "y.lim2",
        "reverse1", 
        "reverse2", 
        "y.threshold1", 
        "y.threshold2", 
        "y.log1", 
        "y.log2", 
        "cute", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the miami.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN miami.R\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
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
# x.lim <- "chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chrY, chrX, chrM" ### "chr1:0-50000, chr3:0-150000"
# vgrid <- "TRUE"
# top.y.column <- "NEG_LOG10_P_VALUE"
# bottom.y.column <- "AF"
# color.column <- "NULL"
# dot.border.color <- "white"
# y.lim1 <- "NULL"
# y.lim2 <- "NULL"
# reverse1 <- "FALSE"
# reverse2 <- "TRUE"
# y.threshold1 <- "1.2"
# y.threshold2 <- "0.5"
# y.log1 <- "FALSE"
# y.log2 <- "FALSE"
# cute <- "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.4.0/cute_little_R_functions.R" 
# log <- "miami_report.txt"


################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    if(run.way == "SCRIPT"){"command"}, 
    "whole", 
    "fisher", 
    "chr.path", 
    "x.lim", 
    "vgrid", 
    "top.y.column", 
    "bottom.y.column",
    "color.column",
    "dot.border.color", 
    "y.lim1", 
    "y.lim2", 
    "reverse1", 
    "reverse2", 
    "y.threshold1", 
    "y.threshold2",
    "y.log1", 
    "y.log2", 
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN miami.R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN miami.R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (miami.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (miami.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
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
    stop(paste0("\n\n============\n\nERROR IN miami.R\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN miami.R\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN miami.R\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else{
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN miami.R: CODE HAS TO BE MODIFIED\n\n============\n\n")
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
    tempo.cat <- paste0("ERROR IN miami.R\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking

################ end import functions from cute little functions toolbox

################ local function: package import

# R Packages required
req.package.list <- c(
    "lubridate", 
    "ggplot2", 
    "scales", 
    "grid", 
    "qqman"
)
for(i in 1:length(req.package.list)){suppressMessages(library(req.package.list[i], character.only = TRUE))}
# fun_pack(req.package = req.package.list, load = TRUE, lib.path = NULL) # packages are imported even if inside functions are written as package.name::function() in the present code

################ end local function: package import

################ other functions

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
tempo <- fun_check(data = fisher, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = chr.path, class = "vector", typeof = "character", length = 1) ; eval(ee)
# tempo <- fun_check(data = cute, class = "vector", typeof = "character", length = 1) ; eval(ee) # check above
if(all(x.lim != "NULL")){
    tempo <- fun_check(data = x.lim, class = "vector", typeof = "character", length = 1) ; eval(ee)
}else{
    x.lim <- NULL
}
tempo <- fun_check(data = vgrid, class = "vector", typeof = "character", length = 1) ; eval(ee)
if(all(top.y.column != "NULL")){
    tempo <- fun_check(data = top.y.column, class = "vector", typeof = "character", length = 1) ; eval(ee)
}else{
    top.y.column <- NULL
}
if(all(bottom.y.column != "NULL")){
    tempo <- fun_check(data = bottom.y.column, class = "vector", typeof = "character", length = 1) ; eval(ee)
}else{
    bottom.y.column <- NULL
}
if(all(color.column != "NULL")){
    tempo <- fun_check(data = color.column, class = "vector", typeof = "character", length = 1) ; eval(ee)
}else{
    color.column <- NULL
}
if(all(dot.border.color != "NULL")){
    tempo <- fun_check(data = dot.border.color, class = "vector", typeof = "character", length = 1) ; eval(ee)
}else{
    dot.border.color <- NULL
}
if(all(y.lim1 != "NULL")){
    tempo <- fun_check(data = y.lim1, class = "vector", typeof = "character", length = 1) ; eval(ee)
}else{
    y.lim1 <- NULL
}
if(all(y.lim2 != "NULL")){
    tempo <- fun_check(data = y.lim2, class = "vector", typeof = "character", length = 1) ; eval(ee)
}else{
    y.lim2 <- NULL
}
tempo <- fun_check(data = reverse1, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = reverse2, class = "vector", typeof = "character", length = 1) ; eval(ee)
if(all(y.threshold1 != "NULL")){
    tempo <- fun_check(data = y.threshold1, class = "vector", typeof = "character", length = 1) ; eval(ee)
}else{
    y.threshold1 <- NULL
}
if(all(y.threshold2 != "NULL")){
    tempo <- fun_check(data = y.threshold2, class = "vector", typeof = "character", length = 1) ; eval(ee)
}else{
    y.threshold2 <- NULL
}
tempo <- fun_check(data = y.log1, class = "vector", typeof = "character", length = 1) ; eval(ee)
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
    "fisher",
    "chr.path", 
    "vgrid", 
    "reverse1", 
    "reverse2", 
    "y.log1", 
    "y.log2", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN miami.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of NULL arguments
# management of ""
tempo.arg <-c(
    "fisher", 
    "chr.path", 
    "x.lim", 
    "vgrid", 
    "top.y.column",
    "bottom.y.column", 
    "color.column", 
    "dot.border.color", 
    "y.lim1", 
    "y.lim2", 
    "reverse1", 
    "reverse2",
    "y.threshold1", 
    "y.threshold2", 
    "y.log1", 
    "y.log2", 
    "cute", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = function(x){any(x == "")})
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN miami.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE \"\"")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of ""
# code that protects set.seed() in the global environment
# end code that protects set.seed() in the global environment
# warning initiation
ini.warning.length <- options()$warning.length
options(warning.length = 8170)
warn <- NULL
# warn.count <- 0 # not required
# end warning initiation
# other checkings
for(i0 in c("vgrid", "reverse1", "reverse2", "y.log1", "y.log2")){
    if(get(i0) == "TRUE"){
        assign(i0, TRUE)
    }else if(get(i0) == "FALSE"){
        assign(i0, FALSE)
    }else{
        tempo.cat <- paste0("ERROR IN miami.R\n", i0, " PARAMETER CAN ONLY BE \"TRUE\" OR \"FALSE\": ", get(i0))
        stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
    }
}
# end other checkings
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking


################################ End pre-ignition checking


################################ Main code


################ Ignition


fun_report(data = paste0("\n\n################################################################ miami PROCESS\n\n"), output = log, path = "./", overwrite = TRUE)
ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING\n\n"), output = log, path = "./", overwrite = FALSE)


################ End ignition


################ Graphical parameter initialization


pdf(file = NULL)
par.ini <- par(no.readonly = TRUE) # to recover the initial graphical parameters if required (reset)
invisible(dev.off()) # close the new window
zone.ini <- matrix(1, ncol=1)
if(erase.graphs == TRUE){
    graphics.off()
}else{
    tempo.warn <- paste0("GRAPHICS HAVE NOT BEEN ERASED. GRAPHICAL PARAMETERS MAY HAVE NOT BEEN REINITIALIZED")
    fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
    warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
}


################ End graphical parameter initialization


################ Data import


if( ! file.exists(fisher)){
    stop(paste0("\n\n============\n\nERROR IN miami.R\nFILE INDICATED IN THE fisher PARAMETER DOES NOT EXISTS: ", fisher, "\n\n============\n\n"), call. = FALSE)
}else{
    obs <- read.table(fisher, sep = "\t", stringsAsFactors = FALSE, header = TRUE, comment.char = "")
    if(length(obs) > 0 & nrow(obs) > 0){
        empty.obs <- FALSE
    }else{
        empty.obs <- TRUE
    }
}
if( ! file.exists(chr.path)){
    stop(paste0("\n\n============\n\nERROR IN miami.R\nFILE INDICATED IN THE chr.path PARAMETER DOES NOT EXISTS: ", chr.path, "\n\n============\n\n"), call. = FALSE)
}else{
    chr <- read.table(chr.path, sep = "\t", stringsAsFactors = FALSE, header = TRUE, comment.char = "")
}


################ end Data import


############ modifications of imported tables

xmin_plot <- 0 # coordinates for plotting 
if(length(obs) > 0 & nrow(obs) > 0){
    # names(obs)[names(obs) == "NEG_LOG10_P_VALUE"] <- "neg.log10.p"
    # names(obs)[names(obs) == "PATIENT_NB"] <- "Nb_of_indiv"
    if(any(grepl(x = obs$CHROM, pattern = "chr"))){
        obs$CHROM <- sub(x = obs$CHROM, pattern = "chr", replacement = "")
    }
    if(any(grepl(x = obs$CHROM, pattern = "^X$"))){
        obs$CHROM[grepl(x = obs$CHROM, pattern = "^X$")] <- "23"
    }
    if(any(grepl(x = obs$CHROM, pattern = "^Y$"))){
        obs$CHROM[grepl(x = obs$CHROM, pattern = "^Y$")] <- "24"
    }
    if(any(grepl(x = obs$CHROM, pattern = "^MT$|^M$"))){
        obs$CHROM[grepl(x = obs$CHROM, pattern = "^MT$|^M$")] <- "25"
    }
    if(any( ! grepl(x = obs$CHROM, pattern = "\\d"))){
        tempo.cat <- paste0("ERROR IN miami.R:\nTHE chr COLUMN of the fisher.tsv FILE HAS LETTERS IN IT, OTHER THAN X, Y and MT:\n", paste0(obs$CHROM[grepl(x = obs$CHROM, pattern = "^\\d")], collapse = "\n"))
        stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
    }else{
        obs$CHROM <- as.integer(obs$CHROM)
    }
    # add the chr info to obs
    obs <- data.frame(obs, coord = vector(mode = "numeric", length = nrow(obs)))
    for(i1 in chr$CHR_NB){
        obs$coord[obs$CHROM == i1] <- obs$POS[obs$CHROM == i1] + chr$LENGTH_CUMUL_TO_ADD[i1]
    }
    # preparation of the x coordinates: three solutions: 1) whole object (see above), 2) single chromo "chr7" or "chr7:0-15", 3) several chromo chr7, chr8" or "chr7:0-15, chr8" or "chr7:0-15, chr8:0-20"
    # The idea is to select rows of chr and potentially restrict some chr limits
    if( ! is.null(x.lim)){
        is.whole <- FALSE
        if(x.lim == whole){ #at that stage, x.lim is a single character
            is.whole <- TRUE
        }
        tempo <- strsplit(x = x.lim, split = ",")[[1]]
        tempo <- gsub(x = tempo, pattern = " ", replacement = "")
        if( ! all(grepl(x = tempo, pattern = "^chr.+"))){
            tempo.cat <- paste0("ERROR IN miami.R:\nTHE x_lim PARAMETER MUST START WITH \"chr\" IF NOT \"none\":\n", paste0(x_lim, collapse = " "))
            stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
        }
        if(any(grepl(x = tempo, pattern = ":"))){
            # means that there are coordinates
            if( ! all(grepl(tempo, pattern = "-"))){# normally no NA with is.null()
                tempo.cat <- paste0("ERROR IN miami.R:\nTHE x_lim PARAMETER MUST BE WRITTEN LIKE THIS \"chr7:0-147000000, chr10:1000000-2000000\" IF COORDINATES ARE SPECIFIED: \n", paste0(x_lim, collapse = " "))
                stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
            }
            tempo2 <- strsplit(x = tempo, split = ":")
            chr_x_lim <- sapply(X = tempo2, FUN = function(x){x[1]})
            chr_x_lim <- gsub(x = chr_x_lim, pattern = " ", replacement = "")
            coord_x_lim <- sapply(X = tempo2, FUN = function(x){x[2]})
            tempo3 <- strsplit(x = coord_x_lim, split = "-")
            xmin_x_lim <- sapply(X = tempo3, FUN = function(x){x[1]})
            xmin_x_lim <- gsub(x = xmin_x_lim, pattern = " ", replacement = "")
            xmax_x_lim <- sapply(X = tempo3, FUN = function(x){x[2]})
            xmax_x_lim <- gsub(x = xmax_x_lim, pattern = " ", replacement = "")
            if(any(grepl(xmin_x_lim, pattern = "\\D")) | any(grepl(xmax_x_lim, pattern = "\\D"))){# normally no NA with is.null()
                tempo.cat <- paste0("ERROR IN miami.R:\nTHE x_lim PARAMETER MUST BE WRITTEN LIKE THIS \"chr7:0-147000000, chr10:1000000-2000000\" IF COORDINATES ARE SPECIFIED: \n", paste0(x_lim, collapse = " "))
                stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
            }else{
                xmin_x_lim <- as.integer(xmin_x_lim)
                xmax_x_lim <- as.integer(xmax_x_lim)
                if(any(xmax_x_lim - xmin_x_lim < 0)){
                    tempo.cat <- paste0("ERROR IN miami.R:\nTHE x_lim PARAMETER MUST BE WRITTEN WITH ORDERED COORDINATES, LIKE THIS \"chr7:0-147000000, chr10:1000000-2000000\", IF COORDINATES ARE SPECIFIED: \n", paste0(x_lim, collapse = " "))
                    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
                }
            }
        }else{
            chr_x_lim <- tempo
            coord_x_lim <- NULL
            xmin_x_lim <- NULL
            xmax_x_lim <- NULL
        }
        # modification of the chr object for restricted plotting
        tempo.coord <- which(chr$CHR %in% chr_x_lim) # which rows of chr to take for plotting
        if(any(chr$BP_LENGTH[tempo.coord] - xmax_x_lim < 0)){
            tempo.cat <- paste0("ERROR IN miami.R:\nTHE x_lim PARAMETER HAS AT LEAST ONE COORDINATE THAT IS ABOVE THE MAX LENGTH OF THE CHROMO.\nCHROMO LENGTH: ", paste0(chr$BP_LENGTH[tempo.coord], collapse = " "), "\nMAX COORDINATE: ", paste0(xmax_x_lim, collapse = " "))
            stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
        }
        if(tempo.coord[1] > 1){
            xmin_plot <- chr$LENGTH_CUMUL[tempo.coord[1] - 1]
        }
        chr <- chr[tempo.coord[1]:tempo.coord[length(tempo.coord)], ]
        if( ! is.null(coord_x_lim)){
            xmin_plot <- xmin_plot + xmin_x_lim[1] # the left boundary of the plot is corrected
            chr$LENGTH_CUMUL[nrow(chr)] <- chr$LENGTH_CUMUL[nrow(chr)] - chr$BP_LENGTH[nrow(chr)] + xmax_x_lim[length(xmax_x_lim)] # the right boundary of the plot is corrected
            chr$CHR_NAME_POS <- (c(xmin_plot, chr$LENGTH_CUMUL[-nrow(chr)]) + chr$LENGTH_CUMUL) / 2 # the positions of names in the x-axis of the plot are corrected
        }
        # restriction of obs
        obs <- obs[obs$coord >= xmin_plot & obs$coord <= chr$LENGTH_CUMUL[nrow(chr)], ]
    }else{
        tempo.warn <- paste0("x.lim is NULL: NO PLOT DRAWN")
        fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
        warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
    }
}

for(i0 in c("y.lim1", "y.lim2")){
    if( ! is.null(get(i0))){
        tempo <- unlist(strsplit(x = get(i0), split = " "))
        if(length(tempo) != 2 | ! all(grepl(tempo, pattern = "^[0123456789.\\-\\+eE]*$"))){
            tempo.cat <- paste0("ERROR IN miami.R:\nTHE ", i0, " PARAMETER MUST BE TWO NUMERIC VALUES SEPARATED BY A SINGLE SPACE\nHERE IT IS: \n", paste0(get(i0), collapse = " "))
            stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
        }else{
            assign(i0, as.numeric(tempo))
        }
    }
}
for(i0 in c("y.threshold1", "y.threshold2")){
    if( ! is.null(get(i0))){
        tempo <- unlist(strsplit(x = get(i0), split = " "))
        if(length(tempo) != 1 | ! all(grepl(tempo, pattern = "^[0123456789.\\-\\+eE]*$"))){
            tempo.cat <- paste0("ERROR IN miami.R:\nTHE ", i0, " PARAMETER MUST BE TWO NUMERIC VALUES SEPARATED BY A SINGLE SPACE\nHERE IT IS: \n", paste0(get(i0), collapse = " "))
            stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
        }else{
            assign(i0, as.numeric(tempo))
        }
    }
}

############ end modifications of imported tables


############ plotting


#fun_open(width = 12, height = 4, pdf.name = paste0("plot_read_length_", kind)) # must be systematically opened for main.nf
png.size <- 1800 # px
png(filename = paste0("miami.png"), width = png.size * 2, height = png.size, units = "px", res = 300)

if(empty.obs == TRUE){
    fun_gg_empty_graph(text = paste0("NO PLOT DRAWN\nTHE region PARAMETER\nMIGHT BE OUTSIDE\nOF THE RANGE OF THE VCF FILE"))
}else if(length(obs) > 0 & nrow(obs) > 0 & ! is.null(x.lim)){
    marging <- (chr$LENGTH_CUMUL[nrow(chr)] - xmin_plot) * 0.005 # chr$LENGTH_CUMUL and xmin_plot have been corrected depending on x.lim boundaries
    y.min.pos <- ifelse(is.null(y.lim1), min(obs[ , top.y.column]), min(y.lim1))
    y.max.pos <- ifelse(is.null(y.lim1), max(obs[ , top.y.column]), max(y.lim1))
    tempo.gg.name <- "gg.indiv.plot."
    tempo.gg.count <- 0
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot(obs, aes_string(x = "coord", y = top.y.column)))
    if(vgrid){
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), geom_vline(
            xintercept = c(xmin_plot, chr$LENGTH_CUMUL),
            size = 0.25,
            color = "grey80"
        ))
    }
    if( ! is.null(y.threshold1)){
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), geom_hline(
            yintercept = y.threshold1,
            linetype = "22", 
            size = 0.25, 
            color = "red"
        ))
    }
    if(is.null(color.column)){
        if(is.null(dot.border.color)){
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), geom_point(
                aes(color = as.factor(CHROM)), 
                alpha = 1, 
                pch = 16, 
                size = 1
            ))
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), scale_color_manual(values = rep(c("grey20", "skyblue"), 25)))
        }else{
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), geom_point(
                aes(fill = as.factor(CHROM)), 
                alpha = 1, 
                color = dot.border.color, 
                pch = 21, 
                size = 1
            ))
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), scale_fill_manual(values = rep(c("grey20", "skyblue"), 25)))
        }
    }else{
        if(is.null(dot.border.color)){
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), geom_point(
                aes(color = color.column), 
                alpha = 1, 
                pch = 16, 
                size = 1
            ))
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), scale_color_gradient2())
        }else{
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), geom_point(
                aes(fill = as.factor(CHROM)), 
                alpha = 1, 
                color = dot.border.color, 
                pch = 21, 
                size = 1
            ))
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), scale_fill_gradient2())
        }
    }
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::ggtitle(
        paste0("x.lim: ", ifelse(is.whole, "whole genome", x.lim), 
        ifelse( ! is.null(y.threshold1), paste0(", top threshold: ", y.threshold1), ""), 
        ifelse( ! (is.null(y.threshold2) & is.null(bottom.y.column)), paste0(", bottom threshold: ", y.threshold2), ""), 
        ifelse(y.log1, ", top y-axis: log10", ""), ifelse(y.log2, ", bottom y-axis: log10", ""))
    ))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), scale_x_continuous(
        name = "CHR", 
        expand = c(0, 0), # remove space after after axis limits
        oob = scales::rescale_none,
        label = chr$CHR_NAME, 
        breaks= chr$CHR_NAME_POS, 
        limits = c(xmin_plot - marging, max(chr$LENGTH_CUMUL) + marging)
    ))
    if(y.log1){
        if(any(obs[ , top.y.column] <= 0)){
            tempo.cat <- paste0("ERROR IN miami.R:\nTHE y_log1 PARAMETER CANNOT BE SET TO \"TRUE\" IF 0 OR NEG VALUES IN THE ", top.y.column, " FIELD OF THE TSV OR VCF")
            stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
        }else{
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), scale_y_continuous(
                expand = c(0, 0), # remove space after after axis limits
                limits = if(reverse1){c(y.max.pos, y.min.pos)}else{c(y.min.pos, y.max.pos)}, # NA indicate that limits must correspond to data limits but ylim() already used
                oob = scales::rescale_none, 
                trans = "log10", 
                breaks = scales::trans_breaks("log10", function(x){10^x}), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))
            ))
            # assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), annotation_logticks(outside = TRUE))
        }
    }else{
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), scale_y_continuous(
            expand = c(0, 0), # remove space after after axis limits
            limits = if(reverse1){c(y.max.pos, y.min.pos)}else{c(y.min.pos, y.max.pos)}, # NA indicate that limits must correspond to data limits but ylim() already used
            oob = scales::rescale_none, 
            trans = "identity"
        ))
    }
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), theme_bw())
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), theme(
        plot.title = ggplot2::element_text(size = 8), 
        legend.position=if(is.null(color.column)){"none"}, 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y.left = element_line(size = 0.25), 
        # axis.line.x.bottom = element_line(size = 0.25), # ugly, thus i added geom_hline below
        axis.line.y.left = element_line(size = 0.25) 
    ))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), geom_hline(
        yintercept = ifelse(
            is.null(y.lim1), 
            ifelse(reverse1, max(obs[ , top.y.column]), min(obs[ , top.y.column])), 
            ifelse(reverse1, max(y.lim1), min(y.lim1))
        ), 
        size = 0.25
    ))
    # add tick lines if vgrid is FALSE
    if( ! vgrid){
        gline = linesGrob(y = c(-0.02, 0),  gp = gpar(col = "black", lwd = 0.5))
        for(i2 in c(xmin_plot, chr$LENGTH_CUMUL)){
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::annotation_custom(gline, xmin = i2, xmax = i2, ymin = -Inf, ymax = Inf))
        }
    }


    fin.plot1 <- suppressMessages(suppressWarnings(eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + ")))))
    # tempo.output <- ggplot2::ggplot_build(fin.plot1)

    if(is.null(bottom.y.column)){
        suppressMessages(suppressWarnings(gridExtra::grid.arrange(fin.plot1, ncol=1, nrow = 1)))
    }else{
        y.min.pos2 <- ifelse(is.null(y.lim2), min(obs[ , bottom.y.column]), min(y.lim2))
        y.max.pos2 <- ifelse(is.null(y.lim2), max(obs[ , bottom.y.column]), max(y.lim2))
        tempo.gg.name2 <- "gg.indiv.plot."
        tempo.gg.count2 <- 0
        assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), ggplot(obs, aes_string(x = "coord", y = bottom.y.column)))
        if(vgrid){
            assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), geom_vline(
                xintercept = c(xmin_plot, chr$LENGTH_CUMUL),
                size = 0.25,
                color = "grey80"
            ))
        }
        if( ! is.null(y.threshold2)){
            assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), geom_hline(
                yintercept = y.threshold2,
                linetype = "22", 
                size = 0.25, 
                color = "red"
            ))
        }
        if(is.null(color.column)){
            if(is.null(dot.border.color)){
                assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), geom_point(
                    aes(color = as.factor(CHROM)), 
                    alpha = 1, 
                    pch = 16, 
                    size = 1
                ))
                assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), scale_color_manual(values = rep(c("grey20", "skyblue"), 25)))
            }else{
                assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), geom_point(
                    aes(fill = as.factor(CHROM)), 
                    alpha = 1, 
                    color = dot.border.color, 
                    pch = 21, 
                    size = 1
                ))
                assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), scale_fill_manual(values = rep(c("grey20", "skyblue"), 25)))
            }
        }else{
            if(is.null(dot.border.color)){
                assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), geom_point(
                    aes(color = color.column), 
                    alpha = 1, 
                    pch = 16, 
                    size = 1
                ))
                assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), scale_color_gradient2())
            }else{
                assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), geom_point(
                    aes(fill = as.factor(CHROM)), 
                    alpha = 1, 
                    color = dot.border.color, 
                    pch = 21, 
                    size = 1
                ))
                assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), scale_fill_gradient2())
            }
        }
        assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), scale_x_continuous(
            expand = c(0, 0), # remove space after after axis limits
            oob = scales::rescale_none,
            label = chr$CHR_NAME, 
            breaks= chr$CHR_NAME_POS, 
            limits = c(xmin_plot - marging, max(chr$LENGTH_CUMUL) + marging)
        ))
        if(y.log2){
            if(any(obs[ , bottom.y.column] <= 0)){
                tempo.cat <- paste0("ERROR IN miami.R:\nTHE y_log2 PARAMETER CANNOT BE SET TO \"TRUE\" IF 0 OR NEG VALUES IN THE ", bottom.y.column, " FIELD OF THE TSV OR VCF")
                stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
            }else{
                assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), scale_y_continuous(
                    expand = c(0, 0), # remove space after after axis limits
                    limits = if(reverse2){c(y.min.pos2, y.max.pos2)}else{c(y.max.pos2, y.min.pos2)}, # NA indicate that limits must correspond to data limits but ylim() already used
                    oob = scales::rescale_none, 
                    trans = "log10", 
                    breaks = scales::trans_breaks("log10", function(x){10^x}), 
                    labels = scales::trans_format("log10", scales::math_format(10^.x))
                ))
                # assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), annotation_logticks(outside = TRUE)) # 
            }
        }else{
            assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), scale_y_continuous(
                expand = c(0, 0), # remove space after after axis limits
                limits = if(reverse2){c(y.min.pos2, y.max.pos2)}else{c(y.max.pos2, y.min.pos2)}, # NA indicate that limits must correspond to data limits but ylim() already used
                oob = scales::rescale_none, 
                trans = "identity" # equivalent to ggplot2::scale_y_reverse() but create the problem of y-axis label disappearance with y.lim decreasing. Thus, do not use. Use ylim() below and after this
            ))
        }
        assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), theme_bw())
        assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), theme(
            legend.position=if(is.null(color.column)){"none"},
            panel.border = element_blank(),
            panel.grid = element_blank(), 
            axis.ticks.y.left = element_line(size = 0.25), 
            # axis.line.x.top = element_line(size = 0.25), # is not displayed. Thus, I add a geom_hline below
            axis.line.y.left = element_line(size = 0.25),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), 
        ))
        # add x-axis line
        assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), geom_hline(
            yintercept = ifelse(
                is.null(y.lim2), 
                ifelse(reverse2, max(obs[ , bottom.y.column]), min(obs[ , bottom.y.column])), 
                ifelse(reverse2, max(y.lim2), min(y.lim2))
            ),
            size = 0.25
        ))
        # add tick lines if vgrid is FALSE
        if( ! vgrid){
            gline = linesGrob(y = c(1, 1.02),  gp = gpar(col = "black", lwd = 0.5))
            for(i2 in c(xmin_plot, chr$LENGTH_CUMUL)){
                assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), ggplot2::annotation_custom(gline, xmin = i2, xmax = i2, ymin = -Inf, ymax = Inf))
            }
        }

        fin.plot2 <- suppressMessages(suppressWarnings(eval(parse(text = paste(paste0(tempo.gg.name2, 1:tempo.gg.count2), collapse = " + ")))))
        gl <- lapply(list(fin.plot1, fin.plot2), ggplotGrob)  
        wd <- do.call(unit.pmax, lapply(gl, "[[", 'widths'))
        gl <- lapply(gl, function(x){x[['widths']] = wd ; x})
        if( ! vgrid){
            gl[[1]]$layout$clip[gl[[1]]$layout$name=="panel"] <- "off"
            gl[[2]]$layout$clip[gl[[2]]$layout$name=="panel"] <- "off"
        }
        suppressMessages(suppressWarnings(gridExtra::grid.arrange(gl[[1]], gl[[2]], ncol=1, nrow = 2)))
    }
}else{
    fun_gg_empty_graph(text = paste0("NO PLOT DRAWN\nTHE x_lim PARAMETER\nMIGHT BE OUTSIDE\nOF THE RANGE OF THE VCF FILE\nOR THE RANGE OF THE region PARAMETER\nOR NULL"))
} #else already dealt above




############ end plotting


################ Pdf window closing


graphics.off()


################ end Pdf window closing


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



################ end Parameter printing


################################ End Main code







