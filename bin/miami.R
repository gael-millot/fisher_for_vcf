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
        "chr", 
        "y.lim1", 
        "y.lim2",
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

# fisher <- "C:/Users/gael/Documents/Git_projects/08002_bourgeron/dataset/fisher.tsv"
# chr <- "C:/Users/gael/Documents/Git_projects/08002_bourgeron/dataset/hg19_grch37p5_chr_size_cumul.txt"
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
    "fisher", 
    "chr", 
    "y.lim1", 
    "y.lim2",
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


################ local function: package import


# R Packages required
req.package.list <- c(
    "lubridate", 
    "ggplot2", 
    "qqman"
)
for(i in 1:length(req.package.list)){suppressMessages(library(req.package.list[i], character.only = TRUE))}
# fun_pack(req.package = req.package.list, load = TRUE, lib.path = NULL) # packages are imported even if inside functions are written as package.name::function() in the present code


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
tempo <- fun_check(data = chr, class = "vector", typeof = "character", length = 1) ; eval(ee)
# tempo <- fun_check(data = cute, class = "vector", typeof = "character", length = 1) ; eval(ee) # check above
if(length(y.lim1) != 1 & any(grepl(y.lim1, pattern = "\\D"))){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN global_logo.R:\nTHE y_lim1 PARAMETER MUST BE A SINGLE INTEGER\nHERE IT IS: \n", paste0(y.lim1, collapse = " "))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}else{
    y.lim1 <- as.integer(y.lim1)
}
if(length(y.lim2) != 1 & any(grepl(y.lim2, pattern = "\\D"))){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN global_logo.R:\nTHE y_lim2 PARAMETER MUST BE A SINGLE INTEGER\nHERE IT IS: \n", paste0(y.lim2, collapse = " "))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}else{
    y.lim2 <- as.integer(y.lim2)
}
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
    "chr", 
    "y.lim1", 
    "y.lim2",
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN miami.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
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


################################ End pre-ignition checking


################################ Main code


################ Ignition


fun_report(data = paste0("\n\n################################################################ miami PROCESS\n\n"), output = log, path = "./", overwrite = FALSE)
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


obs <- read.table(fisher, sep = "\t", stringsAsFactors = FALSE, header = TRUE, comment.char = "")
chr <- read.table(chr, sep = "\t", stringsAsFactors = FALSE, header = TRUE, comment.char = "")


################ end Data import


############ modifications of imported tables


if(length(obs) > 0 & nrow(obs) > 0){
    names(obs)[names(obs) == "X.log10.pval."] <- "neg.log10.p"
    names(obs)[names(obs) == "AN"] <- "Nb_of_indiv"
    if(any(grepl(x = obs$chr, pattern = "chr"))){
        obs$chr <- sub(x = obs$chr, pattern = "chr", replacement = "")
    }
    if(any(grepl(x = obs$chr, pattern = "^X$"))){
        obs$chr[grepl(x = obs$chr, pattern = "^X$")] <- "23"
    }
    if(any(grepl(x = obs$chr, pattern = "^Y$"))){
        obs$chr[grepl(x = obs$chr, pattern = "^Y$")] <- "24"
    }
    if(any(grepl(x = obs$chr, pattern = "^MT$|^M$"))){
        obs$chr[grepl(x = obs$chr, pattern = "^MT$|^M$")] <- "25"
    }
    if(any( ! grepl(x = obs$chr, pattern = "\\d"))){
        tempo.cat <- paste0("ERROR IN miami.R:\nTHE chr COLUMN of the fisher.tsv FILE HAS LETTERS IN IT, OTHER THAN X, Y and MT:\n", paste0(obs$chr[grepl(x = obs$chr, pattern = "^\\d")], collapse = "\n"))
        stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
    }else{
        obs$chr <- as.integer(obs$chr)
    }
    # add the chr info to obs
    obs <- data.frame(obs, coord = vector(mode = "numeric", length = nrow(obs)))
    for(i1 in chr$CHR_NB){
        obs$coord[obs$chr == i1] <- obs$position[obs$chr == i1] + chr$LENGTH_CUMUL_TO_ADD[i1]
    }
}else{
    tempo.warn <- paste0("EMPTY fisher FILE: NO PLOT DRAWN")
    fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
    warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
}


############ end modifications of imported tables


############ plotting


#fun_open(width = 12, height = 4, pdf.name = paste0("plot_read_length_", kind)) # must be systematically opened for main.nf
png(filename = paste0("miami.png"), width = 3600, height = 1800, units = "px", res = 300)


if(length(obs) > 0 & nrow(obs) > 0){
    marging <- 1e7
    tempo.gg.name <- "gg.indiv.plot."
    tempo.gg.count <- 0
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot(obs, aes_string(x = "coord", y = "neg.log10.p")))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), geom_point(aes(color=as.factor(chr)), alpha=0.5, size=1))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), scale_color_manual(values = rep(c("grey", "skyblue"), 25)))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), scale_x_continuous(
        expand = c(0, 0), # remove space after after axis limits
        oob = scales::rescale_none,
        label = chr$CHR_NAME, 
        breaks= chr$CHR_NAME_POS, 
        limits = c(-marging, max(chr$LENGTH_CUMUL) + marging)
    ))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), scale_y_continuous(
        expand = c(0, 0), # remove space after after axis limits
        limits = c(0, y.lim1), # NA indicate that limits must correspond to data limits but ylim() already used
        oob = scales::rescale_none, 
        trans = "identity" # equivalent to ggplot2::scale_y_reverse() but create the problem of y-axis label disappearance with y.lim decreasing. Thus, do not use. Use ylim() below and after this
    ))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), theme_bw())
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank()
    ))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), annotate(
        geom = "segment", 
        x = c(0, chr$LENGTH_CUMUL), 
        xend = c(0, chr$LENGTH_CUMUL), 
        y = 0, 
        yend = -0.05, 
        size = 1
    ))


    fin.plot1 <- suppressMessages(suppressWarnings(eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + ")))))

    tempo.gg.name2 <- "gg.indiv.plot."
    tempo.gg.count2 <- 0
    assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), ggplot(obs, aes_string(x = "coord", y = "OR")))
    assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), geom_point(aes(color=as.factor(chr)), alpha=0.5, size=1))
    assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), scale_color_manual(values = rep(c("grey", "skyblue"), 25)))
    assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), scale_x_continuous(
        expand = c(0, 0), # remove space after after axis limits
        oob = scales::rescale_none,
        label = chr$CHR_NAME, 
        breaks= chr$CHR_NAME_POS, 
        limits = c(-marging, max(chr$LENGTH_CUMUL) + marging)
    ))
    assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), scale_y_continuous(
        expand = c(0, 0), # remove space after after axis limits
        limits = c(y.lim2, 0), # NA indicate that limits must correspond to data limits but ylim() already used
        oob = scales::rescale_none, 
        trans = "reverse" # equivalent to ggplot2::scale_y_reverse() but create the problem of y-axis label disappearance with y.lim decreasing. Thus, do not use. Use ylim() below and after this
    ))
    assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), theme_bw())
    assign(paste0(tempo.gg.name2, tempo.gg.count2 <- tempo.gg.count2 + 1), theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
    ))
    fin.plot2 <- suppressMessages(suppressWarnings(eval(parse(text = paste(paste0(tempo.gg.name2, 1:tempo.gg.count2), collapse = " + ")))))

    suppressMessages(suppressWarnings(gridExtra::grid.arrange(fin.plot1, fin.plot2, ncol=1, nrow = 2)))

}else{
    fun_gg_empty_graph(text = "EMPTY .length FILE: NO PLOT DRAWN")
}




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
fun_report(data = paste0("\n\n################################ JOB END\n\nTIME: ", end.date, "\n\nTOTAL TIME LAPSE: ", total.lapse, "\n"), output = log, path = "./", overwrite = FALSE)


################ end Parameter printing


################################ End Main code







