#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#########################################################################
##                                                                     ##
##     fisher_lod.py                                                   ##
##                                                                     ##
##     Freddy Cliquet                                                  ##
##     GHFC                                                            ##
##     Neuroscience Department                                         ##
##     Institut Pasteur Paris                                          ##
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


################################ End Introduction


################################ Acknowlegments


################################ End Acknowlegments


################################ Initialization


import csv
import sys
import re
import os.path
from cyvcf2 import VCF, Writer
import scipy.stats as stats
import pandas as pd
import numpy as np
import random # for shuffle() function
import warnings # for warnings.catch_warnings()


################################ End Initialization


################################ Parameters that need to be set by the user



################################ End Parameters that need to be set by the user


################################ Config import


vcf_path = sys.argv[1]  # 1st argument: pedigree file name, 1 ID par ligne, sys.argv takes arguments from the bash line command when running a .py script
ped_path = sys.argv[2]
region = sys.argv[3]
vcf_info_field_titles_path = sys.argv[4]
tsv_extra_fields = sys.argv[5]
vcf_csq_subfield_titles_path = sys.argv[6]
fisher_report = sys.argv[7]


################################ End Config import

################################ Test


# vcf_path="/pasteur/zeus/projets/p01/BioIT/gmillot/08002_bourgeron/dataset/Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10000.vcf.gz"
# vcf_path="/mnt/c/Users/gael/Documents/Git_projects/08002_bourgeron/dataset/Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10000.vcf.gz"
# ped_path="/pasteur/zeus/projets/p01/BioIT/gmillot/08002_bourgeron/dataset/Dyslexia.pedigree.txt" # functions for slivar
# ped_path="/mnt/c/Users/gael/Documents/Git_projects/08002_bourgeron/dataset/Dyslexia.pedigree.txt"
# region=[‘chr7:0-147000000’, ‘chr10:1000000-2000000’]


################################ End Test

################################ Recording of the initial parameters




################################ End Recording of the initial parameters


################################ Functions


def fisher(v, status, tsv_columns, tsv_extra_fields_wo_csq, csq_subfield_name, csq_subfield_pos):
    '''
    AIM
        compute fisher for a line of a vcf
    WARNINGS
    ARGUMENTS
        v: a single line of a vcf object from VCF()
        tsv_columns: column names of the final tsv file
        tsv_extra_fields_wo_csq: subfields of INFO field of the vcf to add as column in the tsv (excluding potential CSQ subfields like Polyphen)
        csq_subfield_name: subfields of the CSQ subfield of the INFO field of the vcf to add as column in the tsv
        csq_subfield_pos: positions of csq_subfield_name in the all the subfields of CSQ, indicated by vcf_csq_subfield_titles
    RETURN
        a data frame with a single row
    REQUIRED PACKAGES
        None
    EXAMPLE
        for v in vcf:
            fisher(v = v, tsv_columns = tsv_columns, tsv_extra_fields = tsv_extra_fields)
    DEBUGGING
        # use the container: sudo docker run -ti --entrypoint bash -v /mnt/c/Users/gael/Documents/Git_projects/fisher_for_vcf/dataset:/tmp gmillot/python_v3.9.10_extended_v3.1:gitlab_v8.7)
vcf = VCF("/tmp/Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10.vcf")
count = 1 
for v in vcf:
    if count == 1: break
status = {'C0011JY': 1, 'C0011JZ': 2, 'C0011K1': 1, 'C0011K2': 2, 'C0011K3': 2, 'C0011K5': 1, 'C0011KA': 2, 'C0011KB': 1, 'IP00FNP': 2, 'IP00FNW': 2, 'IP00FNY': 2}
tsv_columns = ['CHROM','POS','REF','ALT', 'INFO', 'GENE','IMPACT','CONSEQUENCE','AFF','UNAFF','OR','P_VALUE','NEG_LOG10_P_VALUE','PATIENT_NB', 'CSQ_TRANSCRIPT_NB']
tsv_extra_fields_wo_csq = ['AC', 'AF']
csq_subfield_name = ['PolyPhen']
csq_subfield_pos = [26]
    '''
    
    # je fais 2 dictionnaire pour stocker mes compte d'atteint/non atteint porteurs/non porteurs
    # dans chaque dictionaire j'associerai a un genotype (clef) un nombre d'individu.
    aff=dict() #dictionary: return of this is {1:0, 2:0, } with key:value each element of the dict (key -> integers by default)
    una=dict()
    df3 = pd.DataFrame(columns = tsv_columns + csq_subfield_name + tsv_extra_fields_wo_csq)


    # je traite tous les variants qui ont un champ CSQ (annotation VEP)
    if v.INFO.get('CSQ') is not None: # https://www.w3schools.com/python/ref_dictionary_get.asp
        # j'initialise une variable qui met permet de garder trace du nombre d'individus pour laquelle j'ai des information pour le variant v courant
        an=0 # nb of indiv
        # pour le variant actuel on parcours la liste des individus (iid) avec les genotype associé (gt), la depth (dp) et la genotyping quality (gq)
        # gt_types: method of VCF that specify the kind of genotype -> written into gt with 0=HOM_REF ; 1=HET ; 2=UNKNOWN ; 3=HOM_ALT, see https://brentp.github.io/cyvcf2/docstrings.html
        for iid, gt, dp, gq in zip(vcf.samples, v.gt_types, v.format('DP'), v.format('GQ')):
            # je filtre un depth>=10 et une genotyping quality >=30
            # on peut decider de baisser le gq a 20 par exemple, c'est juste une suggestion
            # each line of zip() is processed
            with warnings.catch_warnings():
                # https://docs.python.org/3/library/warnings.html#warning-categories
                warnings.filterwarnings("ignore", category=FutureWarning)
                if dp not in ['.', ''] and int(dp)>=10 and gq not in ['.', ''] and float(gq)>=30:
                    # je met a jour le compte d'individus dans mes dictionnaire aff et una en fonction du phenotype de l'individu courant iid
                    if status[iid]==2: # 2 == status "aff" (result from ped)
                        aff[gt]=aff.get(gt,0)+1 # means if the key gt exists in the aff dict, return the value of aff[gt] + 1, otherwise return 0 + 1
                    if status[iid]==1: # 1 == status "unaff" (result from ped)
                        una[gt]=una.get(gt,0)+1
                    an+=1
        # une fois que l'on a lu les information pour tous les individus, nous calculons le Fisher
        # ici c'est porteur (gt 1 ou 3) versus non porteur (gt 0) pour les atteints (aff) versus les non atteint (una)
        oddsratio, pvalue = stats.fisher_exact([[aff.get(1,0)+aff.get(3,0),aff.get(0,0)],[una.get(1,0)+una.get(3,0),una.get(0,0)]])

        # filling a one row data frame with or without adding
        tempo_csq = v.INFO.get('CSQ').split(',') # number of fields in CSQ (comma sep), i.e., nb of rows
        if len(csq_subfield_name) == 0:
            gene = tempo_csq[0].split('|')[3] # See protocole 109: gene taken in position 4 if no SYMBOL field
            consequence = tempo_csq[0].split('|')[1]
            impact = tempo_csq[0].split('|')[2]
            df2=pd.DataFrame([[v.CHROM, v.POS, v.REF, ''.join(v.ALT), ';'.join([i3[0]+"="+str(i3[1]) for i3 in v.INFO]), gene, impact, consequence, aff, una, oddsratio, pvalue, -np.log10(pvalue), an]], columns = tsv_columns)
            if len(tsv_extra_fields_wo_csq) > 0: # add extra columns coming from tsv_extra_fields into the tsv file
                for i4 in tsv_extra_fields_wo_csq:
                    df2[i4] = v.INFO.get(i4)

            df3 = df3.append(df2)

        # get the csq_subfield_name info and if more than 1 -> several lines, per each value, else
        else:
            subfield_pos = list(range(0, len(tempo_csq)))
            gene = []
            impact = []
            consequence = []
            for i3 in csq_subfield_name: # number of subfields of CSQ wanted, i.e., number of columns
                locals()[i3] = []
            for i3 in tempo_csq: # tempo_csq
                gene.append(i3.split('|')[3]) # See protocole 109: gene taken in position 4 if no SYMBOL field
                impact.append(i3.split('|')[2])
                consequence.append(i3.split('|')[1])

                for i4 in list(range(0, len(csq_subfield_name))): # for each CSQ_ field wanted (polyphen, SIFT)
                    # f = open(fisher_report, "a") ; f.write(str(csq_subfield_name) + '\n' + str(i4) + '\n\n\n\n' + str(csq_subfield_name[i4]) + '\n\n\n\n' +  str(csq_subfield_pos) + '\n')
                    locals()[csq_subfield_name[i4]].append(i3.split('|')[csq_subfield_pos[i4]])

            # f = open("test.txt", "a") ; f.write(str(locals()) + '\n\n\n\n' + str(locals()[csq_subfield_name[i4]]) + '\n')
            # from here, we have three lists gene, impact and consequence with n elements, and i csq_subfield_name lists with n elements
            # we can consider the nb of elements as lines (because future lines in the data frame) and the lists as columns
            # removal of "lines" with no values in all the csq_subfield_name
            tempo_empty_log = [] # True if all are empty. Warning, works by lines
            for i3 in list(range(0, len(tempo_csq))):
                tempo_log = []
                for i4 in csq_subfield_name:
                    if locals()[i4][i3] == "":
                        tempo_log.append(True)
                    else:
                        tempo_log.append(False)
                if all([i4 == True for i4 in tempo_log]) is True:
                    tempo_empty_log.append(True)
                else:
                    tempo_empty_log.append(False)

            if any([i3 == False for i3 in tempo_empty_log]): # means at least one non empty
                subfield_pos = [subfield_pos[i4] for i4 in list(range(0, len(subfield_pos))) if tempo_empty_log[i4] == False]
                gene = [gene[i4] for i4 in list(range(0, len(gene))) if tempo_empty_log[i4] == False]
                impact = [impact[i4] for i4 in list(range(0, len(impact))) if tempo_empty_log[i4] == False]
                consequence = [consequence[i4] for i4 in list(range(0, len(consequence))) if tempo_empty_log[i4] == False]
                for i4 in csq_subfield_name: # for each CSQ_ field wanted (polyphen, SIFT)
                    # f = open("test.txt", "a") ; f.write(str(locals()[i4]))
                    tempo = []
                    for i5 in list(range(0, len(locals()[i4]))):
                        if tempo_empty_log[i5] == False:
                            tempo.append(locals()[i4][i5])

                    # f = open("test.txt", "a") ; f.write(str(tempo))
                    locals()[i4] = tempo # [locals()[i4][i5] for i5 in list(reversed(range(0, len(locals()[i4])))) if tempo_empty_log[i5] == False]
                # end removal of "lines" with no values in all the csq_subfield_name
            else: # means everything is empty, thus, I have to shrink these lists
                subfield_pos = ['NA']
                gene = [gene[0]]
                impact = [impact[0]]
                consequence = [consequence[0]]
                for i4 in csq_subfield_name: # for each CSQ_ field wanted (polyphen, SIFT)
                    locals()[i4] = ['NA']

            # f = open("test.txt", "a") ; f.write(str(list(range(0, len(subfield_pos)))))
            for i3 in list(range(0, len(subfield_pos))):
                # f = open("test.txt", "a") ; f.write(str([[v.CHROM, v.POS, v.REF, ''.join(v.ALT), ';'.join([i4[0]+"="+str(i4[1]) for i4 in v.INFO]), gene[i3], impact[i3], consequence[i3], aff, una, oddsratio, pvalue, -np.log10(pvalue), an, subfield_pos[i3]]]) + '\n\n\n\n' + str(i3) + '\n\n\n\n' + str(subfield_pos) + '\n\n\n\n' + str(impact) + '\n\n\n\n')
                tempo = pd.DataFrame([[v.CHROM, v.POS, v.REF, ''.join(v.ALT), ';'.join([i4[0]+"="+str(i4[1]) for i4 in v.INFO]), gene[i3], impact[i3], consequence[i3], aff, una, oddsratio, pvalue, -np.log10(pvalue), an, subfield_pos[i3]]], columns = tsv_columns)

                # f = open("test.txt", "a") ; f.write(str(csq_subfield_name) + '\n')
                for i4 in csq_subfield_name: # for each CSQ_ field wanted (polyphen, SIFT), column added to the data frame
                    # f = open("test.txt", "a") ; f.write(str(csq_subfield_name) + '\n\n\n\n' + str(csq_subfield_pos) + '\n\n\n\n' + str(locals()[i4][i3]) + '\n\n\n\n' + str(i3) + '\n\n\n\n' + str(subfield_pos) + '\n\n\n\n')
                    tempo[i4] = locals()[i4][i3]

                if len(tsv_extra_fields_wo_csq) > 0: # add extra columns coming from tsv_extra_fields into the tsv file
                    for i6 in tsv_extra_fields_wo_csq:
                        tempo[i6] = v.INFO.get(i6)

                if i3 == 0:
                    df2 = tempo

                else:
                    df2 = df2.append(tempo)

            df3 = df3.append(df2)

    return df3


################################ End Functions


################################ Pre-ignition checking


# reserved words
# end reserved words
# argument primary checking
if isinstance(vcf_path, str) is not True:
    sys.exit("\n\n========\n\nError in fisher_lod.py: the vcf_path must be a single character string: \n"+" ".join(vcf_path)+"\n\n========\n\n")

if isinstance(ped_path, str) is not True:
    sys.exit("\n\n========\n\nError in fisher_lod.py: the ped_path must be a single character string: \n"+" ".join(ped_path)+"\n\n========\n\n")

if isinstance(region, str) is not True:
    sys.exit("\n\n========\n\nError in fisher_lod.py: the region must be a single character string: \n"+" ".join(region)+"\n\n========\n\n")

if isinstance(vcf_info_field_titles_path, str) is not True:
    sys.exit("\n\n========\n\nError in fisher_lod.py: the vcf_info_field_titles_path must be a single character string: \n"+" ".join(vcf_info_field_titles_path)+"\n\n========\n\n")

if isinstance(tsv_extra_fields, str) is not True:
    sys.exit("\n\n========\n\nError in fisher_lod.py: the tsv_extra_fields must be a single character string: \n"+" ".join(tsv_extra_fields)+"\n\n========\n\n")

if isinstance(vcf_csq_subfield_titles_path, str) is not True:
    sys.exit("\n\n========\n\nError in fisher_lod.py: the vcf_csq_subfield_titles_path must be a single character string: \n"+" ".join(vcf_csq_subfield_titles_path)+"\n\n========\n\n")

if isinstance(fisher_report, str) is not True:
    sys.exit("\n\n========\n\nError in fisher_lod.py: the fisher_report must be a single character string: \n"+" ".join(fisher_report)+"\n\n========\n\n")
# end argument primary checking
# second round of checking and data preparation
# management of NA arguments
# end management of NA arguments
# management of NULL arguments
# end management of NULL arguments
# code that protects set.seed() in the global environment
random.seed(1)
# end code that protects set.seed() in the global environment
# warning initiation
# warn.count <- 0 # not required
# end warning initiation
# other checkings
if os.path.exists(vcf_path) is not True:
    sys.exit("\n\n========\n\nError in fisher_lod.py: the vcf_path file does not exist: \n"+" ".join(vcf_path)+"\n\n========\n\n")

if os.path.exists(ped_path) is not True:
    sys.exit("\n\n========\n\nError in fisher_lod.py: the ped_path file does not exist: \n"+" ".join(ped_path)+"\n\n========\n\n")

if os.path.exists(vcf_info_field_titles_path) is not True:
    sys.exit("\n\n========\n\nError in fisher_lod.py: the vcf_info_field_titles_path file does not exist: \n"+" ".join(vcf_info_field_titles_path)+"\n\n========\n\n")

tsv_extra_fields = tsv_extra_fields.split(' ')

if os.path.exists(vcf_csq_subfield_titles_path) is not True:
    sys.exit("\n\n========\n\nError in fisher_lod.py: the vcf_csq_subfield_titles_path file does not exist: \n"+" ".join(vcf_csq_subfield_titles_path)+"\n\n========\n\n")
# end other checkings
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking


################################ End pre-ignition checking


################################ Main code


################ Ignition

# if region == "None":
#    region = None
#else:
#    # region = region.strip('[]').replace('"', '').replace(' ', '').split(',') # old version when no parall in nf
#    region = region.replace('"', '').replace(' ', '')

################ End ignition


################ Graphical parameter initialization


################ End graphical parameter initialization


################ Data import


# on crée un dictionnaire qui associe le phenotype a chaque barcode, ici c'est fait un lisant un fichier pedigree
status = dict()
with open(ped_path, 'r') as pin:
    pin = csv.reader(pin, delimiter='\t')
    for row in pin:
        status[row[1]]=int(row[5])
# le dictionnaire ressemble a ca :
# {'IP00FN5': 2, 'IP00FLK': 2, 'IP00FLT': 1, 'IP00FM2': 2, 'IP00FMC': 1}
# f = open("status.txt", "a") ; f.write(str(status))
# le fichier vcf contenant les data des individus que l'on etudie
# si bgzippé (.gz), doit avoir un index tabix (.tbi)
vcf = VCF(vcf_path) # import the VCF as it is. But this tool as the advantage to extract easily, using INFO, etc., see https://brentp.github.io/cyvcf2/docstrings.html
# below is to visualize the vcf if required
# for v in vcf:
#     w = Writer('./vcf.tsv', vcf) # https://brentp.github.io/cyvcf2/docstrings.html
#     w.write_record(v)
#     w.close()
# the INFO field of the VCF isis a complex object


with open(vcf_info_field_titles_path, 'r') as f:
    vcf_info_field_titles = f.readlines()[0].split(' ')

with open(vcf_csq_subfield_titles_path, 'r') as f:
    vcf_csq_subfield_titles = f.readlines()[0].split(' ')

################ end Data import


############ modifications of imported tables

# checking that everything is fine with tsv_extra_fields, and info recovering
tsv_columns=['CHROM','POS','REF','ALT', 'INFO', 'GENE_EXAMPLE','IMPACT_EXAMPLE','CONSEQUENCE_EXAMPLE','AFF','UNAFF','OR','P_VALUE','NEG_LOG10_P_VALUE','PATIENT_NB'] #warning : can be replaced below
csq_subfield_name = []
csq_subfield_pos = []
tsv_extra_fields_wo_csq = []
if all([i0 == 'NULL' for i0 in tsv_extra_fields]) is False:
    tempo_log = [bool(re.search("^CSQ_.*$", i1)) for i1 in tsv_extra_fields] # is there any CSQ_ ?
    if any(i1 for i1 in tempo_log) is True:
        tsv_columns=['CHROM','POS','REF','ALT', 'INFO', 'GENE','IMPACT','CONSEQUENCE','AFF','UNAFF','OR','P_VALUE','NEG_LOG10_P_VALUE','PATIENT_NB', 'CSQ_TRANSCRIPT_NB']
        tempo_pos = []
        for i2 in list(range(0, len(tsv_extra_fields))):
            if bool(re.search("^CSQ_.*$", tsv_extra_fields[i2])) is True:
                tempo = re.sub(pattern = "^CSQ_", repl = "", string = tsv_extra_fields[i2])
                csq_subfield_name.append(tempo)
            else:
                tempo_pos.append(i2)

        if all([i1 in vcf_csq_subfield_titles for i1 in csq_subfield_name]) is not True:
            sys.exit("\n\n========\n\nError in fisher_lod.py: some of the CSQ subfield of the tsv_extra_fields parameter (starting by CSQ_): \n"+" ".join(csq_subfield_name)+"\nare not in the vcf_csq_subfield_titles parameter: \n"+" ".join(vcf_csq_subfield_titles)+"\n\n========\n\n")
        else:
            for i2 in csq_subfield_name:
                for i3 in list(range(0, len(vcf_csq_subfield_titles))):
                    if i2 == vcf_csq_subfield_titles[i3]:
                        csq_subfield_pos.append(i3)

        tsv_extra_fields_wo_csq = [tsv_extra_fields[i1] for i1 in tempo_pos]
        if all([i1 in vcf_info_field_titles for i1 in tsv_extra_fields_wo_csq]) is not True:
            sys.exit("\n\n========\n\nError in fisher_lod.py: not considering CSQ_, some of the tsv_extra_fields parameter values: \n"+" ".join(tsv_extra_fields)+"\nare not in the vcf_info_field_titles parameter: \n"+" ".join(vcf_info_field_titles)+"\n\n========\n\n")

    elif all([i1 in vcf_info_field_titles for i1 in tsv_extra_fields]) is not True:
            sys.exit("\n\n========\n\nError in fisher_lod.py: some of the tsv_extra_fields parameter values: \n"+" ".join(tsv_extra_fields)+"\nare not in the vcf_info_field_titles parameter: \n"+" ".join(vcf_info_field_titles)+"\n\n========\n\n")

    else:
        tsv_extra_fields_wo_csq = tsv_extra_fields




df = pd.DataFrame(columns = tsv_columns + csq_subfield_name + tsv_extra_fields_wo_csq)
with warnings.catch_warnings():
    # https://docs.python.org/3/library/warnings.html#warning-categories
    warnings.filterwarnings("ignore", category=FutureWarning)
#    if region is None:
#        for v in vcf:
#            tempo = fisher(v = v, columns = columns)
#            df = df.append(tempo)
#    else:
#        for i1 in region :
#            for v in vcf(i1):
#                tempo = fisher(v = v, columns = columns)
#                df = df.append(tempo)
    LOG_EVERY_N = 1000
    count = 0
    report = open(fisher_report, "a")
    for v in vcf: # parse each line of vcf with the content of region in it
        count = count + 1
        if v.CHROM == region:
            tempo = fisher(v = v, status = status, tsv_columns = tsv_columns, tsv_extra_fields_wo_csq = tsv_extra_fields_wo_csq, csq_subfield_name = csq_subfield_name, csq_subfield_pos = csq_subfield_pos)
            df = df.append(tempo)

        if (count % LOG_EVERY_N) == 0:
            report.write("NUMBER OF VCF LINES SCANNED: " + str(count) + "\n")

    report.close()

# on ecrit la dataframe dans un fichier
df.to_csv('./fisher.tsv', sep='\t', index=False)



############ end modifications of imported tables


############ plotting


############ end plotting


################ Pdf window closing


################ end Pdf window closing


################ Seeding inactivation
random.seed()

################ end Seeding inactivation


################ Environment saving


################ end Environment saving


################ Warning messages


################ end Warning messages


################ Parameter printing


################ end Parameter printing


################################ End Main code





