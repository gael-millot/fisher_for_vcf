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
ped = sys.argv[2]
region = sys.argv[3]
vcf_info_field_titles_path = sys.argv[4]
tsv_extra_fields = sys.argv[5]
vcf_csq_subfield_titles_path = sys.argv[6]


################################ End Config import

################################ Test


# vcf_path="/pasteur/zeus/projets/p01/BioIT/gmillot/08002_bourgeron/dataset/Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10000.vcf.gz"
# vcf_path="/mnt/c/Users/gael/Documents/Git_projects/08002_bourgeron/dataset/Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10000.vcf.gz"
# ped="/pasteur/zeus/projets/p01/BioIT/gmillot/08002_bourgeron/dataset/Dyslexia.pedigree.txt" # functions for slivar
# ped="/mnt/c/Users/gael/Documents/Git_projects/08002_bourgeron/dataset/Dyslexia.pedigree.txt"
# region=[‘chr7:0-147000000’, ‘chr10:1000000-2000000’]


################################ End Test

################################ Recording of the initial parameters




################################ End Recording of the initial parameters


################################ Functions


def fisher(v, tsv_columns, tsv_extra_fields_wo_csq, csq_subfield_name, csq_subfield_pos):
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
        v = vcf
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


        tempo_csq = v.INFO.get('CSQ').split(',')


################ here, get the csq_subfield_name info and if more than 1 -> several lines, per each value, else 



        if len(csq_subfield_name) == 0:
            gene = tempo_csq[0].split('|')[3] # See protocole 109: gene taken in position 4 if no SYMBOL field
            impact = tempo_csq[0].split('|')[1]
            severity = tempo_csq[0].split('|')[2]
            df2=pd.DataFrame([[v.CHROM, v.POS, v.REF, v.ALT, ';'.join([i3[0]+"="+str(i3[1]) for i3 in v.INFO]), gene, severity, impact, aff, una, oddsratio, pvalue, -np.log10(pvalue), an]], columns = tsv_columns)
            if len(tsv_extra_fields_wo_csq) > 0: # add extra columns coming from tsv_extra_fields into the tsv file
                for i4 in tsv_extra_fields_wo_csq:
                    df2[i4] = v.INFO.get(i4)

            df3 = df3.append(df2)

        else:

            count = 0
            for i2 in tempo_csq:
                count = count + 1 # CSQ subfield nb
                gene = i2.split('|')[3] # See protocole 109: gene taken in position 4 if no SYMBOL field
                impact = i2.split('|')[1]
                severity = i2.split('|')[2]

                df2=pd.DataFrame([[v.CHROM, v.POS, v.REF, v.ALT, ';'.join([i3[0]+"="+str(i3[1]) for i3 in v.INFO]), gene, severity, impact, aff, una, oddsratio, pvalue, -np.log10(pvalue), an, count]], columns = tsv_columns)

            if len(tsv_extra_fields_wo_csq) > 0: # add extra columns coming from tsv_extra_fields into the tsv file
                for i4 in tsv_extra_fields_wo_csq:
                    df2[i4] = v.INFO.get(i4)

                df3 = df3.append(df2)

    return df3


################################ End Functions


################################ Pre-ignition checking


# reserved words
# end reserved words
# argument primary checking
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
tsv_extra_fields = tsv_extra_fields.split(' ')
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
with open(ped, 'r') as pin:
    pin = csv.reader(pin, delimiter='\t')
    for row in pin:
        status[row[1]]=int(row[5])
# le dictionnaire ressemble a ca :
# {'IP00FN5': 2, 'IP00FLK': 2, 'IP00FLT': 1, 'IP00FM2': 2, 'IP00FMC': 1}

# le fichier vcf contenant les data des individus que l'on etudie
# attention doit etre bgzippé (.gz) et avoir un index tabix (.tbi)
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
tsv_columns=['CHROM','POS','REF','ALT', 'INFO', 'GENE_EXAMPLE','SEVERITY_EXAMPLE','IMPACT_EXAMPLE','AFF','UNAFF','OR','P_VALUE','NEG_LOG10_P_VALUE','PATIENT_NB'] #warning : can be replaced below
csq_subfield_name = []
csq_subfield_pos = []
tsv_extra_fields_wo_csq = []
if all([i0 == 'NULL' for i0 in tsv_extra_fields]) is False:
    tempo_log = [bool(re.search("^CSQ_.*$", i1)) for i1 in tsv_extra_fields] # is there any CSQ_ ?
    if any(i1 for i1 in tempo_log) is True:
        tsv_columns=['CHROM','POS','REF','ALT', 'INFO', 'GENE','SEVERITY','IMPACT','AFF','UNAFF','OR','P_VALUE','NEG_LOG10_P_VALUE','PATIENT_NB', 'CSQ_TRANSCRIPT_NB']
        tempo_pos = []
        for i2 in list(range(0, len(tsv_extra_fields))):
            if bool(re.search("^CSQ_.*$", tsv_extra_fields[i2])) is True:
                tempo = re.sub(pattern = "^CSQ_", repl = "", string = tsv_extra_fields[i2])
                csq_subfield_name.append(tempo)
            else:
                tempo_pos.append(i2)

        if all([i1 in vcf_csq_subfield_titles for i1 in csq_subfield_name]) is not True:
            sys.exit("Error in fisher_lod.py: some of the CSQ subfield of the tsv_extra_fields parameter (starting by CSQ_): \n"+" ".join(csq_subfield_name)+"\nare not in the vcf_csq_subfield_titles parameter: \n"+" ".join(vcf_csq_subfield_titles)+"\n")
        else:
            for i2 in list(range(0, len(vcf_csq_subfield_titles))):
                for i3 in csq_subfield_name:
                    if i3 == vcf_csq_subfield_titles[i2]:
                        csq_subfield_pos.append(i2)

        tsv_extra_fields_wo_csq = [tsv_extra_fields[i1] for i1 in tempo_pos]
        if all([i1 in vcf_info_field_titles for i1 in tsv_extra_fields_wo_csq]) is not True:
            sys.exit("Error in fisher_lod.py: not considering CSQ_, some of the tsv_extra_fields parameter values: \n"+" ".join(tsv_extra_fields)+"\nare not in the vcf_info_field_titles parameter: \n"+" ".join(vcf_info_field_titles)+"\n")

    elif all([i1 in vcf_info_field_titles for i1 in tsv_extra_fields]) is not True:
            sys.exit("Error in fisher_lod.py: some of the tsv_extra_fields parameter values: \n"+" ".join(tsv_extra_fields)+"\nare not in the vcf_info_field_titles parameter: \n"+" ".join(vcf_info_field_titles)+"\n")

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
    for v in vcf: # parse each line of vcf with the content of region in it
        if v.CHROM == region:
            tempo = fisher(v = v, tsv_columns = tsv_columns, tsv_extra_fields = tsv_extra_fields, csq_subfield_name = csq_subfield_name, csq_subfield_pos = csq_subfield_pos)
            df = df.append(tempo)

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





