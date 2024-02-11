#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#########################################################################
##                                                                     ##
##     add_fisher.py                                                   ##
##                                                                     ##
##     Freddy Cliquet                                                  ##
##     GHFC                                                            ##
##     Neuroscience Department                                         ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################



################################ Aim


################################ End Aim


################################ Introduction

# https://brentp.github.io/cyvcf2/docstrings.html

################################ End Introduction


################################ Acknowlegments


################################ End Acknowlegments


################################ Initialization




import csv
import sys
import re
import os.path
# sys.path.append('/home/gmillot/anaconda3/lib/python3.11/site-packages') #additional package path added into PYTHONPATH
from cyvcf2 import VCF
import scipy.stats as stats
import pandas as pd
import numpy as np
import random # for shuffle() function
import warnings # for warnings.catch_warnings()
from itertools import compress


################################ End Initialization


################################ Parameters that need to be set by the user



################################ End Parameters that need to be set by the user


################################ Config import


vcf_path = sys.argv[1]  # 1st argument: pedigree file name, 1 ID par ligne, sys.argv takes arguments from the bash line command when running a .py script
ped_path = sys.argv[2]
vcf_info_field_titles_path = sys.argv[3]
tsv_extra_fields = sys.argv[4]
vcf_csq_subfield_titles_path = sys.argv[5]
filter_indiv_DP = sys.argv[6]
filter_indiv_GQ = sys.argv[7]
model = sys.argv[8]


################################ End Config import

################################ Test


# vcf_path="/mnt/c/Users/gmillot/Documents/Git_projects/fisher_for_vcf/dataset/other/res.vcf.gz"
# ped_path="/mnt/c/Users/gmillot/Documents/Git_projects/fisher_for_vcf/dataset/other/pedigree.txt"
# vcf_info_field_titles_path = "/mnt/c/Users/gmillot/Documents/Git_projects/fisher_for_vcf/dataset/other/vcf_info_field_titles.txt"
# tsv_extra_fields =  "AC AF CSQ_SIFT CSQ_PolyPhen"
# vcf_csq_subfield_titles_path = "/mnt/c/Users/gmillot/Documents/Git_projects/fisher_for_vcf/dataset/other/vcf_csq_subfield_titles.txt"
# filter_indiv_DP = "30"
# filter_indiv_GQ = "10"
# model = 'carrier|strict_heterozygous'


################################ End Test

################################ Recording of the initial parameters




################################ End Recording of the initial parameters


################################ Functions


def fisher(v, status, model, tsv_columns, tsv_extra_fields_wo_csq, csq_subfield_name, csq_subfield_pos, filter_indiv_DP, filter_indiv_GQ):
    '''
    AIM
        compute fisher for a line of a vcf
    WARNINGS
        Consider the five first fields as: 'CHROM', 'POS', 'REF', 'ALT', 'INFO'
        Consider the three next fields as: 'GENE', 'IMPACT', 'CONSEQUENCE'
        Consider that gene, impact and consequence are in position 3, 2, and 1 in CSQ subfield of INFO
    ARGUMENTS
        v: a single line of a vcf object from VCF()
        status: dictionnary of the "aff" = 2 or "unaff" = 1 status of each indiv
        model: the model string in the config 
        tsv_columns: column names of the final tsv file
        tsv_extra_fields_wo_csq: subfields of INFO field of the vcf to add as column in the tsv (excluding potential CSQ subfields like Polyphen)
        csq_subfield_name: subfields of the CSQ subfield of the INFO field of the vcf to add as column in the tsv
        csq_subfield_pos: positions of csq_subfield_name in all the subfields of CSQ, indicated by vcf_csq_subfield_titles
    RETURN
        a data frame with a single row
    REQUIRED PACKAGES
        None
    EXAMPLE
        for v in vcf:
            fisher(v = v, tsv_columns = tsv_columns, tsv_extra_fields = tsv_extra_fields)
    DEBUGGING
        # use the container: sudo docker run -ti --entrypoint bash -v /mnt/c/Users/gael/Documents/Git_projects/fisher_for_vcf/dataset:/tmp gmillot/python_v3.9.10_extended_v3.1:gitlab_v8.7)
vcf = VCF("example.vcf.gz")
count = 1 
for v in vcf:
    if count == 1: break

# use dir(v) for the methods of v and str(v) for the string in v and len(str(v).split('\t')) for the fields nb
status = {'P1': 1, 'P2': 2, 'P3': 1, 'P4': 2, 'P5': 2, 'P6': 1, 'P7': 2, 'P8': 1, 'P9': 2, 'P10': 2, 'P11': 2}
model = 'carrier|strict_heterozygous'
tsv_columns = ['CHROM', 'POS', 'REF', 'ALT', 'INFO', 'GENE', 'IMPACT', 'CONSEQUENCE', 'PATIENT_NB', 'AFF', 'UNAFF', 'N_HOM_REF_AFF', 'N_HET_AFF', 'N_HOM_ALT_AFF', 'N_HOM_REF_UNAFF', 'N_HET_UNAFF', 'N_HOM_ALT_UNAFF', 'N11_CARRIER_MODEL', 'N12_CARRIER_MODEL', 'N21_CARRIER_MODEL', 'N22_CARRIER_MODEL', 'OR_CARRIER_MODEL', 'P_VALUE_CARRIER_MODEL', 'NEG_LOG10_P_VALUE_CARRIER_MODEL', 'N11_HETERO_MODEL', 'N12_HETERO_MODEL', 'N21_HETERO_MODEL', 'N22_HETERO_MODEL', 'OR_HETERO_MODEL', 'P_VALUE_HETERO_MODEL', 'NEG_LOG10_P_VALUE_HETERO_MODEL', 'CSQ_TRANSCRIPT_NB']
tsv_extra_fields_wo_csq = ['AC', 'AF']
csq_subfield_name = ['SIFT', 'PolyPhen']
csq_subfield_pos = [25, 26]
filter_indiv_DP = 30.0
filter_indiv_GQ = 10.0
    '''
    
    # je fais 2 dictionnaire pour stocker mes compte d'atteint/non atteint porteurs/non porteurs
    # dans chaque dictionaire j'associerai a un genotype (clef) un nombre d'individu.
    aff=dict() #dictionary: return of this is {1:0, 2:0, } with key:value each element of the dict (key -> integers by default)
    una=dict()
    
    column_names = tsv_columns + csq_subfield_name + tsv_extra_fields_wo_csq # column names of the future data frame # only tsv_columns if csq_subfield_name = [] and  tsv_extra_fields_wo_csq = [], i.e., tsv_extra_fields = 'NULL'
    
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
            if dp not in ['.', ''] and int(dp)>=filter_indiv_DP and gq not in ['.', ''] and float(gq)>=filter_indiv_GQ:
                # je met a jour le compte d'individus dans mes dictionnaire aff et una en fonction du phenotype de l'individu courant iid
                if status[iid]==2: # 2 == status "aff" (result from ped)
                    aff[gt]=aff.get(gt,0)+1 # means if the key gt exists in the aff dict, return the value of aff[gt] + 1, otherwise return 0 + 1
                if status[iid]==1: # 1 == status "unaff" (result from ped)
                    una[gt]=una.get(gt,0)+1
                an+=1
    
    # une fois que l'on a lu les information pour tous les individus, nous calculons le Fisher
    # ici c'est porteur (gt 1 ou 3) versus non porteur (gt 0) pour les atteints (aff) versus les non atteint (una)
    n_hom_ref_aff = aff.get(0,0)
    n_het_aff = aff.get(1,0)
    n_hom_alt_aff = aff.get(3,0)
    n_hom_ref_unaff = una.get(0,0)
    n_het_unaff = una.get(1,0)
    n_hom_alt_unaff = una.get(3,0)
    n_aff = n_hom_ref_aff + n_het_aff + n_hom_alt_aff
    n_unaff = n_hom_ref_unaff + n_het_unaff + n_hom_alt_unaff
    
    list1 = [v.CHROM, v.POS, v.REF, ''.join(v.ALT), ';'.join([i1[0]+"="+str(i1[1]) for i1 in v.INFO])] # length 5
    list3 = [an, aff, una, n_hom_ref_aff, n_het_aff, n_hom_alt_aff, n_hom_ref_unaff, n_het_unaff, n_hom_alt_unaff] # length 9 # between list1 and list3, the 3 positions 5, 6, 7 (start 0) will be filled later. They are for gene, consequence and impact, according to the order in column_names (tsv_columns)
    if bool(re.search("carrier", model)) is True: 
        n11_carrier_model = n_het_aff+n_hom_alt_aff # nb of hetero and HOMO ALT in aff
        n12_carrier_model = n_hom_ref_aff # nb of HOMO REF in aff
        n21_carrier_model = n_het_unaff+n_hom_alt_unaff # nb of hetero and HOMO ALT in unaff 
        n22_carrier_model = n_hom_ref_unaff # nb of HOMO REF in unaff
        oddsratio_carrier_model, pvalue_carrier_model = stats.fisher_exact([[n11_carrier_model, n12_carrier_model],[n21_carrier_model, n22_carrier_model]], alternative='two-sided')
        list3.extend([n11_carrier_model, n12_carrier_model, n21_carrier_model, n22_carrier_model, oddsratio_carrier_model, pvalue_carrier_model, -np.log10(pvalue_carrier_model)]) # length 7
    
    if bool(re.search("strict_heterozygous", model)) is True: 
        n11_hetero_model = n_het_aff # nb of hetero (expected geno) in aff
        n12_hetero_model = n_hom_ref_aff+n_hom_alt_aff # nb of HOMO REF or ALT in aff
        n21_hetero_model = n_het_unaff # nb of hetero (expected geno) in unaff 
        n22_hetero_model = n_hom_ref_unaff+n_hom_alt_unaff # nb of HOMO REF or ALT in unaff
        oddsratio_hetero_model, pvalue_hetero_model = stats.fisher_exact([[n11_hetero_model, n12_hetero_model],[n21_hetero_model, n22_hetero_model]], alternative='two-sided')
        list3.extend([n11_hetero_model, n12_hetero_model, n21_hetero_model, n22_hetero_model, oddsratio_hetero_model, pvalue_hetero_model, -np.log10(pvalue_hetero_model)]) # length 7
    
    if bool(re.search("recessive", model)) is True: 
        n11_recessive_model = n_hom_alt_aff # nb of HOMO ALT in aff
        n12_recessive_model = n_hom_ref_aff+n_het_aff # nb of HOMO REF and hetero in aff
        n21_recessive_model = n_hom_alt_unaff # nb of HOMO ALT in unaff
        n22_recessive_model = n_hom_ref_unaff+n_het_unaff # nb of HOMO REF and hetero in unaff
        oddsratio_recessive_model, pvalue_recessive_model = stats.fisher_exact([[n11_recessive_model, n12_recessive_model],[n21_recessive_model, n22_recessive_model]], alternative='two-sided')
        list3.extend([n11_recessive_model, n12_recessive_model, n21_recessive_model, n22_recessive_model, oddsratio_recessive_model, pvalue_recessive_model, -np.log10(pvalue_recessive_model)]) # length 7
    
    # extra fields required in the tsv in INFO field but not in CSQ of INFO
    list5 = []
    if len(tsv_extra_fields_wo_csq) > 0: # add extra columns coming from tsv_extra_fields into the tsv file
        for i2 in tsv_extra_fields_wo_csq:
            list5.extend([v.INFO.get(i2)])
    
    # je traite tous les champ CSQ (annotation VEP)
    if v.INFO.get('CSQ') is not None: # https://www.w3schools.com/python/ref_dictionary_get.asp
        # Warning: CSQ can have several fields is the same variant fall into different annotations (2 genes for instance). Thus, tempo_csq is used
        # filling a one row data frame with or without adding
        tempo_csq = v.INFO.get('CSQ').split(',') # number of fields in CSQ (comma sep), i.e., nb of rows in the future data frame. Comma separated means that several genes for the same position
        part1 = [list(list1) for _ in range(len(tempo_csq))] # Repeat the list len(tempo_csq) times and convert to a 2D list with len(tempo_csq) that will be len(tempo_csq) rows in the future df. Do not use array because a single type inside
        list3.extend([len(tempo_csq)]) # add the CSQ_TRANSCRIPT_NB field in the tsv file
        part3 = [list(list3) for _ in range(len(tempo_csq))] # Repeat the list len(tempo_csq) times and convert to a 2D list with len(tempo_csq) that will be len(tempo_csq) rows in the future df. Do not use array because a single type inside
        if len(tsv_extra_fields_wo_csq) > 0: # add extra columns coming from tsv_extra_fields into the tsv file
            part5 = [list(list5) for _ in range(len(tempo_csq))]
        
        part2 = []
        for i2 in tempo_csq: # len(tempo_csq) indicates the number of row that will be created for each field of CSQ (comma separated)
            part2.append([i2.split('|')[3], i2.split('|')[2], i2.split('|')[1]]) # take the fields gene, impact and consequence, which are in position 3, 2, and 1, and make a 2D list # See protocole 109: gene taken in position 4 if no SYMBOL field
        
        if len(csq_subfield_name) > 0: # meaning that the user wanted additional fields
            part4 = []
            for i3 in list(range(0, len(tempo_csq))): # number of subfields of CSQ wanted, i.e., number of columns in the future df
                tempo = [] # tempo list
                for i4 in list(range(0, len(csq_subfield_pos))):
                    tempo.append(tempo_csq[i3].split('|')[csq_subfield_pos[i4]])
                
                part4.append(tempo)
    else:
        part1 = [list(list1) for _ in range(1)]
        list3.extend(['']) # add the empty CSQ_TRANSCRIPT_NB field in the tsv file
        part3 = [list(list3) for _ in range(1)]
        part2 = [['', '', '']] # empty gene, impact and consequence because no CSQ
        if len(tsv_extra_fields_wo_csq) > 0:
            part5 = [list(list5) for _ in range(1)]
        
        if len(csq_subfield_name) > 0: # meaning that the user wanted additional fields
            part4 = []
            for i3 in list(range(0, len(csq_subfield_name))):
                part4.extend([''])
            
            part4 = [part4]
    
    # Convert part1, part2, part3 to DataFrame
    dfpart1 = pd.DataFrame(np.array(part1))
    dfpart2 = pd.DataFrame(np.array(part2))
    dfpart3 = pd.DataFrame(np.array(part3))
    
    df = pd.concat([dfpart1, dfpart2, dfpart3], axis=1, join="inner")
    
    if 'part4' in locals(): # or if len(csq_subfield_name) > 0:
        # If part4 exists, concatenate it to result
        dfpart4 = pd.DataFrame(np.array(part4))
        df = pd.concat([df, dfpart4], axis=1, join="inner")
    
    if 'part5' in locals(): # or if len(csq_subfield_name) > 0:
        # If part4 exists, concatenate it to result
        dfpart5 = pd.DataFrame(np.array(part5))
        df = pd.concat([df, dfpart5], axis=1, join="inner")
    
    df = df[~df.astype(str).duplicated()] # remove duplicated rows. astype(str) convert df into string, because duplicated() cannot work if dict columns present. Tild is for not (invert true and false)
    df.columns = column_names
    return df

################################ End Functions


################################ Pre-ignition checking


# reserved words
# end reserved words
# argument primary checking
if isinstance(vcf_path, str) is not True:
    sys.exit("\n\n========\n\nError in add_fisher.py: the vcf_path must be a single character string: \n"+"".join(vcf_path)+"\n\n========\n\n")

if isinstance(ped_path, str) is not True:
    sys.exit("\n\n========\n\nError in add_fisher.py: the ped_path must be a single character string: \n"+"".join(ped_path)+"\n\n========\n\n")

if isinstance(vcf_info_field_titles_path, str) is not True:
    sys.exit("\n\n========\n\nError in add_fisher.py: the vcf_info_field_titles_path must be a single character string: \n"+"".join(vcf_info_field_titles_path)+"\n\n========\n\n")

if isinstance(tsv_extra_fields, str) is not True:
    sys.exit("\n\n========\n\nError in add_fisher.py: the tsv_extra_fields must be a single character string: \n"+"".join(tsv_extra_fields)+"\n\n========\n\n")

if isinstance(filter_indiv_DP, str) is not True:
    sys.exit("\n\n========\n\nError in add_fisher.py: the filter_indiv_DP must be a single character string made of a float number: \n"+"".join(str(filter_indiv_DP))+"\n\n========\n\n")
elif not all(char.isdigit() or char == '.' for char in filter_indiv_DP): # if not isinstance(filter_indiv_DP, (float, int)) :
    sys.exit("\n\n========\n\nError in add_fisher.py: the filter_indiv_DP must be a single character string made of a float number: \n"+"".join(filter_indiv_DP)+"\nType: "+str(type(filter_indiv_DP))+"\n\n========\n\n")
else:
    filter_indiv_DP = float(filter_indiv_DP)

if isinstance(filter_indiv_GQ, str) is not True:
    sys.exit("\n\n========\n\nError in add_fisher.py: the filter_indiv_GQ must be a single character string made of a float number: \n"+"".join(str(filter_indiv_GQ))+"\n\n========\n\n")
elif not all(char.isdigit() or char == '.' for char in filter_indiv_GQ): # if not isinstance(filter_indiv_GQ, (float, int)) :
    sys.exit("\n\n========\n\nError in add_fisher.py: the filter_indiv_GQ must be a single character string made of a float number: \n"+"".join(filter_indiv_GQ)+"\nType: "+str(type(filter_indiv_GQ))+"\n\n========\n\n")
else:
    filter_indiv_GQ = float(filter_indiv_GQ)

if isinstance(vcf_csq_subfield_titles_path, str) is not True:
    sys.exit("\n\n========\n\nError in add_fisher.py: the vcf_csq_subfield_titles_path must be a single character string: \n"+"".join(vcf_csq_subfield_titles_path)+"\n\n========\n\n")

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
    sys.exit("\n\n========\n\nError in add_fisher.py: the vcf_path file does not exist: \n"+"".join(vcf_path)+"\n\n========\n\n")

if os.path.exists(ped_path) is not True:
    sys.exit("\n\n========\n\nError in add_fisher.py: the ped_path file does not exist: \n"+"".join(ped_path)+"\n\n========\n\n")

if os.path.exists(vcf_info_field_titles_path) is not True:
    sys.exit("\n\n========\n\nError in add_fisher.py: the vcf_info_field_titles_path file does not exist: \n"+"".join(vcf_info_field_titles_path)+"\n\n========\n\n")

if all([i0 == 'NULL' for i0 in tsv_extra_fields]) is False:
    tsv_extra_fields = tsv_extra_fields.split(' ') # list

if os.path.exists(vcf_csq_subfield_titles_path) is not True:
    sys.exit("\n\n========\n\nError in add_fisher.py: the vcf_csq_subfield_titles_path file does not exist: \n"+"".join(vcf_csq_subfield_titles_path)+"\n\n========\n\n")

# end other checkings
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking


################################ End pre-ignition checking


################################ Main code


################ Ignition


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
# the INFO field of the VCF is a complex object

if sorted(vcf.samples) != sorted(list(status.keys())):
    sys.exit("\n\n========\n\nError in add_fisher.py: names of indiv must be the same in the pedigree file and in the samples of the VCF file (last sample fields).\nHere it is:\npedigree\n"+" ".join(sorted(list(status.keys())))+"\nvcf\n"+" ".join(sorted(vcf.samples))+"\n\n========\n\n")

with open(vcf_info_field_titles_path, 'r') as f:
    vcf_info_field_titles = f.readlines()[0].split(' ')

with open(vcf_csq_subfield_titles_path, 'r') as f:
    vcf_csq_subfield_titles = f.readlines()[0].split(' ')

################ end Data import


############ modifications of imported tables

# checking that everything is fine with tsv_extra_fields, and info recovering
tsv_columns = ['CHROM', 'POS', 'REF', 'ALT', 'INFO', 'GENE', 'IMPACT', 'CONSEQUENCE', 'PATIENT_NB', 'AFF', 'UNAFF', 'N_HOM_REF_AFF', 'N_HET_AFF', 'N_HOM_ALT_AFF', 'N_HOM_REF_UNAFF', 'N_HET_UNAFF', 'N_HOM_ALT_UNAFF'] #warning : can be replaced below
if bool(re.search("carrier", model)) is True: 
    tsv_columns = tsv_columns + ['N11_CARRIER_MODEL', 'N12_CARRIER_MODEL', 'N21_CARRIER_MODEL', 'N22_CARRIER_MODEL', 'OR_CARRIER_MODEL', 'P_VALUE_CARRIER_MODEL', 'NEG_LOG10_P_VALUE_CARRIER_MODEL']

if bool(re.search("strict_heterozygous", model)) is True: 
    tsv_columns = tsv_columns + ['N11_HETERO_MODEL', 'N12_HETERO_MODEL', 'N21_HETERO_MODEL', 'N22_HETERO_MODEL', 'OR_HETERO_MODEL', 'P_VALUE_HETERO_MODEL', 'NEG_LOG10_P_VALUE_HETERO_MODEL'] 

if bool(re.search("recessive", model)) is True: 
    tsv_columns = tsv_columns + ['N11_RECESS_MODEL', 'N12_RECESS_MODEL', 'N21_RECESS_MODEL', 'N22_RECESS_MODEL', 'OR_RECESS_MODEL', 'P_VALUE_RECESS_MODEL', 'NEG_LOG10_P_VALUE_RECESS_MODEL']

tsv_columns = tsv_columns + ['CSQ_TRANSCRIPT_NB']

csq_subfield_name = []
csq_subfield_pos = []
tsv_extra_fields_wo_csq = []
if all([i0 == 'NULL' for i0 in tsv_extra_fields]) is False:
    tempo_log = [bool(re.search("^CSQ_.*$", i1)) for i1 in tsv_extra_fields] # is there any CSQ_ ?
    if any(i1 for i1 in tempo_log) is True:
        tempo_pos = []
        for i2 in list(range(0, len(tsv_extra_fields))):
            if bool(re.search("^CSQ_.*$", tsv_extra_fields[i2])) is True:
                tempo = re.sub(pattern = "^CSQ_", repl = "", string = tsv_extra_fields[i2])
                csq_subfield_name.append(tempo)
            else:
                tempo_pos.append(i2)
        
        if all([i1 in vcf_csq_subfield_titles for i1 in csq_subfield_name]) is not True:
            sys.exit("\n\n========\n\nError in add_fisher.py: some of the CSQ subfield of the tsv_extra_fields parameter (starting by CSQ_): \n"+"".join(csq_subfield_name)+"\nare not in the vcf_csq_subfield_titles parameter: \n"+"".join(vcf_csq_subfield_titles)+"\nThe wrong fields are:\n"+"".join(compress(vcf_csq_subfield_titles, [ not i1 in vcf_csq_subfield_titles for i1 in csq_subfield_name]))+"\n\n========\n\n")
        else:
            for i2 in csq_subfield_name:
                for i3 in list(range(0, len(vcf_csq_subfield_titles))):
                    if i2 == vcf_csq_subfield_titles[i3]:
                        csq_subfield_pos.append(i3)
        
        tsv_extra_fields_wo_csq = [tsv_extra_fields[i1] for i1 in tempo_pos]
        if all([i1 in vcf_info_field_titles for i1 in tsv_extra_fields_wo_csq]) is not True:
            sys.exit("\n\n========\n\nError in add_fisher.py: not considering CSQ_, some of the tsv_extra_fields parameter values: \n"+"".join(tsv_extra_fields_wo_csq)+"\nare not in the vcf_info_field_titles parameter: \n"+"".join(vcf_info_field_titles)+"\nThe wrong fields are:\n"+"".join(compress(vcf_info_field_titles, [ not i1 in vcf_info_field_titles for i1 in tsv_extra_fields_wo_csq]))+"\n\n========\n\n")
    elif all([i1 in vcf_info_field_titles for i1 in tsv_extra_fields]) is not True:
            sys.exit("\n\n========\n\nError in add_fisher.py: some of the tsv_extra_fields parameter values: \n"+"".join(tsv_extra_fields)+"\nare not in the vcf_info_field_titles parameter: \n"+"".join(vcf_info_field_titles)+"\nThe wrong fields are:\n"+"".join(compress(vcf_info_field_titles, [ not i1 in vcf_info_field_titles for i1 in tsv_extra_fields]))+"\n\n========\n\n")
    else:
        tsv_extra_fields_wo_csq = tsv_extra_fields

df = pd.DataFrame(columns = tsv_columns + csq_subfield_name + tsv_extra_fields_wo_csq) #just tsv_columns if tsv_extra_fields is 'NULL'
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
        tempo = fisher(v = v, status = status, model = model, tsv_columns = tsv_columns, tsv_extra_fields_wo_csq = tsv_extra_fields_wo_csq, csq_subfield_name = csq_subfield_name, csq_subfield_pos = csq_subfield_pos, filter_indiv_DP = filter_indiv_DP, filter_indiv_GQ = filter_indiv_GQ)
        df = pd.concat([df, tempo]) # .append deprecated in panda

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





