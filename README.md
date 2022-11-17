[//]: # "#to make links in gitlab: example with racon https://github.com/isovic/racon"
[//]: # "tricks in markdown: https://openclassrooms.com/fr/courses/1304236-redigez-en-markdown"

| usage | dependencies |
| --- | --- |
| [![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/) | [![Dependencies: Nextflow Version](https://img.shields.io/badge/Nextflow-v21.04.2-blue?style=plastic)](https://github.com/nextflow-io/nextflow) |
| [![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses) | |

<br /><br />
## TABLE OF CONTENTS


   - [AIM](#aim)
   - [CONTENT](#content)
   - [HOW TO RUN](#how-to-run)
   - [OUTPUT](#output)
   - [VERSIONS](#versions)
   - [LICENCE](#licence)
   - [CITATION](#citation)
   - [CREDITS](#credits)
   - [ACKNOWLEDGEMENTS](#Acknowledgements)
   - [WHAT'S NEW IN](#what's-new-in)

<br /><br />
## AIM


Compute the Fisher exact test statistics (score) and p values from a vcf annotated file made of patient and control cases.
<br /><br />
Return a res.tsv file and a Miami plot.
<br /><br />
The tsv file can include fields and sub-fields of the vcf file in different columns. See the tsv_extra_fields parameter of the fisher_for_vcf.config file, as well as the OUTPUT section below.
<br /><br />
Return also a res.vcf file made from the res.tsv file, mimicing a VCF file, i.e., with the inital header of the .vcf and with the fisher results added in the INFO section. Warning: this is not a true VCF file as the results of the FORMAT field and corresponding patients data fields are not anymore present. Thus, this VCF file cannot be used by fisher_for_vcf as initial input.

<br /><br />
## CONTENT


**fisher_for_vcf.nf**: file that can be executed using a CLI (command line interface).

**fisher_for_vcf.config**: parameter settings for the fisher_for_vcf.nf file.

**dataset**: folder containing some datasets than can be used as examples.

| File | Description |
| --- | --- |
| **Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10000.vcf.gz** | First 10,000 lines of /pasteur/zeus/projets/p02/ghfc_wgs_zeus/WGS/Dyslexia/vcf/Dyslexia.gatk-vqsr.splitted.norm.vep.merged.vcf.gz |
| **Dyslexia.pedigree.txt** | Pedigree associated to Dyslexia.gatk-vqsr.splitted.norm.vep.merged.vcf.gz |
| **hg19_grch37p5_chr_size_cumul.txt** | Coordinates of the hg19_grch37p5 Human Genome for the Miami plot |

**example_of_results**: folder containing examples of result obtained with the dataset.
<br /><br />
&nbsp;&nbsp;&nbsp;&nbsp;See the OUTPUT section for the description of the folder and files.


<br /><br />
## HOW TO RUN


### From local using the committed version on gitlab

1) Create the scm file:

```bash
providers {
    pasteur {
        server = 'https://gitlab.pasteur.fr'
        platform = 'gitlab'
    }
}
```

And save it as 'scm' in the .nextflow folder. For instance in:
\\wsl$\Ubuntu-20.04\home\gael\.nextflow

Warning: ssh key must be set for gitlab, to be able to use this procedure (see protocol 44).


2) Mount a server if required:

```bash
DRIVE="C"
sudo mkdir /mnt/share
sudo mount -t drvfs $DRIVE: /mnt/share
```

Warning: if no mounting, it is possible that nextflow does nothing, or displays a message like
```
Launching `fisher_for_vcf.nf` [loving_morse] - revision: d5aabe528b
/mnt/share/Users
```


3) Then run the following command from here \\wsl$\Ubuntu-20.04\home\gael:

```bash
nextflow run -hub pasteur gmillot/fisher_for_vcf -r v1.0.0
```


4) If an error message appears, like:
```
WARN: Cannot read project manifest -- Cause: Remote resource not found: https://gitlab.pasteur.fr/api/v4/projects/gmillot%2Ffisher_for_vcf
```
Make the distant repo public


5) If an error message appears, like:

```
permission denied
```

See chmod in protocol 44.


### From local using local file

Like above but then run the following command from here \\wsl$\Ubuntu-20.04\home\gael:

```bash
nextflow run -c fisher_for_vcf.config fisher_for_vcf.nf
```

with -c to specify the name of the config file used.

If an error message appears, like:
```
Unknown error accessing project `gmillot/fisher_for_vcf` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/fisher_for_vcf
```
Purge using:
```
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
```


### From a cluster using a committed version on gitlab

Start with:

```bash
EXEC_PATH="/pasteur/zeus/projets/p01/BioIT/gmillot/08002_bourgeron" # where the bin folder of the fisher_for_vcf.nf script is located
export CONF_BEFORE=/opt/gensoft/exe # on maestro

export JAVA_CONF=java/13.0.2
export JAVA_CONF_AFTER=bin/java # on maestro
export SINGU_CONF=singularity
export SINGU_CONF_AFTER=bin/singularity # on maestro
export GIT_CONF=git/2.25.0
export GIT_CONF_AFTER=bin/git # on maestro

MODULES="${CONF_BEFORE}/${JAVA_CONF}/${JAVA_CONF_AFTER},${CONF_BEFORE}/${SINGU_CONF}/${SINGU_CONF_AFTER},${CONF_BEFORE}/${GIT_CONF}/${GIT_CONF_AFTER}"
cd ${EXEC_PATH}
chmod 755 ${EXEC_PATH}/bin/*.*
module load ${JAVA_CONF} ${SINGU_CONF} ${GIT_CONF}

```

Then run:

```bash
# distant fisher_for_vcf.nf file
HOME="${ZEUSHOME}/fisher_for_vcf/" ; trap '' SIGINT ; nextflow run --modules ${MODULES} -hub pasteur gmillot/fisher_for_vcf -r v1.0 -c $HOME/fisher_for_vcf.config ; HOME="/pasteur/appa/homes/gmillot/"  ; trap SIGINT

# local fisher_for_vcf.nf file ($HOME changed to allow the creation of .nextflow into /$ZEUSHOME/fisher_for_vcf/. See NFX_HOME in the nextflow soft script)
HOME="${ZEUSHOME}/fisher_for_vcf/" ; trap '' SIGINT ; nextflow run --modules ${MODULES} -c fisher_for_vcf.config fisher_for_vcf.nf ; HOME="/pasteur/appa/homes/gmillot/" ; trap SIGINT
```

If an error message appears, like:
```
Unknown error accessing project `gmillot/fisher_for_vcf` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/fisher_for_vcf
```
Purge using:
```
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
```


### To get the miamiplot only

Copy-paste this into a linux console

```
PWD=$(pwd)
# see the fisher_for_vcf.config file for info about the parameters
echo -e '
fisher="${PWD}/dataset/res_fisher.tsv.gz"
chr="${PWD}/dataset/hg19_grch37p5_chr_size_cumul.txt"
x_lim="chr1"
vgrid <- "TRUE"
top.y.column <- "NEG_LOG10_P_VALUE"
bottom.y.column <- "AF"
color.column <- "NULL"
dot.border.color <- "white"
y.lim1 <- "NULL"
y.lim2 <- "NULL"
reverse1 <- "FALSE"
reverse2 <- "TRUE"
y.threshold1 <- "1.2"
y.threshold2 <- "0.5"
y_log1="FALSE"
y_log2="TRUE"
cute="https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.4.0/cute_little_R_functions.R"

Rscript ${PWD}/bin/miami.R ${fisher} ${chr} "${x_lim}" "${vgrid}" "${top_y_column}" "${bottom_y_column}" "${color_column}" "${dot_border_color}" "${y_lim1}" "${y_lim2}" "${y_reverse1}" "${y_reverse2}" "${y_threshold1}" "${y_threshold2}" "${y_log1}" "${y_log2}" "${cute}" "miami_report.txt"
' | sudo docker run --workdir /tmp/ -i --mount "type=bind,src=${PWD},dst=/tmp/" --entrypoint bash gmillot/r_v4.1.2_extended_v2.1:gitlab_v8.8
```

The outputs files are in $PWD, i.e., where the code has been executed.

<br /><br />
## OUTPUT


**reports**: folder containing all the reports of the different processes including the *fisher_for_vcf.config* file used, *vcf_csq_subfield_titles.txt* (control file indicating the names of the CSQ subfields, as indicated in the header of the VCF file analyzed), *vcf_info_field_titles.txt* (control file indicating the names of the INFO fields, as indicated in the header of the VCF file analyzed)
<br /><br />
**Miami.png**: miami plot in the .png format.
<br /><br />
**res_fisher.tsv.gz**: table
<br /><br />
&nbsp;&nbsp;&nbsp;&nbsp;Rows:
<br />
<ul><li>Each row representing a different variant if the tsv_extra_fields parameter of the fisher_for_vcf.config file does not contain the CSQ field (VEP).
<br />
</li><li>Several lines per variant otherwise, depending on the number of subfields (comma separated) in the CSQ field (VEP) of the INFO field of the VCF file.
</li><br />
&nbsp;&nbsp;&nbsp;&nbsp;Columns:
<br />

| File | Description |
| --- | --- |
| **CHROM** | chromosome |
| **POS** | position in the chromosome (hg19)| 
| **REF** | nucleotide on the reference sequence (hg19) |
| **ALT** | alternative sequence |
| **INFO** | INFO field of the initial VCF file |
| **GENE** | name of the gene where the POS is (SYMBOL field of the CSQ field (VEP) of INFO field of the VCF) |
| **CONSEQUENCE** | ALT consequence (consequence field of the CSQ field (VEP) of INFO field of the VCF) |
| **IMPACT** | severity of the ALT consequence (IMPACT field of the CSQ field (VEP) of INFO field of the VCF) |
| **AFF** | Count of the number of 0=HOM_REF ; 1=HET ; 2=UNKNOWN ; 3=HOM_ALT in the affected cases. Example: {0:4, 1:2} for 4 cases HOM_REF and 2 cases HET |
| **UNAFF** | as in *aff* in the unaffected cases |
| **OR** | Odds ratio (n11/n12)/(n21/n22) = (n11\*n22)/(n12\*n21) with:<br /><ul><li>n11 = nHET_aff + nHOM_ALT_aff<br /></li><li>n12 = nHOM_REF_aff<br /></li><li>n21 = nHET_unaff + nHOM_ALT_unaff<br /></li><li>n22 = nHOM_REF_unaff<br /></li>OR > 1 meaning OR in favor of HET+HOM_ALT/aff versus HET+HOM_ALT/unaff |
| **P_VALUE** | p-value of the exact fisher test |
| **NEG_LOG10_P_VALUE** | -log10 of the p-value |
| **PATIENT_NB** | Number of AFF and UNAFF used for the fisher data |
| Optional colums | |
| **CSQ_TRANSCRIPT_NB** | number of fieds in the CSQ field (comma separated). Present only if "CSQ" is present in the tsv_extra_fields parameter |
| ***\<NAME\>*** | name of the fields of INFO field of the vcf or subfield of CSQ, as indicated in the tsv_extra_fields parameter |

<br /><br />
**res_fisher.vcf.gz**: file made from the res.tsv file, mimicing a VCF file, i.e., with the inital header of the .vcf and with the fisher results added in the INFO section. Warning: this is not a true VCF file as the results of the FORMAT field and corresponding patients data fields are not anymore present. Thus, this VCF file cannot be used by fisher_for_vcf as initial input.
.

<br /><br />
## VERSIONS


The different releases are tagged [here](https://gitlab.pasteur.fr/gmillot/fisher_for_vcf/-/tags)

<br /><br />
## LICENCE


This package of scripts can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchandability or fitness for a particular purpose.
See the GNU General Public License for more details at https://www.gnu.org/licenses.

<br /><br />
## CITATION


Not yet published


<br /><br />
## CREDITS


[Freddy Cliquet](https://gitlab.pasteur.fr/fcliquet), GHFC, Institut Pasteur, Paris, France

[Gael A. Millot](https://gitlab.pasteur.fr/gmillot), Hub-CBD, Institut Pasteur, Paris, France

<br /><br />
## ACKNOWLEDGEMENTS


The mentioned softwares and packages developers & maintainers

Gitlab developers

<br /><br />
## WHAT'S NEW IN


### v2.4

miamiplot: new features added and bug removed


### v2.3

miamiplot: log10 scale and alignment of graphics


### v2.2

bug fixed in miamiplot
Now possible to get the miamiplot only (see the README file)


### v2.1

Miamiplot improved for right display of the color legend


### v2.0

Code debugged because it was overwritting the ipunt file if named res.vcf.gz


### v1.10

Parall debugged<br />
Check in python part improved


### v1.9

Zipping added


### v1.8

Code improved


### v1.7

Code improved


### v1.6

miami plot updated

### v1.5

fields of the VCF (like AC or subfield of CSQ like SIFT) can be added into the tsv file


### v1.4

now the script deals with .gz and non .gz vcf


### v1.3

INFO column added in the .tsv file


### v1.2

Now the Miami plot can be zoomed on a particular region thanks to the region and x_lim parameters. Remain to get the correct INFO column in the .py code


### v1.1

Miami plot improvement


### v1.0

everything




