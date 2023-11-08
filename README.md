[//]: # "#to make links in gitlab: example with racon https://github.com/isovic/racon"
[//]: # "tricks in markdown: https://openclassrooms.com/fr/courses/1304236-redigez-en-markdown"

| usage | dependencies |
| --- | --- |
| [![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/) | [![Dependencies: Nextflow Version](https://img.shields.io/badge/Nextflow-v21.04.2-blue?style=plastic)](https://github.com/nextflow-io/nextflow) |
| [![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses) | |

<br /><br />
## TABLE OF CONTENTS


   - [AIM](#aim)
   - [WARNING](#warning)
   - [CONTENT](#content)
   - [INPUT](#input)
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


- Compute the Fisher exact test statistics (score) and p values from a vcf annotated file made of patient and control cases.
- Return a res.tsv file and a Miami plot.
- The tsv file can include fields and sub-fields of the vcf file in different columns. See the tsv_extra_fields parameter of the fisher_for_vcf.config file, as well as the OUTPUT section below.
- Return also a res.vcf file made from the res.tsv file, mimicing a VCF file, i.e., with the inital header of the .vcf and with the fisher results added in the INFO section.

<br /><br />
## WARNINGS

The returned res.vcf is not a true VCF file as the results of the FORMAT field and corresponding patients data fields are not anymore present. Thus, this VCF file cannot be used by fisher_for_vcf as initial input.

<br /><br />
## CONTENT

| fisher_for_vcf folder | Description |
| --- | --- |
| **main.nf** | File that can be executed using a linux terminal, a MacOS terminal or Windows 10 WSL2. |
| **nextflow.config** | Parameter settings for the *main.nf* file. Users have to open this file, set the desired settings and save these modifications before execution. |
| **bin folder** | Contains files required by the *main.nf* file. |


<br /><br />
## INPUT

| Required files |
| --- |
| Variant Calling Format (VCF) file (zipped or not). |
| Pedigree file. |
| Human chromo infos file. |


<br /><br />
The dataset used in the *nextflow.config* file, as example, is available at https://zenodo.org/records/10084500/files/fisher_for_vcf.zip


| Dataset folder | Description |
| --- | --- |
| **Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10000.vcf.gz** | VCF file. |
| **Dyslexia.pedigree.txt** | Pedigree file. |
| **hg19_grch37p5_chr_size_cumul.txt** | Coordinates of the hg19_grch37p5 Human Genome for the Miami plot |


<br /><br />
## HOW TO RUN

### 1. Prerequisite

Installation of:<br />
[nextflow DSL1](https://github.com/nextflow-io/nextflow)<br />
[Graphviz](https://www.graphviz.org/download/), `sudo apt install graphviz` for Linux ubuntu<br />
[Apptainer](https://github.com/apptainer/apptainer)<br />

<br /><br />
### 2. Local running (personal computer)


#### 2.1. *main.nf* file in the personal computer

- Mount a server if required:

<pre>
DRIVE="Z" # change the letter to fit the correct drive
sudo mkdir /mnt/share
sudo mount -t drvfs $DRIVE: /mnt/share
</pre>

Warning: if no mounting, it is possible that nextflow does nothing, or displays a message like:
<pre>
Launching `main.nf` [loving_morse] - revision: d5aabe528b
/mnt/share/Users
</pre>

- Run the following command from where the *main.nf* and *nextflow.config* files are (example: \\wsl$\Ubuntu-20.04\home\gael):

<pre>
nextflow run main.nf -c nextflow.config
</pre>

with -c to specify the name of the config file used.

<br /><br />
#### 2.3. *main.nf* file in the public gitlab repository

Run the following command from where you want the results:

<pre>
nextflow run -hub pasteur gmillot/fisher_for_vcf -r v1.0.0
</pre>

<br /><br />
### 3. Distant running (example with the Pasteur cluster)

#### 3.1. Pre-execution

Copy-paste this after having modified the EXEC_PATH variable:

<pre>
EXEC_PATH="/pasteur/zeus/projets/p01/BioIT/gmillot/fisher_for_vcf" # where the bin folder of the main.nf script is located
export CONF_BEFORE=/opt/gensoft/exe # on maestro

export JAVA_CONF=java/13.0.2
export JAVA_CONF_AFTER=bin/java # on maestro
export APP_CONF=apptainer/1.2.3
export APP_CONF_AFTER=bin/apptainer # on maestro
export GIT_CONF=git/2.39.1
export GIT_CONF_AFTER=bin/git # on maestro
export GRAPHVIZ_CONF=graphviz/2.42.3
export GRAPHVIZ_CONF_AFTER=bin/graphviz # on maestro

MODULES="${CONF_BEFORE}/${JAVA_CONF}/${JAVA_CONF_AFTER},${CONF_BEFORE}/${APP_CONF}/${APP_CONF_AFTER},${CONF_BEFORE}/${GIT_CONF}/${GIT_CONF_AFTER}/${GRAPHVIZ_CONF}/${GRAPHVIZ_CONF_AFTER}"
cd ${EXEC_PATH}
chmod 755 ${EXEC_PATH}/bin/*.* # not required if no bin folder
module load ${JAVA_CONF} ${APP_CONF} ${GIT_CONF} ${GRAPHVIZ_CONF}
</pre>

<br /><br />
#### 3.2. *main.nf* file in a cluster folder

Modify the second line of the code below, and run from where the *main.nf* and *nextflow.config* files are (which has been set thanks to the EXEC_PATH variable above):

<pre>
HOME_INI=$HOME
HOME="${ZEUSHOME}/fisher_for_vcf/" # $HOME changed to allow the creation of .nextflow into /$ZEUSHOME/fisher_for_vcf/, for instance. See NFX_HOME in the nextflow software script
trap '' SIGINT
nextflow run --modules ${MODULES} main.nf -c nextflow.config
HOME=$HOME_INI
trap SIGINT
</pre>

<br /><br />
#### 3.3. *main.nf* file in the public gitlab repository

Modify the first and third lines of the code below, and run (results will be where the EXEC_PATH variable has been set above):

<pre>
VERSION="v1.0"
HOME_INI=$HOME
HOME="${ZEUSHOME}/fisher_for_vcf/" # $HOME changed to allow the creation of .nextflow into /$ZEUSHOME/fisher_for_vcf/, for instance. See NFX_HOME in the nextflow software script
trap '' SIGINT
nextflow run --modules ${MODULES} -hub pasteur gmillot/fisher_for_vcf -r $VERSION -c $HOME/nextflow.config
HOME=$HOME_INI
trap SIGINT
</pre>

<br /><br />
### 4. Error messages and solutions

#### Message 1
```
Unknown error accessing project `gmillot/fisher_for_vcf` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/fisher_for_vcf
```

Purge using:
<pre>
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
</pre>

#### Message 2
```
WARN: Cannot read project manifest -- Cause: Remote resource not found: https://gitlab.pasteur.fr/api/v4/projects/gmillot%2Ffisher_for_vcf
```

Contact Gael Millot (distant repository is not public).

#### Message 3

```
permission denied
```

Use chmod to change the user rights. Example linked to files in the bin folder: 
```
chmod 755 bin/*.*
```


<br /><br />
## OUTPUT

An example of results obtained with the dataset is present at this address: https://zenodo.org/records/10084500/files/fisher_for_vcf.zip
<br /><br />


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

[Gael A. Millot](https://gitlab.pasteur.fr/gmillot), Hub, Institut Pasteur, Paris, France

<br /><br />
## ACKNOWLEDGEMENTS


The developers & maintainers of the mentioned softwares and packages, including:

- [Slivar](https://github.com/brentp/slivar)
- [Nextflow](https://www.nextflow.io/)
- [Apptainer](https://apptainer.org/)
- [Docker](https://www.docker.com/)
- [Gitlab](https://about.gitlab.com/)
- [Bash](https://www.gnu.org/software/bash/)
- [Ubuntu](https://ubuntu.com/)


<br /><br />
## WHAT'S NEW IN

### v2.5

- .config improved
- README improved. Dataset and results are in zenodo

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




