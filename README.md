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
Return a tsv file and a Miami plot.

<br /><br />
## CONTENT

**fisher_for_vcf.nf** file that can be executed using a CLI (command line interface)

**fisher_for_vcf.config** parameter settings for the fisher_for_vcf.nf file

**dataset** folder containing some datasets than can be used as examples

Dyslexia.gatk-vqsr.splitted.norm.vep.merged_first_10000.vcf.gz: first 10,000 lines of /pasteur/zeus/projets/p02/ghfc_wgs_zeus/WGS/Dyslexia/vcf/Dyslexia.gatk-vqsr.splitted.norm.vep.merged.vcf.gz

Dyslexia.pedigree.txt: pedigree associated to Dyslexia.gatk-vqsr.splitted.norm.vep.merged.vcf.gz

hg19_grch37p5_chr_size_cumul.txt: coordinates of the hg19_grch37p5 Human Genome for the Miami plot

**example_of_results** folder containing examples of result obtained with the dataset


<br /><br />
## HOW TO RUN

See Protocol 136 (ask me).


### If error message

If an error message appears, like:
```
Unknown error accessing project `gmillot/fisher_for_vcf` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/fisher_for_vcf
```
Purge using:
```
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
```


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
Launching `main.nf` [loving_morse] - revision: d5aabe528b
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
nextflow run vcf_slivar_extraction.nf -C vcf_slivar_extraction.config
```

with -C to specify the name of the config file used.


### From a cluster using a committed version on gitlab

Start with:

```bash
EXEC_PATH="/pasteur/zeus/projets/p01/BioIT/gmillot/fisher_for_vcf" # where the bin folder of the main.nf script is located
export CONF_BEFORE=/opt/gensoft/exe # on maestro

export JAVA_CONF=java/13.0.2
export JAVA_CONF_AFTER=bin/java # on maestro
export SINGU_CONF=singularity/3.8.3
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
# distant main.nf file
/fisher_for_vcf/" ; nextflow run --modules ${MODULES} -hub pasteur gmillot/fisher_for_vcf -r v1.0 -c $HOME/nextflow.config ; HOME="/pasteur/appa/homes/gmillot/"

# local main.nf file ($HOME changed to allow the creation of .nextflow into /$ZEUSHOME/fisher_for_vcf/. See NFX_HOME in the nextflow soft script)
HOME="$ZEUSHOME/fisher_for_vcf/" ; nextflow run --modules ${MODULES} fisher_for_vcf.nf -c fisher_for_vcf.config ; HOME="/pasteur/appa/homes/gmillot/"
```

If an error message appears, like:
```
Unknown error accessing project `gmillot/fisher_for_vcf` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/fisher_for_vcf
```
Purge using:
```
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
```

<br /><br />
## OUTPUT

**reports** folder containing all the reports of the different processes including the **fisher_for_vcf.config** file used

**Miami.png** miami plot in the .png format

**fisher.tsv** 

| Column | Description |
| --- | --- |
| **chr** | Chromosome name |
| **position** | Position on the chromosome (in bp) |
| **ref** | Reference allele sequence |
| **alt** | Alternative allele (ALT) sequence |
| **gene** | Gene name |
| **severity** | Severity of the alternative allele from VEP (*consequence* field) |
| **impact** | Impact of the alternative allele from VEP (*IMPACT* field) |
| **aff** | Count of the number of 0=HOM_REF ; 1=HET ; 2=UNKNOWN ; 3=HOM_ALT in the affected cases. Example: {0:4, 1:2} for 4 cases HOM_REF and 2 cases HET |
| **una** | as in *aff* in the unaffected cases |
| **OR** | Odds ratio (n11/n12)/(n21/n22) with:<br /><ul><li>n11 = nHET_aff + nHOM_ALT_aff<br /></li><li>n12 = nHOM_REF_aff<br /></li><li>n21 = nHET_unaff + nHOM_ALT_unaff<br /></li><li>n22 = nHOM_REF_unaff<br /></li>
OR > 1 meaning OR in favor of HET+HOM_ALT/aff |
| **p-value** | p-value of the exact fisher test |
| **-log10(pval)** | -log10 of the p-value |
| **AN** |  |


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

### v1.0

1) everything




