# GOTCHA: an automated workflow for eDNA target capture bait design.

This is the repository for the automated workflow for environmental DNA target capture bait design - “<i>Gotcha</i>”. This tool is leveraging a phylogenetic approach, which is expected to be more solid with skewed datasets and incomplete taxon sampling. A manuscript containing a detail description of its functioning will be available here.

Developed by Kevin Nota (kevin_nota@eva.mpg.de) & Giobbe Forni (giobbe.forni@gmail.com)

#### Quick links
-  [Install dependencies](https://github.com/Kevinnota/gotcha/blob/main/documentation/1.md)
-  [Custom inputs usage](https://github.com/Kevinnota/gotcha/blob/main/documentation/3.md)
-  [Parameters tweaking - flags](https://github.com/Kevinnota/gotcha/blob/main/documentation/4.md)
-  [Tutorials](https://github.com/Kevinnota/gotcha/blob/main/documentation/5.md)
<br>
<br>

## Short description of <i>Gotcha</i>
By default, Gotcha is using <i>BOLD-CLI</i> for downloading standard barcoding genes such as COI, <i>rbc</i>L, and <i>mat</i>K, followed by filtering and selection of the fragment of the barcoding gene which is covered the most. This will create a “clean” multiple sequence alignment (MSA) from which a gene tree is inferred which is then used to reconstruct ancestral state sequences. These sequences are then processed to only the node/tip sequences required for capturing the fast genetic diversity in the original MSA. The tool allows costume inputs for most major steps, which allows the use of non-standard barcoding genes or manually improve the MSA before building gene trees etc. See below the schematic for <i>gotcha</i>.

<p align="center">
<img src="https://github.com/Kevinnota/gotcha/blob/main/documentation/workflow.jpg" data-canonical-src="https://github.com/Kevinnota/gotcha/blob/main/documentation/workflow.jpg" width="650" height="650" />
</p>
<br>
<br>

# Install <i>Gotcha</i>

The easiest way to install gotcha is using the gotcha_env.yml file either with conda or mamba.

```
mamba env create -f gotcha_env.yaml
```
or
```
conda env create -f gotcha_env.yaml
```

 Here you will find a list of all the dependencies and a guide to install gotcha manually - [install dependencies](https://github.com/Kevinnota/gotcha/blob/main/documentation/1.md)
<br>
<br>
# Quick start to <i>Gotcha</i>

The most default way to run <i>gotcha</i> only requires the -t (--taxa), -m (--marker) and -c (--code) flags. This will download the specified marker from BOLD for the specified taxonomic name, and the [genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) belonging to the marker. The available markers are ["COI-5P", "COI-3P", "rbcL", "matK", and "ITS"]. See the example below, and more detailed for more tutorials for downloading sequences using BOLD [bold download usage](https://github.com/Kevinnota/gotcha/blob/main/documentation/2.md).
<br>
<br>

><i>`The code below will download all sequences from the family Strigiformes (owls).`</i>
 
```
python3 gotcha.py -m "COI-5P" -t Strigiformes -c 2 -o Strigiformes
```
><i>A short summary of the progress of gotcha is printed in the terminal, see example below. The output of gotcha will be stored in a directory "probes" (default, can be specified with -o). The most important statistics are stored in the log.txt file, the baits are printed into the .bai, the clean filtered alignment .aln, and basic statistic and gene tree in the summary_stats.pdf.</i>

```
Run Gotcha Main

analysis started on 2023-03-27 16:15:39

downloading COI-5P sequences for taxa Strigiformes
#        downloaded:  794
skipping taxonomic filtering
#        kept:  743
collapsing identical sequences with a percent identitiy of 1
#        kept:  478
removing sequences with non ATGC nucleotides
#        kept:  360
filtering sequences with a minimum length of 100
#        kept:  360
removing sequences with stop codons
#        kept:  356

selecting optimal fragment size based on 356 :
#       kept:			 335
#       omitted:         21
#       fragment size:   611


inferring tree

inferring ancestral sequences using genetic code 2

finding baits of 80 nt with max distance of 0.09 % and tiling 30
Scanning tile : 18it [00:02,  6.23it/s]


Gotch run successfully

a totall of 397 baits are needed to capture the selected taxa


analysis finished on 2023-03-27 16:16:33

```
Gotcha is producing a number of files. The final baits are stored in <i>*_filtered_uniq_clustered_baits.fasta</i>. The sequence header starts with the node in the <i>*.nwk</i>

```
>node #668 0:80 GC=0.5375; GC_flag=good; dust_flag=pass;
GCCCTCAGCCTGCTCATCCGAGCTGAACTAGGCCAACCAGGCACACTACTCGGCGATGACCAAATCTACAACGTAATTGT
```

>`output dir:`
```
ls -la probes
drwxrwxr-x  2 public staff    4096 Mar 27 16:16 ./
drwxrwxr-x 50 public staff    4096 Mar 27 16:15 ../
-rw-rw-r--  1 public staff  211033 Mar 27 16:16 Strigiformes.aln
-rw-rw-r--  1 public staff   13341 Mar 27 16:16 Strigiformes.nwk 
-rw-rw-r--  1 public staff   55981 Mar 27 16:16 Strigiformes_filtered_baits.fasta
-rw-rw-r--  1 public staff   55841 Mar 27 16:16 Strigiformes_filtered_uniq_baits.fasta
-rw-rw-r--  1 public staff   74969 Mar 27 16:16 Strigiformes_raw_baits.fasta
-rw-rw-r--  1 public staff     800 Mar 27 16:16 log.txt
-rw-rw-r--  1 public staff 5150703 Mar 27 16:16 node_tree.svg
-rw-rw-r--  1 public staff   42691 Mar 27 16:16 summary_stats.pdf
```
<br>

>`log.txt`
```
# analysis started on 2023-03-27 16:15:39


         INPUTS:

         BOLD TAXONOMY          Strigiformes
         GEOGRAPHY              N
         CUSTOM SEQUENCES       N
         CUSTOM ALIGNMENT       N
         CUSTOM TREE            N
         MARKER                 COI-5P
         GEN CODE               2
         FAST MODE              N


         PARAMETERS:

         TAXONOMIC FILTER       N
         MIN MARKER LENGTH      100
         MAX MARKER IDENTITY    1
         BAITS LENGTH           80
         BAITS TILING           30
         TRIMAL                 N
         DUST TRESHOLD          1


         FILTERING:

         794     starting sequences
         743     after taxonomic filtering
         478     after identity collapse
         360     removing sequences with non ATGC nucleotides
         360     after length filter
         356     the removal of sequences with stop codons
         335     after size selection


         BAITS:

         513     baits
         398     baits passing gc+complexity filter
         397     unique baits
         397     unique clustered baits

# analysis finished on 2023-03-27 16:16:33
```

>`summary_stats.pdf`

<p align="center">
<img src="https://github.com/Kevinnota/gotcha/blob/main/documentation/Plot_A_1250_350.svg" data-canonical-src="https://github.com/Kevinnota/gotcha/blob/main/documentation/Plot_A_1250_350.svg" width="800" height="250" />

<img src="https://github.com/Kevinnota/gotcha/blob/main/documentation/Plot_B_1250_350.svg" data-canonical-src="https://github.com/Kevinnota/gotcha/blob/main/documentation/Plot_B_1250_350.svg" width="800" height="250" />

<img src="https://github.com/Kevinnota/gotcha/blob/main/documentation/Plot_C_1250_350.svg" data-canonical-src="https://github.com/Kevinnota/gotcha/blob/main/documentation/Plot_C_1250_350.svg" width="800" height="250" />

<img src="https://github.com/Kevinnota/gotcha/blob/main/documentation/Plot_E_1250_1250.svg" data-canonical-src="https://github.com/Kevinnota/gotcha/blob/main/documentation/Plot_E_1250_1250.svg" width="800" height="800" />
</p>
<br>
<br>
# <i>Gotcha</i> with costume input files 
In some cases it might be required to run <i>gotcha</i> with costume input files, such as non-standard marker gene, or in house reference sequences. See here for the costume files that can be used - [custom inputs usage](https://github.com/Kevinnota/gotcha/blob/main/documentation/3.md)
<br>
<br>
# Parameters tweaking
<i>Gotcha</i> has a lot of parameters that can be changed to allow for general usage using for example larger bait/probe sizes. For more details see here - [parameters tweaking](https://github.com/Kevinnota/gotcha/blob/main/documentation/4.md)

<br>
<br>
# Tutorial
For multiple tutorials on how to run gotcha click here.