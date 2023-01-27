# GOTCHA: an automated workflow for eDNA target capture bait design.

This is the repository for the automated workflow for environmental DNA target capture bait design - “<i>Gotcha</i>”. This tool is leveraging a phylogenetic approach, which is expected to be more solid with skewed datasets and incomplete taxon sampling. A manuscript containing a detail description of its functioning will be available here.

Developed by Kevin Nota (kevin_nota@eva.mpg.de) & Giobbe Forni (giobbe.forni@gmail.com)

## Short description of <i>Gotcha</i>
By default, Gotcha is using <i>BOLD-CLI</i> for downloading standard barcoding genes such as COI, <i>rbc</i>L, and <i>mat</i>K, followed by filtering and selection of the fragment of the barcoding gene which is covered the most. This will create a “clean” multiple sequence alignment (MSA) from which a gene tree is inferred which is then used to reconstruct ancestral state sequences. These sequences are then processed to only the node/tip sequences required for capturing the fast genetic diversity in the original MSA. The tool allows costume inputs for most major steps, which allows the use of non-standard barcoding genes or manually improve the MSA before building gene trees etc. See below the schematic for <i>gotcha</i>.

<p align="center">
<img src="https://github.com/Kevinnota/gotcha/blob/main/documentation/workflow.jpg" data-canonical-src="https://github.com/Kevinnota/gotcha/blob/main/documentation/workflow.jpg" width="650" height="650" />
</p>

# Install <i>Gotcha</i>

Here you will find the guide to install gotcha - [install dependencies](https://github.com/Kevinnota/gotcha/blob/main/documentation/1.md)

# Quick start to <i>Gotcha</i>

The most default way to run <i>gotcha</i> only requires the -t (--taxa), -m (--marker) and -c (--code) flags. This will download the specified marker from BOLD for the specified taxonomic name, and the [genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) belonging to the marker. The available markers are ["COI-5P", "COI-3P", "rbcL", "matK", and "ITS"]. See the example below, and more detailed for more tutorials for downloading sequences using BOLD [bold download usage](https://github.com/Kevinnota/gotcha/blob/main/documentation/2.md).
<br>
<br>

>`The code below will download all sequences from the genus Circus (Harriers).`
 
```
python3 gotcha.2.2.2.py -t Circus -m COI-5P -c 2 
```
><i>A short summary of the progress of gotcha is printed in the terminal, see example below. The output of gotcha will be stored in a directory "probes" (default, can be specified with -o). The most important statistics are stored in the log.txt file, the baits are printed into the .bai, the clean filtered alignment .aln, and basic statistic and gene tree in the summary_stats.pdf.</i>

<br>

>`printed in terminal`:

```
analysis started on 2023-01-24 16:07:06

downloading COI-5P sequences for taxa Circus
#	downloaded:  61
skipping taxonomic filtering
#	kept:  60
collapsing identical sequences with a percent identitiy of 1
#	kept:  43
removing sequences with non ATGC nucleotides
#	kept:  33
filtering sequences with a minimum length of 100
#	kept:  33
aligning sequences
removing sequences with stop codons
#	kept:  33
selecting optimal fragment size of 33 sequences with max length of 1548
calculating missing data
100%|██████████████████████████████| 1548/1548 [00:00<00:00, 1855.85it/s]
finding maximal sequence fragment
100%|████████████████████████████| 101/101 [00:00<00:00, 173.42it/s]
optimal missing data threshold = 4%
number of sequences retained = 32
number of sequences omitted = 1
selected fragment size = 618
start position = 84
end position = 702
inferring tree
inferring ancestral sequences using genetic code 2
finding baits of 80 nt with max distance of 0.09 % and tiling 10
```


<br>

>`log.txt`
```
# analysis started on 2023-01-24 16:07:06


	INPUTS:

	 BOLD TAXONOMY 		Circus
	 GEOGRAPHY 		N
	 CUSTOM SEQUENCES 	N
	 CUSTOM ALIGNMENT 	N
	 CUSTOM TREE 		N
	 MARKER 		COI-5P
	 GEN CODE 		2
	 FAST MODE 		N


	 PARAMETERS:

	 TAXONOMIC FILTER 	N
	 MIN MARKER LENGTH 	100
	 MAX MARKER IDENTITY 	1
	 BAITS LENGTH 		80
	 BAITS TILING 		10
	 TRIMAL 		N
	 DUST TRESHOLD 		1


	 FILTERING:

	 61	 starting sequences
	 60	 after taxonomic filtering
	 43	 after identity collapse
	 33	 after removing sequences with non ATGC nucleotides
	 33	 after length filter
	 33	 after the removal of sequences with stop codons
	 32	 after size selection


	 BAITS:

	 162	 baits
	 128	 baits passing gc+complexity filter
	 96	 unique baits

```

>`bait.fasta`
```
>BISE004-07 0:80 GC=0.55; GC_flag=good; dust_flag=pass;
ATAGCCGGCACCGCCCTCAGTCTACTCCTTCGTGCAGAACTCGGTCAACCAGGCACCCTTCTAGGTGATGACCAAATCTA
>node #61 0:80 GC=0.5375; GC_flag=good; dust_flag=pass;
ATAGCCGGCACCGCCCTCAGTCTACTCATTCGTGCAGAACTCGGTCAACCAGGCACCCTTCTAGGTGATGACCAAATCTA
>GBIR9581-19 0:80 GC=0.5375; GC_flag=good; dust_flag=pass;
ATAGCCGGCACCGCCCTCAGTCTACTCATTCGTGCAGAACTCGGTCAACCAGGCACCCTTCTAGGTGATGACCAAATCTA
>node #36 0:80 GC=0.525; GC_flag=good; dust_flag=pass;
ATAGTCGGCACCGCCCTTAGCCTACTCATTCGCGCAGAACTTGGTCAACCAGGCACACTCCTAGGTGATGACCAAATCTA
```

>`summary_stats.pdf`

<p align="center">
<img src="https://github.com/Kevinnota/gotcha/blob/main/example_files/Plot_E_850_350.png.svg" data-canonical-src="https://github.com/Kevinnota/gotcha/blob/main/example_files/Plot_A_850_350.svg" width="650" height="150" />

<img src="https://github.com/Kevinnota/gotcha/blob/main/example_files/Plot_E_850_350.png.svg" data-canonical-src="https://github.com/Kevinnota/gotcha/blob/main/example_files/Plot_B_850_350.svg" width="650" height="150" />

<img src="https://github.com/Kevinnota/gotcha/blob/main/example_files/Plot_E_850_350.svg" data-canonical-src="https://github.com/Kevinnota/gotcha/blob/main/example_files/Plot_C_850_350.svg" width="650" height="150" />

<img src="https://github.com/Kevinnota/gotcha/blob/main/example_files/Plot_E_950_950..svg" data-canonical-src="https://github.com/Kevinnota/gotcha/blob/main/example_files/Plot_E_850_350.svg" width="650" height="650" />
</p>

# <i>Gotcha</i> with costume input files 
In some cases it might be required to run <i>gotcha</i> with costume input files, such as non-standard marker gene, or in house reference sequences. See here for the costume files that can be used - [custom inputs usage](https://github.com/Kevinnota/gotcha/blob/main/documentation/3.md)

# Parameters tweaking
<i>Gotcha</i> has a lot of parameters that can be changed to allow for general usage using for example larger bait/probe sizes. For more details see here - [parameters tweaking](https://github.com/Kevinnota/gotcha/blob/main/documentation/4.md)
