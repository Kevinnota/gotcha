# GOTCHA: an automated workflow for eDNA target capture bait design.

This is the repository for the automated workflow for environmental DNA target capture bait design - “<i>Gotcha</i>”. By default, Gotcha is using <i>BOLD-CLI</i> for downloading standard barcoding genes such as COI, <i>rbc</i>L, and <i>mat</i>K, followed by filtering and alignment steps, and selection of the fragment of the barcoding gene which is covered the most. This will create a “clean” multiple sequence alignment (MSA) which is used to infer a gene tree which is used to reconstruct ancestral state sequences. These sequences are then processed to only the node/tip sequences required for capturing the fast genetic diversity in the original MSA. The tool allows costume inputs for most major steps, which allows the use of non-standard barcoding genes to be used or manually improve the MSA before building gene trees etc.

## Install <i>Gotcha</i>
Think/hope this will become
  
```{bash}
conda install gotcha #not working at the moment
```

## How to run <i>Gotcha</i>

The most default way to run <i>gotcha</i> only requires the -t (--taxa), -m (--marker) and -c (--code) flags. This will download the specified marker from BOLD for the specified taxonomic name, and the genetic code belonging to the marker. The available markers are ["COI-5P", "COI-3P", "rbcL", "matK", and "ITS"]. 

<i>The code below will download all sequences from the Pilosa order containing anteaters and sloths.</i>
```{bash}
python3 gotcha.2.2.2.py -t Pilosa -m COI-5P -c 5 
```
<i>The output of the run code looks like this.  </i>
```{bash}
analysis started on 2023-01-18 15:58:28

downloading COI-5P sequences for taxa Pilosa
downloaded 246 sequences
skipping taxonomic filtering
#	kept:  145
collapsing identical sequences with a percent identitiy of 1
#	kept:  76
aligning sequences
removing sequences with stop codons
#	kept:  76
selecting optimal fragment size of 69 sequences with max length of 1554
calculating missing data
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 1554/1554 [00:00<00:00, 1772.82it/s]
finding maximal sequence fragment
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 101/101 [00:00<00:00, 174.54it/s]
optimal missing data threshold = 58%
number of sequences retained = 29
number of sequences omitted = 40
selected fragment size = 1530
start position = 0
end position = 1530
#	kept:  29
inferring tree
inferring ancestral sequences using genetic code 5
finding baits of 80 nt with max distance of 0.09 % and tiling 10
analysis finished on 2023-01-18 15:58:45

```
