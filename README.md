# <i>Gotcha</i>
<div style="text-align: justify"> This is the repository for the automated workflow for environmental DNA target capture bait design - “<i>Gotcha</i>”. By default, Gotcha is using <i>BOLD-CLI</i> for downloading standard barcoding genes such as COI, <i>rbc</i>L, and <i>mat</i>K, followed by filtering and alignment steps, and selection of the fragment of the barcoding gene which is covered the most. This will create a “clean” multiple sequence alignment (MSA) which is used to infer a gene tree which is used to reconstruct ancestral state sequences. These sequences are then processed to only the node/tip sequences required for capturing the fast genetic diversity in the original MSA. The tool allows costume inputs for most major steps, which allows the use of non-standard barcoding genes to be used or manually improve the MSA before building gene trees etc.</div>

# How to run <i>Gotcha</i>

