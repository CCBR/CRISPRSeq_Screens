# CCBR CRISPRSeq Screen Framework
This is a repository for analyzing CRISPR screening data generated from CRISPR KO libraries like the GeCKO library (https://doi.org/10.1038/nmeth.3047) and the TKO library (https://doi.org/10.1016/j.celrep.2019.02.041). 

Initial alignment to the reference library is performed using MAGeCK, and tests for gene essentiality can be conducted using either `MAGeCK` or `BAGEL2`. Further analysis with drug treatment can performed with the `drugZ` package.

Note: Any text contained within diagonal brackets `<like this>` indicates a field that the user fills out.

## Setup to run on the NIH Biowulf HPC Cluster via MacOS 
Open a terminal window and log into Biowulf 

```ssh -Y biowulf.nih.gov```

Change to a desired directory:

```cd <directory-path>```

For example

```cd /data/username```

Initialize an interactive node with a tunnel

```sinteractive --mem=64g --cpus-per-task=8 --tunnel --time=24:00:00```

Load required packages directly from Biowulf

```
module load python
module load mageck
module load mageck-vispr 
```
 
If necessary, install BAGEL2 and drugz from Github

```
git clone https://github.com/hart-lab/bagel
git clone https://github.com/hart-lab/drugz
```


## Running CRISPR screening tools independently

### Running MAGeCK to count sgRNAs from FASTQ files
Locate the following files and and note the path to the files
* FASTQ raw reads file
* sgRNA library file

The sgRNA library file must be 1-to-1 (1 gRNA for each PAM binding site), with a total of 3 columns separated by tabs:
1. sgRNA label
2. PAM sequence
3. Gene target

Example:

```
s_10007	TGTTCACAGTATAGTTTGCC	CCNA1
s_10008	TTCTCCCTAATTGCTTGCTG	CCNA1
s_10027	ACATGTTGCTTCCCCTTGCA	CCNC
s_10035	AGAGACCAGCCCGCTGACCG	CCND2
s_10164	GCAGGCGGTACTCAAGGGCA	CCS
s_10200	TTAGAGAAGATCCATCATTC	CCT7
s_10232	AACACGACAGACTTCTGTTC	CD164
s_10264	GAGTCACAGGACGCCCTGAT	CD1D
s_10340	CACGGCTCTGTCACCATCAC	CD276
s_1035	ACACTTGTCATCCGCCTTCA	ADAMTS14
```

Sample FASTQ files may be compressed as `tar.gz`. 

Run the following command: 

```
mageck count -l <library_file> -n <output_label> --sample-label <comma,list,of,labels>  \
    --fastq <fastqFile1> <fastqFile2> … <fastqFileN>
```

The outputs will be contained in the current working directory unless a different path is specified in the output label, e.g. `<path/to/label>`. All output files will start with the label, and include count summary R files, a count summary text file, the raw counts file and normalized counts file.

### Running MAGeCK to test sgRNAs for differential expression between conditions
This step can be run from a counts file, as long as the counts file matches the output format from the previous alignment step, which is typically a tab-delimited text format with sgRNA, gene, and sample labels. Example:

```
sgRNA		Gene		LX	CTRL
s_47512	    RNF111	    1	0
s_24835	    HCFC1R1	    1	0
s_14784	    CYP4B1	    4	0
s_51146	    SLC18A1	    1	0
s_58960	    TRIM5		1	0
s_48256	    RPRD2		1	0
s_30297	    KRTAP5-5	1	0
s_14555	    CYB5B		1	0
s_39959	    PAAF1		1	1
```

After identifying the counts file, run the command: 

```
mageck test -k <countFile.txt> -t <TreatmentSampleLabel> -c <ControlSampleLabel> -n <output_label>
```

The treatment and control labels need to match the labels used in the counts file; in the example above, the treatment label is `LX` and the control label is `CTRL`. It is strongly recommended to make the output label different from the one used for the counts file.

The output files generated include `<output_label>.gene_summary.txt` and `<output_label>.sgrna_summary.txt`. The sgRNA summary table contains the sgRNA ID, the target gene, and p-values for each sgRNA for positive or negative expression. The gene summary contains similar data, but also evaluates the cumulative effects of the sgRNAs targeting the same gene.

### Running MAGeCK to test sgRNA targets for essentiality using MLE
MAGeCK can use a maximum likelihood estimation (MLE) model to determine how essential genes are for a survival study. This code uses a raw counts file, which can also be a comma-separated value (CSV) file, and a design matrix file. A sample design matrix file is shown here:

```
Samples		baseline	HL60_HAEMATOPOIETIC_LYMPHOID_TISSUE	KBM7
HL60.initial	1		0							0
KBM7.initial	1 		0							0
HL60.final		1		1							0
KBM7.final		1		0							1
```

The baseline needs to be set at 1 for all samples. Each sample then has a single baseline and a flag for a different stage to be evaluated against the baseline.

The command to run the MLE for essentiality is:	

```
mageck mle -k <counts_file> -d <design_matrix.txt> -n <output_label>
```

The command returns a sgRNA summary and a gene summary. The sgRNA summary includes an estimate of gRNA efficiency. The gene summary is more detailed and contains additional information. The main addition is the inclusion of the beta score, which indicates the direction of selection for a gene (positive beta scores are associated with positive selection). These statistics are calculated for each gene and each comparison listed in the columns of the design matrix.

More complex designs are described here: https://sourceforge.net/p/mageck/wiki/advanced_tutorial/#tutorial-4-make-full-use-of-mageck-mle-for-more-complicated-experimental-design-eg-paired-samples-time-series

### Running BAGEL for gene essentiality
BAGEL requires two steps, starting from a counts file. The authors recommended using MAGeCK to create the counts file, or possibly Bowtie1 or Bowtie2. The fold change file can be created using the following command:

```
python BAGEL.py fc -i <counts_file.txt> -o <outputFileLabel> -c <controlLabel>
```

The resulting fold change file, with the label <outputFileLabel>.foldchange, is not as detailed as the one provided by MAGeCK, as it contains three columns: 

```
REAGENT_ID	GENE		LX
s_47512	    RNF111	    0.265
s_24835	    HCFC1R1	    0.265
s_14784	    CYP4B1	    0.850
s_51146	    SLC18A1	    0.265
s_58960	    TRIM5		0.265
s_48256	    RPRD2		0.265
s_30297	    KRTAP5-5	0.265
s_14555	    CYB5B		0.265
s_39959	    PAAF1		0.002
```

BAGEL will also generate a normalized read count file, which I will be examining for consistency against the MAGeCK normalized read counts:

```
sgRNA		LX		CTRL
s_47512	4224.46	3515.93
s_24835	4224.46	3515.93
s_14784	6336.69	3515.93
s_51146	4224.46	3515.93
s_58960	4224.46	3515.93
s_48256	4224.46	3515.93
s_30297	4224.46	3515.93
s_14555	4224.46	3515.93
s_39959	4224.46	4219.11
```

The essential gene function calculates the log2 Bayes factor (BF) for each gene, through this command:

```
python BAGEL.py bf -i <fold_change_file> -o <outputFileName> -e <essentialGeneList> -n <nonessentialGeneList> -c <whichColumnsToTest>
```

The essential and nonessential gene lists are used to train the dataset prior to running the Bayesian analysis and are provided by the BAGEL team within their directories as CEGv2.txt and NEGv1.txt, respectively. The columns to test are the numeric column IDs for the individual fold change comparisons and can be a comma-separated list to collapse multiple replicates. 

The resulting file contains Bayes factors for each gene, where positive BFs indicate that the gene is essential:

```
GENE		BF
TTLL6		-4.654
PAX3		-4.072
NEURL4	    3.122
TRA2B		3.104
ZNF506	    -3.364
ZMYM3		-3.326
PSMB8		-3.041
PLEKHG5	    -3.671
NUP98		-3.048
MCOLN3	    -3.671
```

The Bayes factors can then be inserted into a precision-recall curve to evaluate accuracy, using the essential and nonessential genes as the values.

```
python BAGEL.py pr -i <bf_file> -o <output_pr_fileName> -e <essentialGeneList> -n <nonessentialGeneList>
```

The output is a file with precision and recall, which can then be used as graphical inputs.

```
Gene		BF		    Recall	    Precision	FDR
NEURL4	    3.122		0.000		1.000		0.000
TRA2B		3.104		0.000		1.000		0.000
ST8SIA4	    2.855		0.000		1.000		0.000
RUFY3		2.805		0.000		1.000		0.000
FABP3		2.629		0.000		1.000		0.000
PPAP2C	    2.593		0.000		1.000		0.000
TMPRSS11E	2.440		0.000		1.000		0.000
ZNF611	    2.090		0.000		1.000		0.000
TNFSF12	    2.067		0.000		1.000		0.000
```

More details are included tutorial in the BAGEL directory: BAGEL-v2-tutorial.html.


### Running drugZ
drugZ is another program created by the Hart group at UToronto. The primary goal of drugZ is to evaluate the difference in sgRNA counts and if the targeted gene operates synergistically with the drug being tested or suppresses its effects. The base command is as follows: 

```
python drugz.py -i <countsFile> -r <genesToRemove> -c <controlLabels> -x <experimentalDrugLabels> -o <outputFileName> --half_window_size=<num>
```

The counts file can contain multiple replicates and can indicate paired (default) or unpaired experiments (addressed by using the `--paired` flag with TRUE or FALSE). A sample file is shown here (with the header row flowing to the next line):

```
sgRNA	Gene	T0	T15_A_control	T15_B_control	T15_C_control	T15_A_olaparib	T15_B_olaparib	T15_C_olaparib
A1BG_CACCTTCGAGCTGCTGCGCG	A1BG	313	235	47	337	428	115	340
A1BG_AAGAGCGCCTCGGTCCCAGC	A1BG	99	8	1	13	26	5	28
A1BG_TGGACTTCCAGCTACGGCGC	A1BG	650	336	74	185	392	193	304
A1BG_CACTGGCGCCATCGAGAGCC	A1BG	718	192	34	296	178	69	185
A1BG_GCTCGGGCTTGTCCACAGGA	A1BG	180	230	29	122	394	148	364
A1BG_CAAGAGAAAGACCACGAGCA	A1BG	428	300	158	294	366	184	489
A1CF_CGTGGCTATTTGGCATACAC	A1CF	677	452	74	423	585	446	434
A1CF_GGTATACTCTCCTTGCAGCA	A1CF	138	69	43	109	96	184	127
A1CF_GACATGGTATTGCAGTAGAC	A1CF	396	183	38	106	193	120	198
```

The genes to remove, as shown in the example in the tutorial, include LacZ, luciferase, and EGFR and are input as comma-separated list (e.g. `LacZ,luciferase,EGFR`). The controls and experimental samples should match the column headers and be input as a comma separated list as well (e.g. `-c T15_A_control,T15_B_control,T15_C_control`). The half_window_size is the width of the variance window (how many sgRNAs used to calculate variance at any given time) and is set at a default value of 500, but can be adjusted to match dataset size; the recommendation is to make this value the “size of the first bin and half the size of the initial window.”  There is also the option to include a fold change text file, if desired.

The outputs are in a single text file:

```
GENE	sumZ	numObs	normZ	pval_synth	rank_synth	fdr_synth	pval_supp	rank_supp	fdr_supp
A1CF	1.64	9		-1.00	0.159		1		0.317		0.841		2		0.841
A1BG	3.10	18		1.00	0.841		2		0.841		0.159		1		0.317
```

The values include the total Z-score for the individual genes, the number of observations, the normalized Z-score, and significance statistics for the gene for either synergistic synthetic lethality (`synth`) or suppressive (`supp`) effects.

More details are included in the tutorial in the directory: drugZ_in_jupyter_notebook_tutorial.html


## References
MAGeCK
* Tutorial: https://sourceforge.net/p/mageck/wiki/demo/
* Biowulf HPC: https://hpc.nih.gov/apps/MAGeCK.html 
* Reference paper: https://doi.org/10.1186/s13059-014-0554-4 

MAGeCK-VISPR
* Tutorial: https://bitbucket.org/liulab/mageck-vispr/src/master/ 
* Biowulf HPC: https://hpc.nih.gov/apps/mageck-vispr.html 
* Reference paper: https://doi.org/10.1186/s13059-015-0843-6 

BAGEL2
* Tutorial: Included with package
* Github: https://github.com/hart-lab/bagel 
* Reference paper: https://doi.org/10.1186/s13073-020-00809-3

drugZ
* Tutorial: Included with package
* Github: https://github.com/hart-lab/drugz 
* Reference paper: https://doi.org/10.1186/s13073-019-0665-3 

