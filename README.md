# DiscoDivas

## Disambiguation
This page is for "a Genetic ***Dis***tance-assisted PRS ***Co***mbination Pipeline for ***Div***erse Genetic ***A***ncestrie***s***, an R script for calculating PRS for diverse populations, especially admixed populations. 

## Required R packages
- data.table
- dplyr
- stringr
- rio
- optparse

## Command line to run DiscoDivas.R
Here is an example command line to run DiscoDivas.R

```
Rscript DiscoDivas.R -m file/med.g1000G.4pop.tsv \
-p ukbb.pca.txt.gz \
--prs.list afr.prs.tsv.gz,${wd}/eas.prs.tsv.gz,eur.prs.tsv.gz,sas.prs.tsv.gz \
-s IID,PRS \
-A 1,0.9,1,1 \
--regress.PCA T \
-o combine.linear.prs
```
- `-m ukbb.pca.med.1000G.tsv`
  
  `med.g1000G.4pop.tsv`can be found in the `files`. **Please notice that the order of the list after `-A` and `--prs.list` should has the same the order as the order of cohorts/ populations in this file.** User should replace the file by other file that contains the median value of top PC of the validation cohorts actually being used in the analysis.
  
  The format should be: 
  - must have a header; does't contain row name or row index;
  - the first column is the name of cohorts or populations; the seconds column and other columns after that are PCs
```
POP	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10
AFR	-10358.4	837.239	193.871	-103.469	195.704	64.7022	-400.128	-49.6729	0.0890048	106.919
EAS	2707.46	6490.42	-166.124	-1408.4	134.848	-0.94705	-193.3735	-362.328	8.36484	598.289
EUR	1280.17	-5075.18	-794.419	-1903.47	172.565	224.566	-208.995	-321.358	-5.28843	113.007
SAS	1276.87	-1124.26	4017.77	1557.19	200.125	25.1516	-195.331	-349.297	-417.462	45.8134
```
    
- `-p ukbb.pca.txt.gz`

  `ukbb.pca.txt.gz` should be replaced by user's file contains the PC of the testing individual. It can be plain text file or `.gz` file. **Please notice that in all the list, items should be separated by `,` and no space should be inserted.**

  The format is:
  - must have header; does't contain row name or row index;
  - the first column is the individual ID; the seconds column and other columns after that are PCs

- `--prs.list afr.prs.tsv.gz,${wd}/eas.prs.tsv.gz,eur.prs.tsv.gz,sas.prs.tsv.gz`

  The list of PRS files of the testing cohort that are based on the weight fine-tuned in each of the validation cohorts. They can be plain text file or `.gz` file. The format is:
  - must have header; does't contain row name or row index;
  - must have a column of individual ID and a column of PRS
  - The individual ID in the PRS files should be the same as the individual ID in the PCA file

- `-s IID,PRS`
  The column names of individual ID and PRS in the PRS files. If the column names are the different in different files, please list them in the same order of `prs.list` : `ID1,PRS1,ID2,PPRS2,ID3,PRS3,ID4,PRS4`

- `-A 1,0.9,1,1`

  The list of shrinkage parameter for each PRS based on their quality. The default is a list of 1

- `--regress.PCA T`

  Choose whether combining PRS after correcting the ancestry information or not. Default and recommanded value is TRUE/T

- `-o`

  The prefix of output file. The format of output file is `.tsv.gz`


## Supporting files
Please visit `files` for the suppporting files







