`g1k_hm3_maf5_woamb_wolr.pca.weight` is the PC loading (SNP effect size to generate the PC) used in the DiscoDivas paper. 

`med.1000G.allpop.tsv` is the median of top 10 PCs of QC'ed 1000 Genomes samples of AFR, AMR, EAS, EUR, and SAS. which could be used as the approximated PC values when the actual PC of the fine-tuning sample is unavailable. The PC is calculated based on the PCA loading file `g1k_hm3_maf5_woamb_wolr.pca.weight`

`med.1000G.4pop.tsv` is the median of top 10 PCs of QC'ed 1000 Genomes samples of AFR, EAS, EUR, and SAS, which could be used as the approximated PC values when the actual PC of the fine-tuning sample is unavailable.  This file is the input of the example command line of running DiscoDivas



