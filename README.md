# UV_syntrichia

This repository contains the code and final data (transcript count data) for the manuscript tiled "The transcriptomic effects of acute UV radiation exposure on two Syntrichia mosses" published in X.

The script ```clean_rna.sh``` uses Trimmomatic vs. 0.30 to quality rim raw reads

The script ```tophat_paired-end_acute.sh``` assembles the read to the S. caninervis genome/transcriptome

The scripts ```htseq-count_acute-BR.sh``` and ```htseq-count_acute-SC.sh``` estimate transcript counts for each genome for S. caninervis (SC) and S. ruralis (BR)

Using the count data in R, the script ```timeseries.R``` performs differential abundance/expression analyses and produces summary tables and plots

Finally, the R script ```test_OGs_DPs_ELIPs_LEAs.R``` counts ELIPs and LEAs in the summary tables from the differential abundance analyses. 