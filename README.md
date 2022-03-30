# Celegans_sensingDyf1
This repository contains code and data necessary to replicate analysis of Dyf1 C.elegans mutants and generates figures from A.Segref manuscript: 
 "Thermosensation is linked /coupled to ubiquitin-dependent protein turnover through insulin and calcineurin signaling"

* **dataAnalysisANDviz_dyf1_publication.r** <br />
Code for replicating analysis and paper figures.

* **allVSdyf1_masspecDEresults.tsv**  <br />
Result of differential expression analysis of the MassSpec data, it is the output from running DEP::run_app("LFQ") on proteinGroups.txt input (see raw_data). 3 contrasts are tested: Dyf1 VS control, Dyf1 VS Dyf1/unc13 and Dyf1 VS Dyf1/unc31.

* **dyf1_marray_DEall.csv** <br />
Result of Microarray differential expression analysis: Dyf1 mutant VS control. This is the output from running limma::lmfit() on rma normalised input.

* **raw_data** <br /> 
A folder that contains Microarray and MassSpec input data.
