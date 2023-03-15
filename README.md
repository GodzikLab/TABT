# TABT

<pre>
Therapeutic AntiBody Targets Score

Supporting data and scripts for the manuscript:
Integration of omics data suggests new antibody targets in autoimmune diseases
by Lukasz Jaroszewski, Carl Ware, and Adam Godzik

Contents:

data/ - all input data for TABT

data/abs/ - lists of therapeutic antibodies used to evaluate accuracy of TABT

data/expression/ - differential expression results
data/expression/_disease_/_geo_accession_.top.table.tsv - limma results for individual GEO series
data/expression/_disease__consensus.tsv - consensus lists for differential expression

data/gwas/ - results from GenomeWide Association Studies
data/gwas/_disease_.tsv - association data dowloaded from GWAS catalog https://www.ebi.ac.uk/gwas/

data/gwas/_disease__genes.tsv - lists of genes with polymorphisms strongly associated with a disease (p-value < 1e-20)

data/interactions/interactions.tsv - list of protein-protein interactions prepared based on data from the STRING database.
	The input files (currently 9606.protein.links.v11.5.txt and 9606.protein.aliases.v11.5.txt)
	are not included in the github repository due to the github's size limits -
        they are downloaded from STRING by the parse_interactions.py script instead.

data/location/uniprot.tsv - SWISSPROT data for human proteins
	(downloaded from https://www.uniprot.org/uniprotkb) 

data/tissues/E-MTAB-513-query-results.tpms.tsv - data on expression in human tissues
	from publicly available Illumina Body Map (https://www.ebi.ac.uk/gxa/experiments/E-MTAB-513/Results)

results/ - folder where TABT scores and TABT evaluation are saved
results/_disease__tabt.tsv - full list of intermediate results, TABT score and ranking for _disease_
results/evaluation.tsv - evaluation of TABT score based on rankings of verified antobody targest for ibd,pso,ra and sle

scripts/ - all python programs used by TABT
scripts/prep_*.py - preprocess expression, GWAS and interactions data used by TABT
scripts/tabt.py - calculates TABT scores and evaluates them based on testing sets of therapeutic antibodies

To run TABT scripts, please create the environment variable TABT_DIR containing the name of the
current TABT folder e.g., by adding the following line to your .bashrc in your home folder:

export TABT_DIR="_full_path_of_tabt_folder_"
</pre>

