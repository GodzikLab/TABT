# TABT

Therapeutic AntiBody Targets Score

for details see the manuscript:
Integration of expression, genomic, and protein interaction data suggests new target possibilities in autoimmune diseases
by Lukasz Jaroszewski, Carl Ware, and Adam Godzik

Software and data

data/ - all input data for TABT

data/abs/ - lists of therapeutic antibodies used to evaluate accuracy of TABT
data/expression/ - results from differential expression 
data/expression/*/*.top.table.tsv - limma results for individual GEO series
data/expression/*_consensus.tsv - expression consensus lists

data/gwas/ - lists of results from GenomeWide Association Studies
data/gwas/*.tsv - association data dowloaded from GWAS catalog https://www.ebi.ac.uk/gwas/
data/gwas/*_genes.tsv - lists of genes with polymorphisms strongly associated with a disease (p-vale < 1e-20)

data/interactions/interactions.tsv - list of protein-protein interactions prepared based on data
	dowloaded from STRING database (https://string-db.org/cgi/download)
	The downloaded input files (currently 9606.protein.links.v11.5.txt and 9606.protein.aliases.v11.5.txt)
	are not included in the github repository due to size limits.

data/location/uniprot.tsv - SWISSPROT data for human proteins
	(downloaded from https://www.uniprot.org/uniprotkb) 

data/tissues/E-MTAB-513-query-results.tpms.tsv - data on expression in human tissues
	from publicly available Illumina Body Map (https://www.ebi.ac.uk/gxa/experiments/E-MTAB-513/Results)

results/ - directory where TABT scores and TABT evaluation are saved

scripts/ - all python programs used by TABT
scripts/prep_*.py - scripts to prepare different types of data used by TABT
scripts/tabt.py - the scripts which calculates tabt scores and evaluates them based on testing set of
	of therapeutic antibodies

To run TABT scripts, please create the environment variable TABT_DIR with a name of the
current TABT directory e.g., by adding the following line to your .bashrc in your home directory:

export TABT_DIR="<full_path_of_your_tabt_directory>"
