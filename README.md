# phylodb
Scripts for maintaining JCVI phylodb database

Most of these scripts are written in Ruby. Some of them are in Julia
for performance reasons. I design the scripts to give usage
information if run without arguments.



writePhyloDB		Writes flatfile from mysql phylodb database
splitPhyloDB.jl		script for splitting phylodb into sections for sub-blasts	
loadPhyloDB		script for adding organism(s) to phylodb given contigs.txt, proteins.txt, optionally transcripts.txt, geneorders.txt rrnas.txt
generateTaxStrings.jl	generate phylodb taxonomy from scratch, using silva, pr2, viral taxonomy. Generally only used for major updates as minor tweaks lost


== needed files ===
The eukaryotic taxonomy is taken from pr2 -- http://ssu-rrna.org/pr2 using the qiime formatted taxonomy

The prokaryotic taxonomy is taken from SILVA (eg. http://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_119.txt)
