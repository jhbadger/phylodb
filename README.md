# phylodb
Scripts for maintaining JCVI phylodb database

Most of these scripts are written in Ruby. Some of them are in Julia
for performance reasons. I design the scripts to give usage
information if run without arguments.



writePhyloDB		Writes flatfile from mysql phylodb database
splitPhyloDB.jl		script for splitting phylodb into sections for sub-blasts	
loadPhyloDB		script for adding organism(s) to phylodb given contigs.txt, proteins.txt, optionally transcripts.txt, geneorders.txt rrnas.txt