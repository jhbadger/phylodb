# phylodb
Scripts for maintaining JCVI phylodb database

Most of these scripts are written in Ruby. Some of them are in Julia
for performance reasons. I design the scripts to give usage
information if run without arguments.

Requirements
============

Ruby (1.9+) gems (if not present type in bash: gem install <gem name>)
----------------------------------------------------------------------
bio
digest/sha1
mysql
net/ftp
open-uri
trollop

Julia (0.3+) packages (if not present type in Julia prompt: Pkg.add("<package name>") 
-------------------------------------------------------------------------------------
ArgParse

Files
-----
The eukaryotic taxonomy is taken from pr2
http://5.196.17.195/pr2/download/entire_database/qiime_gb203_taxo.txt.gz

(uncompress; in fact currently need to uncompress twice as they've gzipped a gzipped file)

Prokaryotic taxonomy from SILVA
http://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_nr_119.txt
