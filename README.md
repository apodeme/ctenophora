# ctenophora
Whilst working on the Ctenophora root project, I developed numerous scripts in both Python and Bash. Some of these have general applications, such as batch converting many fasta files to one line format, whilst others are more specific, such as the the Agalma pipeline cleaner, which removes redunant sequences produced by postassemble routine in Agalma after translation of DNA sequences to amino acids.

dotreplace.py - input: newick gene tree(s) in text file - give filename as argument. Replace "." in titles of genes with "DOTHERE" -> ALE does not function correctly otherwise

mcldist.py - input: mcl tab file - give filename as argument. gives mcl cluster distribution, outputs into text file. copy into the included excel file.

faoneline.sh - batch convert fasta files to one line. requires ./oneline/wolineone folders. also requires list of files to be converted, use ls -1 *.afa > taxafilenames.txt

move.sh - moves fasta files with less than four sequences (e.g. after trimming or removing duplicates). requires ./lessthanfour directory and fasta files must be in one line format.

removeasterisk.py - input: fasta file. Removes asterisk from end of sequences - otherwise Muscle throws a warning.

removedupeseqs.py - input: fasta file. give filename as argument (can be ran in parallel). fasta headers must contain species name before the gene followed by a forward slash, e.g. "abylopsis/gene1". removes duplicate sequences only if species is same.

agalmacleaner.py - input: fasta file. cleans output of Agalma pipeline such that there is only one copy of each gene in the resulting file. non-overlapping genes are concatenated together if same sense. 


