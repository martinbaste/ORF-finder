##ORF finder

This program takes a fasta nucleotide file and outputs the list of possible ORF's with info in gff3 format.
It includes: GC content in all frames, Shine-Dalgarno sequence presence, and codon bias.

Currently it doesn't take into account ambiguity bases.

The output from this program could be used to predict genes from procariotic genomes.

Usage: python3 blastParser.py [options] -i FILE
Options:
    -h : print this help
    -v : verbose mode
    -i [FILE] : path to input file (in fasta)
    -o [FILE] : name of the output file (gff format)
    -m [INT] : Shortest ORF threshold (default = 200)
    -n [INT] : Longest ORF threshold (default = 3000)
    -s : Don't filter ORFs with Shine-Dalgarno probability of 0