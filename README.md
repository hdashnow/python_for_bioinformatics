Python for Bioinformatics
=========================

Some examples of using Python for Bioinformatics

BioPython
---------

###Sequence objects
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec17



###SeqRecord objects
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec32

The SeqRecord class is a more sophisticated way to represent a sequence. It has a number of attributes, but the most useful are usually:  
.seq  
– The sequence itself, typically a Seq object.  
.id  
– The primary ID used to identify the sequence – a string. In most cases this is something like an accession number.  

You can create a SeqRecord object by hand like this:  



However usually, you would obtain a SeqRecord object by reading in a file, such as fasta file.
http://en.wikipedia.org/wiki/FASTA_format

An example sequence in FASTA format is:

>AB000263 |acc=AB000263|descr=Homo sapiens mRNA for prepro cortistatin like peptide, complete cds.|len=368
ACAAGATGCCATTGTCCCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCGCTGCCCTGCC
CCTGGAGGGTGGCCCCACCGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAATAAGGAAAAGCAGC
CTCCTGACTTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGAGAGG
AAGCTCGGGAGGTGGCCAGGCGGCAGGAAGGCGCACCCCCCCAGCAATCCGCGCGCCGGGACAGAATGCC
CTGCAGGAACTTCTTCTGGAAGACCTTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAG
TTTAATTACAGACCTGAA

###Reading and writing sequencing files
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec51

####Resources:
* http://rosalind.info/
* http://biopython.org/DIST/docs/tutorial/Tutorial.html

####Advanced resources:
* http://pysam.readthedocs.org/en/latest/
