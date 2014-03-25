Python for Bioinformatics
=========================

Some examples of using Python for Bioinformatics

Biopython
---------

###Sequence objects
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec17

    from Bio.Seq import Seq
    my_seq = Seq("AGTACACTGGT")

However, Biopython doesn't know if your sequence is DNA. Seq objects can be DNA, RNA or protein. We can use the IUPAC standards to define what kind of sequence this is.

    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)

Now that Biopython knows we are dealing with DNA, you can use some handy functions:  
.complement()  
.reverse_complement()  
.transcribe()  
...and many more.  

Find the reverse complement of this sequence:  
ACAAGATGCCATTGTCCCCCGGCCTCCTG  

###SeqRecord objects
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec32

The SeqRecord class is a more sophisticated way to represent a sequence. It has a number of attributes, but the most useful are usually:  
.seq  
– The sequence itself, typically a Seq object.  
.id  
– The primary ID used to identify the sequence – a string. In most cases this is something like an accession number.  

You can create a SeqRecord object by hand like this:  
  
    from Bio.Seq import Seq
    simple_seq = Seq("GATC")
    from Bio.SeqRecord import SeqRecord
    simple_seq_r = SeqRecord(simple_seq)
    simple_seq_r.id

**Try this**  

However usually, you would obtain a SeqRecord object by reading in a file, such as fasta file.
http://en.wikipedia.org/wiki/FASTA_format

Here is an example DNA sequence in FASTA format:

\>AB000263 |acc=AB000263|descr=Homo sapiens mRNA for prepro cortistatin like peptide, complete cds.|len=368
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
