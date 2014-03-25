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

**Try this**  
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

Bio.SeqIO.parse() takes a file handle (or filename) and format string, and returns a SeqRecord iterator.

    from Bio import SeqIO
    for record in SeqIO.parse("example.fasta", "fasta") :
        print record.id

For writing records to a file use the function Bio.SeqIO.write(), which takes a SeqRecord iterator (or list of SeqRecords), output handle (or filename) and a format string:

    from Bio import SeqIO
    list_of_SeqRecords = [record1, record2]
    SeqIO.write(list_of_SeqRecords, "example.fasta", "fasta")

**Try these**  

resistance_sample.fasta contains 10 genes known to cause resistance to anibiotics in bacteria  
Write a script that prints the id of every sequence in this file. Then edit it so that it also prints the GC content of each sequence (the proportion of bases that are G or C).  
Note: you can treat a Seq object like a string and maniplate it in the ususal ways.  

resistance.fasta contains the full set of bacterial resistance genes.  
Write a script that reads this file, then writes a new fasta file containing only those sequences with "blaOXA" in the file name (there should be 202).


####Want more Bioinformatics problems? Try these resources:
* http://rosalind.info/
* http://biopython.org/DIST/docs/tutorial/Tutorial.html
And just in case you ever need to read a BAM file:  
* http://pysam.readthedocs.org/en/latest/
