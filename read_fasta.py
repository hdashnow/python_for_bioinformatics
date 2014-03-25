from Bio import SeqIO

infilename = "resistance_sample.fasta"
outfilename = "resistance_sample.blaOXA.fasta"

blaOXA_records = []

for seq_record in SeqIO.parse(infilename, "fasta"):
    #print seq_record.id
    #print seq_record.seq
    if "blaOXA" in seq_record.id:
        blaOXA_records.append(seq_record)

SeqIO.write(blaOXA_records, outfilename, "fasta")
