from Bio import SeqIO

infilename = "resistance.fasta"
outfilename = "resistance.blaOXA.fasta"

blaOXA_records = list()

count = 0
for seq_record in SeqIO.parse(infilename, "fasta"):
    #print seq_record.id
    #print seq_record.seq
    if "blaOXA" in seq_record.id:
        blaOXA_records.append(seq_record)
        count += 1

SeqIO.write(blaOXA_records, outfilename, "fasta")
print count
