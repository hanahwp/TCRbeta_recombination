import time
t0=time.time()
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
from data_arrange import *
from recombination import *

extract_data("rawdata/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP.txt")
l_match("rawdata/lregion1.txt")
#lead seq is included in TRBV_F
TRBLV=[genedict(list(SeqIO.parse("extdata/TRBLV_F.fasta", "fasta"))).next()]
TRBD=[genedict(list(SeqIO.parse("extdata/TRBD1_F.fasta", "fasta"))).next(),genedict(list(SeqIO.parse("extdata/TRBD2_F.fasta", "fasta"))).next()]
TRBJ=[genedict(list(SeqIO.parse("extdata/TRBJ1_F.fasta", "fasta"))).next(),genedict(list(SeqIO.parse("extdata/TRBJ2_F.fasta", "fasta"))).next()]
TRBC=[genedict(list(SeqIO.parse("extdata/TRBC1_F.fasta", "fasta"))).next(),genedict(list(SeqIO.parse("extdata/TRBC2_F.fasta", "fasta"))).next()]
t1=time.time()
print "import, arrangement time: ", t1-t0

tcrb = dict()
nsense = dict()
frameout = dict()
loop_no = 0
tcrb_no = 0
nsense_no = 0
frameout_no = 0
while len(tcrb) <20:
    t2=time.time()
    #returns no., name, seq, aa
    recomb = recombinate(TRBLV, TRBD, TRBJ, TRBC)
    picked = checkseq(recomb, tcrb, nsense, frameout)
    if picked == "tcrb":
        tcrb_no += 1
    elif picked == "frameout":
        frameout_no += 1
    elif picked == "nsense":
        nsense_no += 1
    loop_no +=1
print "Total:", loop_no, "\nUnique:", len(tcrb)+len(frameout)+len(nsense)
print "Productive_Total:", tcrb_no, "\nProductive_Unique:",  len(tcrb)
print "Out_of_Frame_Total:", frameout_no, "\nOut_of_Frame_Unique:", len(frameout)
print "Has_Stop_Total:", nsense_no, "\nHas_Stop_Unique:", len(nsense)

writetcrb= [SeqRecord(Seq(str(tcrb[item]), SingleLetterAlphabet()),id=item, dbxrefs=[]) for item in tcrb]
SeqIO.write(writetcrb, "extdata/trb.fasta", "fasta")

file = open("extdata/explanation.txt", "w")
file.write("1.'A00000' is IMGT/LIGM-DB accession number. Four accession numbers of each V, D, J, and C genes are serially aligned\n")
file.write("2. gene and allele names are serially written with semicolon in between\n")
file.write("3. combined sequence")
file.close()
t4=time.time()

print 'run time :', t4-t0
