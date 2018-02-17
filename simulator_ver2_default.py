
## this program generates translatable rearranged trb genes
## input file is "IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP.txt", downloaded from IMGT
## output file is [explanation.txt], tcrb.fasta
## date 2018. 2. 17.
## by HWPark at CSBL, Korea University

import time
t0=time.time()
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import numpy 
from gene_arrange import *
from recombination import *

wanted_no =3000 #input("Type in the number of unique tcrb you want: ")
v_insmax =12 #input("Type in the number of maximum insertion for VD junction: ")
j_insmax = 12 #input("Type in the number of maximum insertion for DJ junction: ")
vdel = 12#input("Type in the number of maximum deletion for V 3': ")
jdel = 12#input("Type in the number of maximum deletion for 5 3': ")



extract_data("rawdata/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP.txt")
l_match("rawdata/LVEXON.txt")
#lead seq is included in TRBV_F
TRBLV=[genedict(list(SeqIO.parse("extdata/TRBLV_F.fasta", "fasta"))).next()]
TRBD=[genedict(list(SeqIO.parse("extdata/TRBD1_F.fasta", "fasta"))).next(),genedict(list(SeqIO.parse("extdata/TRBD2_F.fasta", "fasta"))).next()]
TRBJ=[genedict(list(SeqIO.parse("extdata/TRBJ1_F.fasta", "fasta"))).next(),genedict(list(SeqIO.parse("extdata/TRBJ2_F.fasta", "fasta"))).next()]
TRBC=[genedict(list(SeqIO.parse("extdata/TRBC1_F.fasta", "fasta"))).next(),genedict(list(SeqIO.parse("extdata/TRBC2_F.fasta", "fasta"))).next()]
cys2 = list(SeqIO.parse("rawdata/cys2.txt", "fasta"))
jphe =list(SeqIO.parse("rawdata/jphe.txt",  "fasta"))

    
tcrb = dict()
nsense = dict()
frameout = dict()
unconserved = dict()
loop_no = 0
tcrb_no = 0
nsense_no = 0
frameout_no = 0
unconserved_no = 0

while len(tcrb) < wanted_no:
    #returns no., name, seq, aa
    rearranged = recomb_n_check(TRBLV, TRBD, TRBJ, TRBC, v_insmax, j_insmax, vdel, jdel, cys2,jphe )
    if rearranged[1] == "tcrb":
        keynamevar_dict("%s|%s"%(rearranged[0]["assno"], rearranged[0]["name"]),str(rearranged[0]["seq"]),tcrb )
        tcrb_no += 1
    elif rearranged[1] == "frameout":
        keynamevar_dict("%s|%s"%(rearranged[0]["assno"], rearranged[0]["name"]),str(rearranged[0]["seq"]),frameout )
        frameout_no += 1
    elif rearranged[1]  == "nsense":
        keynamevar_dict("%s|%s"%(rearranged[0]["assno"], rearranged[0]["name"]),str(rearranged[0]["seq"]),nsense )
        nsense_no += 1
    elif rearranged[1] == "unconserved":
        keynamevar_dict("%s|%s"%(rearranged[0]["assno"], rearranged[0]["name"]),str(rearranged[0]["seq"]),unconserved )
        unconserved_no+=1
    loop_no +=1

print "Total:", loop_no, "\nUnique:", len(tcrb)+len(frameout)+len(nsense)
print "Productive_Total:", tcrb_no, "\nProductive_Unique:",  len(tcrb)
print "Out_of_Frame_Total:", frameout_no, "\nOut_of_Frame_Unique:", len(frameout)
print "Has_Stop_Total:", nsense_no, "\nHas_Stop_Unique:", len(nsense)
print "Unconserved_Total:", unconserved_no, "\nUnconserved_Unique:", len(unconserved)
print "Productive_Total/Total % : ", numpy.around((100*float(tcrb_no))/float(loop_no),decimals=2)

writetcrb= [SeqRecord(Seq(str(tcrb[item][0]), SingleLetterAlphabet()),id=item, dbxrefs=[]) for item in tcrb]
SeqIO.write(writetcrb, "extdata/tcrb.fasta", "fasta")

explanation = open("extdata/explanation_of_trb.fasta.txt", "w")
readme = ("1.'A00000' is IMGT/LIGM-DB accession number. Four accession numbers of each V, D, J, and C genes are serially aligned\n",
          "2. gene and allele names are serially written with semicolon in between\n",
          "3. number after ';;;' shows the ordinal number of variations generated from the same combination of genes\n",
          "4. combined sequence\n",
          "#by HWPark at CSBL, Korea University")
explanation.write(''.join(readme))
explanation.close()

nonsensefile = [SeqRecord(Seq(str(nsense[item][0]), SingleLetterAlphabet()),id=item, dbxrefs=[]) for item in nsense]
SeqIO.write(nonsensefile, "extdata/nonsensefile.fasta", "fasta")

t4=time.time()

print 'run time :', t4-t0
