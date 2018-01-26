## please put a comment
## this program is...
## input file is..
## output file is ...
## date 2018. 1. 26
## by HWPark at CSBL, Korea University

#|F| |ORF| |(F)| |P| |(P)| |[ORF]| |[F]| |[P]|

import time
t0=time.time()
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import numpy as np
import random
from random import choice
import os.path
from collections import defaultdict as ddict
directory = [os.getcwd(), '/extdata']
final_directory="".join(directory)
if not os.path.exists(final_directory):
    os.makedirs(final_directory)
rawdata = list(SeqIO.parse("rawdata/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP.txt", "fasta"))
homo=[]
for element in rawdata:
    if 'Homo' in element.id:
        homo.append(element)
TRB =[]
for element in homo:
    if 'TRB' in element.id:
        TRB.append(element)
TRBV=[]
TRBD1=[]
TRBD2=[]
TRBJ1=[]
TRBJ2=[]
TRBC1=[]
TRBC2=[]
for element in TRB:
    if 'TRBV' in element.id:
        TRBV.append(element)
    if 'TRBD1' in element.id:
        TRBD1.append(element)
    if 'TRBD2' in element.id:
        TRBD2.append(element)
    if 'TRBJ1' in element.id:
        TRBJ1.append(element)
    if 'TRBJ2' in element.id:
        TRBJ2.append(element)
    if 'TRBC1' in element.id:
        TRBC1.append(element)
    if 'TRBC2' in element.id:
        TRBC2.append(element)
    #V
TRBV_F = []
TRBV_P = []
for i in range(len(TRBV)):
    if any(c in TRBV[i].description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
        TRBV_F.append(TRBV[i])
    if any(c in TRBV[i].description for c in ('|P|', '|(P)|','|[P]|')) :
        TRBV_P.append(TRBV[i])
#print "# of V : ", len(TRBV), "\n# of fuctional V : ", len(TRBV_F), "\n# of pseudo V : ", len(TRBV_P)
SeqIO.write(TRBV_F, "extdata/TRBV_F.fasta", "fasta")

    #D
TRBD1_F = []
TRBD2_F = []
TRBD1_P = []
TRBD2_P = []
for element in TRBD1:
    if any(c in element.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
        TRBD1_F.append(element)
    if any(c in element.description for c in ('|P|', '|(P)|','|[P]|')) :
        TRBD1_P.append(element)
for element in TRBD2:
    if any(c in element.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
        TRBD2_F.append(element)
    if any(c in element.description for c in ('|P|', '|(P)|','|[P]|')) :
        TRBD2_P.append(element)
#print "\n# of D : ", len(TRBD1)+len(TRBD2), "\n# of D1 :",len(TRBD1), "\n# of D2 :",len(TRBD2),"\n# of fuctional D : ", len(TRBD1_F)+len(TRBD2_F)
#print "# of pseudo D : ", len(TRBD1_P)+len(TRBD2_P), "\n# of fuctional TRBD1 : ", len(TRBD1_F), "\n# of fuctional TRBD2 : ", len(TRBD2_F)
SeqIO.write(TRBD1_F, "extdata/TRBD1_F.fasta", "fasta")
SeqIO.write(TRBD2_F, "extdata/TRBD2_F.fasta", "fasta")

    #J
TRBJ1_F = []
TRBJ2_F = []
TRBJ1_P = []
TRBJ2_P = []
for element in TRBJ1:
    if any(c in element.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
        TRBJ1_F.append(element)
    if any(c in element.description for c in ('|P|', '|(P)|','|[P]|')) :
        TRBJ1_P.append(element)
for element in TRBJ2:
    if any(c in element.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
        TRBJ2_F.append(element)
    if any(c in element.description for c in ('|P|', '|(P)|','|[P]|')) :
        TRBJ2_P.append(element)
#print "\n# of J : ", len(TRBJ1)+len(TRBJ2), "\n# of J1 :",len(TRBJ1), "\n# of J2 :",len(TRBJ2),"\n# of fuctional J : ", len(TRBJ1_F)+len(TRBJ2_F)
#print "# of pseudo J : ", len(TRBJ1_P)+len(TRBJ2_P), "\n# of fuctional TRBJ1 : ", len(TRBJ1_F), "\n# of fuctional TRBJ2 : ", len(TRBJ2_F)
SeqIO.write(TRBJ1_F, "extdata/TRBJ1_F.fasta", "fasta")
SeqIO.write(TRBJ2_F, "extdata/TRBJ2_F.fasta", "fasta")


    #C
TRBC1_F = []
TRBC2_F = []
TRBC1_P = []
TRBC2_P = []

for element in TRBC1:
    if any(c in element.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
        TRBC1_F.append(element)
    if any(c in element.description for c in ('|P|', '|(P)|','|[P]|')) :
        TRBC1_P.append(element)
for element in TRBC2:
    if any(c in element.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
        TRBC2_F.append(element)
    if any(c in element.description for c in ('|P|', '|(P)|','|[P]|')) :
        TRBC2_P.append(element)
#print "\n# of C : ", len(TRBC1)+len(TRBC2), "\n# of C1 :",len(TRBC1), "\n# of C2 :",len(TRBC2),"\n# of fuctional C : ", len(TRBC1_F)+len(TRBC2_F)
#print "# of pseudo C : ", len(TRBC1_P)+len(TRBC2_P), "\n# of fuctional TRBC1 : ", len(TRBC1_F), "\n# of fuctional TRBC2 : ", len(TRBC2_F)
SeqIO.write(TRBC1_F, "extdata/TRBC1_F.fasta", "fasta")
SeqIO.write(TRBC2_F, "extdata/TRBC2_F.fasta", "fasta")


#appending exons of C
f1= SeqRecord(Seq(str(TRBC1_F[0].seq+TRBC1_F[1].seq+TRBC1_F[2].seq+TRBC1_F[3].seq), SingleLetterAlphabet()),id=TRBC1_F[0].id, name=TRBC1_F[0].name, dbxrefs=TRBC1_F[0].dbxrefs)
f2= SeqRecord(Seq(str(TRBC1_F[4].seq+TRBC1_F[5].seq+TRBC1_F[6].seq+TRBC1_F[7].seq), SingleLetterAlphabet()),id=TRBC1_F[4].id, name=TRBC1_F[4].name, dbxrefs=TRBC1_F[4].dbxrefs)
TRBC1_F =[f1,f2]
f1= SeqRecord(Seq(str(TRBC2_F[0].seq+TRBC2_F[1].seq+TRBC2_F[2].seq+TRBC2_F[3].seq), SingleLetterAlphabet()),id=TRBC2_F[0].id, name=TRBC2_F[0].name, dbxrefs=TRBC2_F[0].dbxrefs)
f2= SeqRecord(Seq(str(TRBC2_F[4].seq+TRBC2_F[5].seq+TRBC2_F[6].seq+TRBC2_F[7].seq), SingleLetterAlphabet()),id=TRBC2_F[4].id, name=TRBC2_F[4].name, dbxrefs=TRBC2_F[4].dbxrefs)
TRBC2_F =[f1,f2]

def todict(arrlist):    
#todict(TRBJ1_F) -> {'TRBJ1-1':[*1,*2], 'TRBJ1-2':[*1,*2, *3]}
    a = {}
    for j in arrlist:
        a1 =[]
        name = j.id[j.id.find('TRB'):j.id.find('*')]
        for i in arrlist:
            match = i.id[i.id.find('TRB'):i.id.find('*')]
            if name == match:
                a1.append(i)
                a[name] = a1
    yield a

TRBCdict = [todict(TRBC1_F).next(), todict(TRBC2_F).next()]
TRBJdict = [todict(TRBJ1_F).next(), todict(TRBJ2_F).next()]
TRBDdict = [todict(TRBD1_F).next(), todict(TRBD2_F).next()]
TRBVdict = [todict(TRBV_F).next()]


trb = dict()
degrade = dict()
loopno = 0
nsensepick = 0
sensepick = 0
t1=time.time()
print "import, arrangement time: ", t1-t0
while len(trb) <900:
    t2=time.time()
    loopno +=1
    #choose 1or2
    n = np.random.randint(0,2)
    D = TRBDdict[n] 
    J = TRBJdict[n]
    C = TRBCdict[n]
    V = TRBVdict[0]
    
    #choose gene
    vgene = V[random.choice(V.keys())]
    dgene = D[random.choice(D.keys())]
    jgene = J[random.choice(J.keys())]
    cgene = C[random.choice(C.keys())]

    #choose allele
    call=  choice(cgene)
    dall = choice(dgene)
    jall = choice(jgene)
    vall = choice(vgene)

    #get accesssion no., allele name, and seq
    creg  = [call.id[:6], call.id[call.id.find('TRB'):call.id.find('*')+3], call.seq]
    dreg  = [dall.id[:6], dall.id[dall.id.find('TRB'):dall.id.find('*')+3], dall.seq]
    jreg  = [jall.id[:6], jall.id[jall.id.find('TRB'):jall.id.find('*')+3], jall.seq]
    vreg  = [vall.id[:6], vall.id[vall.id.find('TRB'):vall.id.find('*')+3], vall.seq]

    #recombination
    recomb_assno= vreg[0]+dreg[0]+jreg[0]+creg[0]
    recomb_name= vreg[1]+";"+dreg[1]+";"+jreg[1]+";"+creg[1] 
    recomb_seq= vreg[2]+dreg[2]+jreg[2]+creg[2]
    
    # translate
    aa = recomb_seq.translate()
    #assume that shifted frame(not using the conventional stop codon) is degraded
    if '*' in aa:
        nsensepick += 1
        degrade["%s|%s"%(recomb_assno, recomb_name)] = recomb_seq
    else:
        if (len(recomb_seq) % 3) == 0:
            trb["%s|%s"%(recomb_assno, recomb_name)] = recomb_seq
            sensepick += 1
        else:
            nsensepick += 1
            degrade["%s|%s"%(recomb_assno, recomb_name)] = recomb_seq
    t3=time.time()
    
#    print "loop %d"%loopno, t3-t2

trblist= [SeqRecord(Seq(str(trb[item]), SingleLetterAlphabet()),id=item, dbxrefs=[]) for item in trb]
SeqIO.write(trblist, "extdata/trb.fasta", "fasta")

file = open("extdata/explanation.txt", "w")
file.write("1.'A00000' is IMGT/LIGM-DB accession number. Four accession numbers of each V, D, J, and C genes are serially aligned\n")
file.write("2. gene and allele names are serially written with semicolon in between\n")
file.write("3. combined sequence")
file.close()
t4=time.time()

print 'total_pick#:', sensepick+nsensepick, ',total_trb#:', sensepick, ',unique_trb#:', len(trb)
print "total_trb/total_pick%:", (100*sensepick/loopno),
print 'unique_tcr/unique_pick%:',100*len(trb)/(len(trb)+len(degrade))
print 'run time :', t4-t0

