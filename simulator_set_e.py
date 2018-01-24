# -*- coding: utf-8 -*-
import time
t0 = time.time()
from Bio import SeqIO
import numpy as np
from random import choice
from Bio.Alphabet import IUPAC,  generic_dna
import random
from Bio.Seq import Seq

TRBC1_F = list(SeqIO.parse("TRBC1_F.fasta", "fasta"))
TRBC2_F = list(SeqIO.parse("TRBC2_F.fasta", "fasta"))
TRBD1_F = list(SeqIO.parse("TRBD1_F.fasta", "fasta"))
TRBD2_F = list(SeqIO.parse("TRBD2_F.fasta", "fasta"))
TRBJ1_F = list(SeqIO.parse("TRBJ1_F.fasta", "fasta"))
TRBJ2_F = list(SeqIO.parse("TRBJ2_F.fasta", "fasta"))
TRBV_F = list(SeqIO.parse("TRBV_F.fasta", "fasta"))
def arrasdic_gene(seqlist):    
#arrange('TRBC1_F') -> TRBC1:[,], TRBC2:[,]
    a = {}
    for j in range(len(seqlist)):
        a1 =[]
        name = seqlist[j].id[seqlist[j].id.find('TRB'):seqlist[j].id.find('*')]
        for i in range(len(seqlist)):
            match = seqlist[i].id[seqlist[i].id.find('TRB'):seqlist[i].id.find('*')]
            if name == match:
                a1.append(seqlist[i])
                a[name] = a1
    yield a

TRBC = [arrasdic_gene(TRBC1_F).next(), arrasdic_gene(TRBC2_F).next()]
TRBJ = [arrasdic_gene(TRBJ1_F).next(), arrasdic_gene(TRBJ2_F).next()]
TRBD = [arrasdic_gene(TRBD1_F).next(), arrasdic_gene(TRBD2_F).next()]
TRBV = [arrasdic_gene(TRBV_F).next()]
tcrs = set()
degrade = set()


t1 = time.time()
loopno = 0
nsensepick = 0
sensepick = 0
while len(tcrs) <200:
    loopno +=1
    t2 = time.time()
    #choose 1or2
    n = np.random.randint(0,2)
    D = TRBD[n] 
    J = TRBJ[n]
    C = TRBC[n]
    V = TRBV[0]
    #choose gene
    vgene = V[random.choice(V.keys())]
    dgene = D[random.choice(D.keys())]
    jgene = J[random.choice(J.keys())]
    cgene = C[random.choice(C.keys())]

    #choose allele
    c_exon_intgr = [cgene[0].seq+cgene[1].seq+cgene[2].seq+cgene[3].seq, cgene[4].seq+cgene[5].seq+cgene[6].seq+cgene[7].seq]
    call=  choice(c_exon_intgr)
    dall = choice(dgene).seq
    jall = choice(jgene).seq
    vall = choice(vgene).seq

    #recombination
    recomb= vall+dall+jall+call

    # translate
    aa = recomb.translate()
    #assume that shifted frame(not using the conventional stop codon) is degraded
    if '*' in aa:
        nsensepick += 1
        degrade.update([recomb])
    else:
        if (len(recomb) % 3) == 0:
            sensepick += 1
            tcrs.update([recomb])
        else:
            nsensepick += 1
            degrade.update([recomb])
    t3 = time.time()
#    print "loop %d" % loopno, t3 - t2
        
file = open("tcrs.txt", "w")
for item in tcrs:
    file.write("%s\n" % item)
file.close()

#file = open("degrade.txt", "w")
#for item in degrade:
#    file.write("%s\n" % item)
#file.close()

print 'total_pick#:', sensepick+nsensepick, ',total_tcrs#:', sensepick, ',unique_tcrs#:', len(tcrs)
print "total_tcrs/total_pick%:", (100*sensepick/loopno),
print 'unique_pick/total_pick%:', (100*(len(tcrs)+len(degrade))/loopno), 'unique_tcr/unique_pick%:',100*len(tcrs)/(len(tcrs)+len(degrade))
t4=time.time()
print t4-t0



